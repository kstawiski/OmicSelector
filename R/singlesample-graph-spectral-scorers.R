#' @title Fit a graph-Fourier discriminator on a frozen co-expression graph
#'
#' @description
#' Learns a single-sample circulating-miRNA discriminator by building a frozen
#' feature-feature co-expression graph from the training set, extracting the
#' low-frequency graph-Fourier modes of its symmetric-normalized Laplacian, and
#' fitting a frozen linear head on per-specimen graph-spectral coordinates. Each
#' specimen is represented by a robust centered log-ratio (rCLR) transform
#' computed from that specimen alone: positive entries are log-transformed and
#' centered by the mean log abundance of the specimen's own positive entries,
#' while zeros remain neutral at 0.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{cor_method}, one of
#'   \code{"spearman"} or \code{"pearson"} (default \code{"spearman"});
#'   \code{graph_k}, positive integer number of nearest neighbours retained per
#'   feature before symmetrization (default \code{min(10L, D - 1L)});
#'   \code{n_modes}, positive integer cap on retained low-frequency graph
#'   Fourier modes (default \code{min(10L, D)}); \code{min_features}, positive
#'   integer scoring overlap floor (default \code{3L}); and \code{eps}, positive
#'   finite numerical floor (default \code{1e-8}).
#'
#' @return A plain list of class \code{gsp_gft_model} containing the frozen
#'   \code{feature_universe}, low-frequency basis \code{U_low}, retained
#'   eigenvalues \code{eigvals_low}, graph settings, frozen linear head, and the
#'   resolved hyperparameters.
#'
#' @details
#' The co-expression graph is learned from training rCLR signals only. The
#' absolute training feature-feature correlation matrix is kNN-sparsified with a
#' deterministic feature-index tie break, symmetrized by union, and converted to
#' the symmetric-normalized graph Laplacian
#' \deqn{L = I - D^{-1/2} W D^{-1/2}.}
#' Its eigenvectors, ordered by increasing eigenvalue, define graph-Fourier
#' coordinates. The retained low-frequency coordinates are standardized by
#' frozen training means and scales and then passed to a frozen linear head. A
#' binomial GLM is used when it is finite and converged; otherwise the fit falls
#' back to a diagonal-LDA mean-difference head in the same standardized
#' coordinates. Larger scores are oriented to be more case-like.
#'
#' At inference, no statistic is estimated from the scoring matrix. The graph,
#' basis, coordinate center/scale, and head coefficients are all frozen at fit
#' time from the training data. Each row's rCLR vector uses only that row's own
#' positive present features, with absent and zero features encoded as 0.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(60 * 30, shape = 3), nrow = 60,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 30)
#' X[y == 1, 1:6] <- X[y == 1, 1:6] + 20
#' model <- fit_gsp_gft(X, y)
#' score <- score_gsp_gft(model, X)
#' }
#'
#' @references
#' Shuman DI, Narang SK, Frossard P, Ortega A, Vandergheynst P. (2013) The
#' emerging field of signal processing on graphs. \emph{IEEE Signal Processing
#' Magazine} 30: 83-98.
#'
#' Chung FRK. (1997) \emph{Spectral Graph Theory}. American Mathematical
#' Society.
#'
#' @export
fit_gsp_gft <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_gsp_gft", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_gsp_gft", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_gsp_gft")
  hp <- .gsp_gft_resolve_hp(hp, ncol(X_train))
  if (ncol(X_train) < hp$min_features) {
    stop("fit_gsp_gft: X_train must contain at least hp$min_features features")
  }

  feature_universe <- colnames(X_train)
  Z <- .gsp_gft_rclr_matrix(X_train)
  graph <- .gsp_gft_graph_basis(Z, hp)
  C <- Z %*% graph$U_low
  head <- .gsp_gft_fit_head(C, y, hp)

  model <- list(
    feature_universe = feature_universe,
    U_low = graph$U_low,
    eigvals_low = graph$eigvals_low,
    graph = list(cor_method = hp$cor_method, graph_k = hp$graph_k),
    head = head,
    hp = hp
  )
  class(model) <- "gsp_gft_model"
  model
}


#' @title Score specimens with a frozen graph-Fourier discriminator
#'
#' @description
#' Scores each row of \code{X} independently with a frozen
#' \code{gsp_gft_model}. For each specimen, the scorer computes a per-sample
#' rCLR vector over the model features that are present in \code{X}, aligns that
#' vector to the frozen training feature universe with absent features set to 0,
#' projects it onto the frozen low-frequency graph-Fourier basis, applies the
#' frozen coordinate standardization, and evaluates the frozen linear head.
#'
#' @param model A \code{gsp_gft_model} object returned by
#'   \code{\link{fit_gsp_gft}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Larger values
#'   are more case-like. If fewer than \code{model$hp$min_features} model
#'   features are present, the documented neutral score \code{0} is returned for
#'   every row.
#'
#' @details
#' The score of row \eqn{i} depends only on \code{X[i, ]} and the frozen model.
#' The scored batch contributes no centering, normalization, graph update,
#' bandwidth, basis, or head statistic. The rCLR representation is exactly
#' invariant to multiplying one specimen by a positive scalar, because
#' \eqn{\log(c v_j) - mean_k \log(c v_k) = \log(v_j) - mean_k \log(v_k)}
#' over that specimen's own positive entries.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(60 * 30, shape = 3), nrow = 60,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 30)
#' X[y == 1, 1:6] <- X[y == 1, 1:6] + 20
#' model <- fit_gsp_gft(X, y)
#' score_gsp_gft(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Shuman DI, Narang SK, Frossard P, Ortega A, Vandergheynst P. (2013) The
#' emerging field of signal processing on graphs. \emph{IEEE Signal Processing
#' Magazine} 30: 83-98.
#'
#' Chung FRK. (1997) \emph{Spectral Graph Theory}. American Mathematical
#' Society.
#'
#' @export
score_gsp_gft <- function(model, X, meta = NULL) {
  if (!inherits(model, "gsp_gft_model")) {
    stop("score_gsp_gft: model must have class gsp_gft_model")
  }
  X <- .reo_check_matrix(X, "score_gsp_gft", "X")
  .reo_check_meta(meta, nrow(X), "score_gsp_gft", "meta")

  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  Z_align <- .gsp_gft_aligned_rclr_matrix(
    X = X,
    feature_universe = model$feature_universe,
    present = present
  )
  out <- vapply(seq_len(nrow(Z_align)), function(i) {
    c_i <- as.numeric(Z_align[i, , drop = FALSE] %*% model$U_low)
    c_std <- (c_i - model$head$center) / model$head$scale
    sum(c_std * model$head$coef) + model$head$intercept
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_gsp_gft: scorer produced non-finite or wrong-length output")
  }
  out
}

.gsp_gft_resolve_hp <- function(hp, p) {
  if (!is.list(hp)) stop("fit_gsp_gft: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(is.na(nm)) || any(!nzchar(nm))) {
      stop("fit_gsp_gft: hp fields must be named")
    }
  }
  allowed <- c("cor_method", "graph_k", "n_modes", "min_features", "eps")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_gsp_gft: unknown hp field(s): ", paste(unknown, collapse = ", "))
  }

  cor_method <- hp$cor_method
  if (is.null(cor_method)) cor_method <- "spearman"
  if (!is.character(cor_method) || length(cor_method) != 1L ||
      is.na(cor_method) || !cor_method %in% c("pearson", "spearman")) {
    stop("fit_gsp_gft: hp$cor_method must be 'pearson' or 'spearman'")
  }

  graph_k <- hp$graph_k
  if (is.null(graph_k)) graph_k <- min(10L, max(1L, p - 1L))
  graph_k <- .gsp_gft_check_positive_integer(
    graph_k, "graph_k", "fit_gsp_gft"
  )
  if (graph_k >= p) {
    stop("fit_gsp_gft: hp$graph_k must be less than ncol(X_train)")
  }

  n_modes <- hp$n_modes
  if (is.null(n_modes)) n_modes <- min(10L, p)
  n_modes <- .gsp_gft_check_positive_integer(
    n_modes, "n_modes", "fit_gsp_gft"
  )

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  min_features <- .gsp_gft_check_positive_integer(
    min_features, "min_features", "fit_gsp_gft"
  )

  eps <- hp$eps
  if (is.null(eps)) eps <- 1e-8
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_gsp_gft: hp$eps must be a positive finite number")
  }

  list(
    cor_method = cor_method,
    graph_k = graph_k,
    n_modes = n_modes,
    min_features = min_features,
    eps = as.numeric(eps)
  )
}

.gsp_gft_check_positive_integer <- function(x, field, fn) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1L ||
      x > .Machine$integer.max || x != floor(x)) {
    stop(fn, ": hp$", field, " must be a positive integer")
  }
  as.integer(x)
}

.gsp_gft_rclr_row <- function(v) {
  z <- numeric(length(v))
  pos <- which(v > 0)
  if (length(pos) > 0L) {
    g <- mean(log(v[pos]))
    z[pos] <- log(v[pos]) - g
  }
  z
}

.gsp_gft_rclr_matrix <- function(X) {
  Z <- do.call(rbind, lapply(seq_len(nrow(X)), function(i) {
    .gsp_gft_rclr_row(X[i, ])
  }))
  dimnames(Z) <- dimnames(X)
  storage.mode(Z) <- "double"
  Z
}

.gsp_gft_aligned_rclr_matrix <- function(X, feature_universe, present) {
  X_present <- X[, present, drop = FALSE]
  Z_present <- .gsp_gft_rclr_matrix(X_present)
  Z_align <- matrix(
    0,
    nrow = nrow(X),
    ncol = length(feature_universe),
    dimnames = list(rownames(X), feature_universe)
  )
  Z_align[, match(present, feature_universe)] <- Z_present
  Z_align
}

.gsp_gft_graph_basis <- function(Z, hp) {
  D <- ncol(Z)
  S <- suppressWarnings(stats::cor(Z, method = hp$cor_method))
  S <- abs(S)
  S[!is.finite(S)] <- 0
  diag(S) <- 0

  S_knn <- matrix(0, nrow = D, ncol = D, dimnames = dimnames(S))
  for (i in seq_len(D)) {
    others <- setdiff(seq_len(D), i)
    ord <- order(-S[i, others], others, method = "radix")
    keep <- others[ord[seq_len(min(hp$graph_k, length(others)))]]
    S_knn[i, keep] <- S[i, keep]
  }
  W <- pmax(S_knn, t(S_knn))
  diag(W) <- 0

  d <- rowSums(W)
  d_safe <- pmax(d, hp$eps)
  Dm12 <- 1 / sqrt(d_safe)
  L <- diag(D) - sweep(sweep(W, 1L, Dm12, "*"), 2L, Dm12, "*")
  L <- (L + t(L)) / 2

  e <- eigen(L, symmetric = TRUE)
  idx <- order(e$values, seq_along(e$values), method = "radix")
  K <- min(hp$n_modes, D)
  U_low <- e$vectors[, idx[seq_len(K)], drop = FALSE]
  U_low <- .gsp_gft_fix_eigen_signs(U_low)
  rownames(U_low) <- colnames(Z)
  colnames(U_low) <- paste0("gft", seq_len(K))

  list(
    U_low = U_low,
    eigvals_low = as.numeric(e$values[idx[seq_len(K)]])
  )
}

.gsp_gft_fix_eigen_signs <- function(U) {
  for (j in seq_len(ncol(U))) {
    idx <- which.max(abs(U[, j]))
    if (is.finite(U[idx, j]) && U[idx, j] < 0) {
      U[, j] <- -U[, j]
    }
  }
  U
}

.gsp_gft_fit_head <- function(C, y, hp) {
  center <- colMeans(C)
  scale <- vapply(seq_len(ncol(C)), function(j) stats::sd(C[, j]), numeric(1))
  scale[!is.finite(scale) | scale < hp$eps] <- hp$eps
  C_std <- sweep(sweep(C, 2L, center, "-"), 2L, scale, "/")
  colnames(C_std) <- paste0("gft", seq_len(ncol(C_std)))

  warn <- character()
  dat <- data.frame(y = y, C_std, check.names = FALSE)
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm(y ~ ., data = dat, family = stats::binomial()),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )
  bad_glm <- is.null(fit) ||
    !isTRUE(fit$converged) ||
    any(!is.finite(stats::coef(fit))) ||
    any(grepl("algorithm did not converge|fitted probabilities numerically",
              warn))

  if (!bad_glm) {
    beta <- stats::coef(fit)
    return(list(
      type = "glm",
      coef = as.numeric(beta[-1]),
      intercept = as.numeric(beta[1]),
      center = as.numeric(center),
      scale = as.numeric(scale)
    ))
  }

  lda <- .gsp_gft_lda_head(C_std, y, hp$eps)
  list(
    type = "diag_lda",
    coef = as.numeric(lda$coef),
    intercept = as.numeric(lda$intercept),
    center = as.numeric(center),
    scale = as.numeric(scale)
  )
}

.gsp_gft_lda_head <- function(C_std, y, eps) {
  case <- C_std[y == 1L, , drop = FALSE]
  control <- C_std[y == 0L, , drop = FALSE]
  mu_case <- colMeans(case)
  mu_control <- colMeans(control)
  var_case <- .gsp_gft_col_var(case)
  var_control <- .gsp_gft_col_var(control)
  n_case <- nrow(case)
  n_control <- nrow(control)
  denom <- max(n_case + n_control - 2L, 1L)
  pooled_var <- ((n_case - 1L) * var_case + (n_control - 1L) * var_control) /
    denom
  coef <- (mu_case - mu_control) / (pooled_var + eps)
  midpoint <- (mu_case + mu_control) / 2
  intercept <- -sum(midpoint * coef)
  list(coef = coef, intercept = intercept)
}

.gsp_gft_col_var <- function(X) {
  if (nrow(X) < 2L) return(rep(0, ncol(X)))
  as.numeric(vapply(seq_len(ncol(X)), function(j) stats::var(X[, j]),
                    numeric(1)))
}
