#' @title Fit a Nystrom diffusion-map discriminator on frozen landmarks
#'
#' @description
#' Learns a single-sample circulating-miRNA discriminator by selecting frozen
#' training landmarks, building a diffusion-map geometry on their per-sample
#' rCLR coordinates, embedding every training specimen by the Nystrom
#' out-of-sample formula, and fitting a frozen linear head. The rCLR transform is
#' computed independently for each specimen: positive entries are log
#' transformed and centered by that specimen's own mean log-positive abundance,
#' while zero entries remain at 0.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of
#'   non-negative abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{n_landmarks} positive integer
#'   landmark cap (default \code{min(200L, nrow(X_train))});
#'   \code{n_components} positive integer diffusion-coordinate cap (default
#'   \code{10L}); \code{t} positive integer diffusion time (default \code{1L});
#'   \code{gamma} positive RBF bandwidth, or \code{NULL} for the landmark
#'   median heuristic (default \code{NULL}); \code{eig_floor} non-negative
#'   eigenvalue floor for retained non-trivial modes (default \code{1e-6});
#'   \code{min_features} positive integer scoring overlap floor (default
#'   \code{3L}); \code{eps} positive numerical floor (default \code{1e-8}); and
#'   \code{seed} non-negative integer used only for deterministic landmark
#'   subsampling while restoring the global RNG state (default \code{1L}).
#'
#' @return A plain list of class \code{man_nystrom_model} containing the frozen
#'   training \code{feature_universe}, rCLR \code{landmarks},
#'   \code{landmark_idx}, frozen \code{gamma}, diffusion settings, cached
#'   full-universe \code{geometry}, frozen linear \code{head}, and resolved
#'   \code{hp}.
#'
#' @details
#' The landmark graph is learned once from training data. For landmarks
#' \eqn{L_j}, the RBF kernel \eqn{W_{ij}=\exp(-\gamma\|L_i-L_j\|^2)} is
#' symmetrically degree-normalized and eigendecomposed. The trivial top
#' diffusion pair is dropped; retained right eigenvectors of the row-stochastic
#' transition matrix define the diffusion coordinates. A new specimen \eqn{z}
#' is embedded from only its own affinities to the frozen landmarks:
#' \deqn{\Psi_k(z)=\lambda_k^{t-1}\sum_j p_j(z)\psi_{jk},}
#' where \eqn{p_j(z)} is the row-normalized RBF affinity of \eqn{z} to landmark
#' \eqn{j}. These coordinates are standardized by frozen training means/scales
#' and passed to a frozen binomial-GLM head, with a diagonal-LDA
#' mean-difference fallback if the GLM fit is numerically unsafe. Larger scores
#' are oriented to be more case-like.
#'
#' No scored-batch statistic is used. At scoring, the present feature set is
#' \code{intersect(model$feature_universe, colnames(X))}; the landmark geometry
#' is re-derived from frozen landmark coordinates restricted to those features,
#' and each scored row contributes only its own per-sample rCLR vector and its
#' own affinity row to the frozen landmarks. If fewer than
#' \code{hp$min_features} model features are present, the documented neutral
#' score \code{0} is returned for every row.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(60 * 30, shape = 3), nrow = 60,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 30)
#' X[y == 1, 1:6] <- X[y == 1, 1:6] + 20
#' X[y == 0, 7:12] <- X[y == 0, 7:12] + 20
#' model <- fit_man_nystrom(X, y)
#' score <- score_man_nystrom(model, X)
#' }
#'
#' @references
#' Coifman RR, Lafon S. (2006) Diffusion maps. \emph{Applied and Computational
#' Harmonic Analysis} 21: 5-30.
#'
#' Bengio Y, Paiement JF, Vincent P, Delalleau O, Le Roux N, Ouimet M. (2004)
#' Out-of-sample extensions for LLE, Isomap, MDS, Eigenmaps, and Spectral
#' Clustering. \emph{Advances in Neural Information Processing Systems} 16.
#'
#' @export
fit_man_nystrom <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_man_nystrom", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_man_nystrom", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_man_nystrom")
  hp <- .man_nystrom_resolve_hp(hp, nrow(X_train))
  if (ncol(X_train) < hp$min_features) {
    stop("fit_man_nystrom: X_train must contain at least hp$min_features features")
  }

  feature_universe <- colnames(X_train)
  Z <- .man_nystrom_rclr_matrix(X_train)
  landmark_idx <- .man_nystrom_landmark_indices(
    n = nrow(Z),
    n_landmarks = hp$n_landmarks,
    seed = hp$seed
  )
  L <- Z[landmark_idx, , drop = FALSE]

  gamma <- hp$gamma
  if (is.null(gamma)) {
    gamma <- .man_nystrom_median_gamma(L, hp$eps)
    hp$gamma <- gamma
  }

  geom <- .man_nystrom_geometry(L, gamma, hp)
  E_train <- .man_nystrom_embed_matrix(Z, L, gamma, geom, hp)
  head <- .man_nystrom_fit_head(E_train, y, hp)

  model <- list(
    feature_universe = feature_universe,
    landmarks = L,
    landmark_idx = landmark_idx,
    gamma = gamma,
    t = hp$t,
    eig_floor = hp$eig_floor,
    n_components = ncol(E_train),
    geometry = geom,
    head = head,
    hp = hp
  )
  class(model) <- "man_nystrom_model"
  model
}


#' @title Score specimens with a frozen Nystrom diffusion-map discriminator
#'
#' @description
#' Scores each row of \code{X} independently with a frozen
#' \code{man_nystrom_model}. For each specimen, the scorer computes a
#' per-sample rCLR vector over the model features present in \code{X}, derives a
#' diffusion-map geometry from the frozen landmarks restricted to those same
#' present features, embeds that specimen by its own Nystrom affinity row to the
#' frozen landmarks, applies frozen coordinate standardization, and evaluates
#' the frozen linear head.
#'
#' @param model A \code{man_nystrom_model} object returned by
#'   \code{\link{fit_man_nystrom}}.
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
#' At inference, no statistic is estimated from the scored rows / test batch: no
#' test-batch centering, renormalization, bandwidth update, graph update, or head
#' refit is performed. When the present feature set is a strict subset of the
#' training universe, the landmark diffusion geometry IS re-derived over the
#' present columns -- but solely from the frozen training landmarks, never from
#' any scored row -- so row-equivariance is preserved; at full feature overlap
#' this re-derivation reproduces the fit-time geometry exactly. Partial-overlap
#' scoring is therefore an approximate graceful degradation (the scored row's
#' per-sample rCLR is re-centred over the present features while the frozen
#' landmark coordinates are sub-columned, and the re-derived eigenbasis need not
#' align with the fit-time head), not an exact projection. The rCLR
#' representation is exactly invariant to multiplying one specimen by a positive
#' scalar because
#' \eqn{\log(c v_j)-mean_k\log(c v_k)=\log(v_j)-mean_k\log(v_k)} over that
#' specimen's own positive present features.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(60 * 30, shape = 3), nrow = 60,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 30)
#' X[y == 1, 1:6] <- X[y == 1, 1:6] + 20
#' X[y == 0, 7:12] <- X[y == 0, 7:12] + 20
#' model <- fit_man_nystrom(X, y)
#' score_man_nystrom(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Coifman RR, Lafon S. (2006) Diffusion maps. \emph{Applied and Computational
#' Harmonic Analysis} 21: 5-30.
#'
#' Bengio Y, Paiement JF, Vincent P, Delalleau O, Le Roux N, Ouimet M. (2004)
#' Out-of-sample extensions for LLE, Isomap, MDS, Eigenmaps, and Spectral
#' Clustering. \emph{Advances in Neural Information Processing Systems} 16.
#'
#' @export
score_man_nystrom <- function(model, X, meta = NULL) {
  if (!inherits(model, "man_nystrom_model")) {
    stop("score_man_nystrom: model must have class man_nystrom_model")
  }
  X <- .reo_check_matrix(X, "score_man_nystrom", "X")
  .reo_check_meta(meta, nrow(X), "score_man_nystrom", "meta")

  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  X_use <- X[, present, drop = FALSE]
  L_present <- model$landmarks[, present, drop = FALSE]
  geom <- .man_nystrom_geometry(L_present, model$gamma, model$hp)

  out <- vapply(seq_len(nrow(X_use)), function(i) {
    x_i <- X_use[i, ]
    names(x_i) <- present
    z_i <- .man_nystrom_rclr_row(x_i)
    coord <- .man_nystrom_embed_row(z_i, L_present, model$gamma, geom,
                                    model$hp)
    coord <- .man_nystrom_align_coords(coord, model$head$center)
    coord_std <- (coord - model$head$center) / model$head$scale
    sum(coord_std * model$head$coef) + model$head$intercept
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_man_nystrom: scorer produced non-finite or wrong-length output")
  }
  out
}

.man_nystrom_resolve_hp <- function(hp, n) {
  if (!is.list(hp)) stop("fit_man_nystrom: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(is.na(nm)) || any(!nzchar(nm))) {
      stop("fit_man_nystrom: hp fields must be named")
    }
  }
  allowed <- c("n_landmarks", "n_components", "t", "gamma", "eig_floor",
               "min_features", "eps", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_man_nystrom: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  n_landmarks <- hp$n_landmarks
  if (is.null(n_landmarks)) n_landmarks <- min(200L, n)
  n_landmarks <- .man_nystrom_check_positive_integer(
    n_landmarks, "n_landmarks", "fit_man_nystrom"
  )
  if (n_landmarks < 2L) {
    stop("fit_man_nystrom: hp$n_landmarks must be at least 2")
  }

  n_components <- hp$n_components
  if (is.null(n_components)) n_components <- 10L
  n_components <- .man_nystrom_check_positive_integer(
    n_components, "n_components", "fit_man_nystrom"
  )

  t <- hp$t
  if (is.null(t)) t <- 1L
  t <- .man_nystrom_check_positive_integer(t, "t", "fit_man_nystrom")

  gamma <- hp$gamma
  if (!is.null(gamma) &&
      (!is.numeric(gamma) || length(gamma) != 1L ||
       !is.finite(gamma) || gamma <= 0)) {
    stop("fit_man_nystrom: hp$gamma must be NULL or a positive finite number")
  }

  eig_floor <- hp$eig_floor
  if (is.null(eig_floor)) eig_floor <- 1e-6
  if (!is.numeric(eig_floor) || length(eig_floor) != 1L ||
      !is.finite(eig_floor) || eig_floor < 0) {
    stop("fit_man_nystrom: hp$eig_floor must be a non-negative finite number")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  min_features <- .man_nystrom_check_positive_integer(
    min_features, "min_features", "fit_man_nystrom"
  )

  eps <- hp$eps
  if (is.null(eps)) eps <- 1e-8
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_man_nystrom: hp$eps must be a positive finite number")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
    stop("fit_man_nystrom: hp$seed must be a non-negative integer")
  }

  list(
    n_landmarks = n_landmarks,
    n_components = n_components,
    t = t,
    gamma = if (is.null(gamma)) NULL else as.numeric(gamma),
    eig_floor = as.numeric(eig_floor),
    min_features = min_features,
    eps = as.numeric(eps),
    seed = as.integer(seed)
  )
}

.man_nystrom_check_positive_integer <- function(x, field, fn) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1L ||
      x > .Machine$integer.max || x != floor(x)) {
    stop(fn, ": hp$", field, " must be a positive integer")
  }
  as.integer(x)
}

.man_nystrom_landmark_indices <- function(n, n_landmarks, seed) {
  m <- min(n_landmarks, n)
  if (m >= n) return(seq_len(n))

  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)

  set.seed(seed)
  sort(sample(seq_len(n), m, replace = FALSE))
}

.man_nystrom_rclr_row <- function(v) {
  z <- numeric(length(v))
  pos <- which(v > 0)
  if (length(pos) > 0L) {
    g <- mean(log(v[pos]))
    z[pos] <- log(v[pos]) - g
  }
  z
}

.man_nystrom_rclr_matrix <- function(X) {
  Z <- do.call(rbind, lapply(seq_len(nrow(X)), function(i) {
    .man_nystrom_rclr_row(X[i, ])
  }))
  dimnames(Z) <- dimnames(X)
  storage.mode(Z) <- "double"
  Z
}

.man_nystrom_pairwise_d2 <- function(L) {
  ss <- rowSums(L^2)
  D2 <- outer(ss, ss, "+") - 2 * tcrossprod(L)
  D2 <- pmax(D2, 0)
  dimnames(D2) <- list(rownames(L), rownames(L))
  D2
}

.man_nystrom_median_gamma <- function(L, eps) {
  D2 <- .man_nystrom_pairwise_d2(L)
  off <- D2[upper.tri(D2)]
  med <- stats::median(off, na.rm = TRUE)
  1 / (2 * max(med, eps))
}

.man_nystrom_geometry <- function(L, gamma, hp) {
  D2 <- .man_nystrom_pairwise_d2(L)
  W <- exp(-gamma * D2)
  diag(W) <- 1
  q <- rowSums(W)
  Qm12 <- 1 / sqrt(pmax(q, hp$eps))
  A <- sweep(sweep(W, 1L, Qm12, "*"), 2L, Qm12, "*")
  A <- (A + t(A)) / 2

  e <- eigen(A, symmetric = TRUE)
  idx <- order(-e$values, seq_along(e$values), method = "radix")
  lambda <- as.numeric(e$values[idx])
  V <- e$vectors[, idx, drop = FALSE]

  if (length(lambda) <= 1L) {
    V_keep <- matrix(0, nrow = nrow(L), ncol = 0L)
    lambda_keep <- numeric(0)
  } else {
    lambda_rest <- lambda[-1L]
    V_rest <- V[, -1L, drop = FALSE]
    keep <- which(lambda_rest > hp$eig_floor)
    if (length(keep) > hp$n_components) keep <- keep[seq_len(hp$n_components)]
    V_keep <- V_rest[, keep, drop = FALSE]
    V_keep <- .man_nystrom_fix_eigen_signs(V_keep)
    lambda_keep <- lambda_rest[keep]
  }

  # No retained diffusion modes (a high eig_floor, m<=1 landmarks, or
  # near-degenerate -- e.g. collinear -- landmarks where every non-trivial
  # eigenvalue is below eig_floor) yields a 0-column psi. Build it explicitly so
  # the dimension is never dropped, and guard the colnames assignment (the R
  # gotcha paste0("dm", seq_len(0)) is "dm", length 1, not character(0)). The
  # downstream Nystrom embed (K==0 -> numeric(0)) and head (ncol(C)==0 ->
  # constant head) then produce the intended neutral-0 score, at fit AND at
  # score time (this helper is re-derived over the present subset when scoring).
  if (ncol(V_keep) == 0L) {
    psi <- matrix(numeric(0), nrow = nrow(L), ncol = 0L)
  } else {
    psi <- sweep(V_keep, 1L, Qm12, "*")
    colnames(psi) <- paste0("dm", seq_len(ncol(psi)))
  }
  names(lambda_keep) <- colnames(psi)
  list(
    q = as.numeric(q),
    psi = psi,
    lambda = as.numeric(lambda_keep),
    Qm12 = as.numeric(Qm12)
  )
}

.man_nystrom_fix_eigen_signs <- function(V) {
  for (j in seq_len(ncol(V))) {
    idx <- which.max(abs(V[, j]))
    if (is.finite(V[idx, j]) && V[idx, j] < 0) {
      V[, j] <- -V[, j]
    }
  }
  V
}

.man_nystrom_embed_row <- function(z, L, gamma, geom, hp) {
  K <- length(geom$lambda)
  if (K == 0L) return(numeric(0))
  d2 <- rowSums(sweep(L, 2L, z, "-")^2)
  w <- exp(-gamma * pmax(d2, 0))
  s <- sum(w)
  p <- if (!is.finite(s) || s <= 0) rep(1 / nrow(L), nrow(L)) else w / s
  as.numeric(p %*% geom$psi) * (geom$lambda^(hp$t - 1L))
}

.man_nystrom_embed_matrix <- function(Z, L, gamma, geom, hp) {
  K <- length(geom$lambda)
  if (K == 0L) {
    return(matrix(0, nrow = nrow(Z), ncol = 0L,
                  dimnames = list(rownames(Z), character(0))))
  }
  E <- t(vapply(seq_len(nrow(Z)), function(i) {
    .man_nystrom_embed_row(Z[i, ], L, gamma, geom, hp)
  }, numeric(K)))
  dimnames(E) <- list(rownames(Z), paste0("dm", seq_len(K)))
  storage.mode(E) <- "double"
  E
}

.man_nystrom_align_coords <- function(coord, center) {
  K <- length(center)
  if (K == 0L) return(numeric(0))
  out <- center
  n <- min(length(coord), K)
  if (n > 0L) out[seq_len(n)] <- coord[seq_len(n)]
  out
}

.man_nystrom_fit_head <- function(C, y, hp) {
  if (ncol(C) == 0L) {
    return(list(
      type = "constant",
      coef = numeric(0),
      intercept = 0,
      center = numeric(0),
      scale = numeric(0)
    ))
  }

  center <- colMeans(C)
  scale <- vapply(seq_len(ncol(C)), function(j) stats::sd(C[, j]), numeric(1))
  scale[!is.finite(scale) | scale < hp$eps] <- hp$eps
  C_std <- sweep(sweep(C, 2L, center, "-"), 2L, scale, "/")
  colnames(C_std) <- paste0("dm", seq_len(ncol(C_std)))

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
      coef = as.numeric(beta[-1L]),
      intercept = as.numeric(beta[1L]),
      center = as.numeric(center),
      scale = as.numeric(scale)
    ))
  }

  lda <- .man_nystrom_lda_head(C_std, y, hp$eps)
  list(
    type = "diag_lda",
    coef = as.numeric(lda$coef),
    intercept = as.numeric(lda$intercept),
    center = as.numeric(center),
    scale = as.numeric(scale)
  )
}

.man_nystrom_lda_head <- function(C_std, y, eps) {
  case <- C_std[y == 1L, , drop = FALSE]
  control <- C_std[y == 0L, , drop = FALSE]
  mu_case <- colMeans(case)
  mu_control <- colMeans(control)
  var_case <- .man_nystrom_col_var(case)
  var_control <- .man_nystrom_col_var(control)
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

.man_nystrom_col_var <- function(X) {
  if (nrow(X) < 2L) return(rep(0, ncol(X)))
  as.numeric(vapply(seq_len(ncol(X)), function(j) stats::var(X[, j]),
                    numeric(1)))
}
