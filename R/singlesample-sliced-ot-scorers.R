#' @title Fit sliced random-projection Gaussian LRT discriminator
#'
#' @description
#' Learns a frozen two-class discriminator from training data by projecting each
#' specimen's within-sample robust centered log-ratio (rCLR) vector onto
#' deterministic random one-dimensional slices. Along every slice, the case and
#' control training anchors are summarized by frozen Gaussian moments, and
#' scoring uses the average Gaussian log-likelihood ratio. Larger scores are
#' more case-like.
#'
#' The rCLR transform is computed separately within each specimen:
#' \code{z[pos] = log(v[pos]) - mean(log(v[pos]))} for positive abundances and
#' \code{z[!pos] = 0}. This makes the representation invariant to multiplying a
#' specimen by a positive scalar, without estimating any statistic from the
#' scoring batch.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of
#'   non-negative abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{n_slices} positive integer
#'   number of random projection directions (default \code{100L});
#'   \code{max_anchors_per_class} positive integer anchor cap per class (default
#'   \code{200L}); \code{min_features} positive integer feature-overlap floor
#'   for fitting and scoring (default \code{3L}); \code{eps} positive finite
#'   variance floor (default \code{1e-6}); and \code{seed} non-negative integer
#'   used for deterministic anchor subsampling and random directions while
#'   restoring the global RNG state (default \code{1L}).
#'
#' @return A plain list of class \code{ot_slicedlrt_model} containing the frozen
#'   \code{feature_universe}, raw \code{case_anchors}, raw
#'   \code{control_anchors}, the full-universe sliced Gaussian representation
#'   \code{repr}, and resolved \code{hp}.
#'
#' @details
#' Fitting stores raw training anchors rather than only fitted slice moments so
#' that scoring can re-derive the same representation over the feature subset
#' present at inference. The helper \code{.ot_slicedlrt_repr()} is used both at
#' fit time over the full training feature universe and at score time over
#' \code{intersect(model$feature_universe, colnames(X))}. Its only inputs are
#' frozen anchors, frozen hyperparameters, and the present-feature dimension.
#' Consequently, no test-batch centering, quantile, variance, bandwidth, or
#' other cross-row statistic is estimated from scored specimens. If the present
#' overlap contains fewer than \code{hp$min_features} features,
#' \code{\link{score_ot_slicedlrt}} returns the documented neutral score
#' \code{0} for every row.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(80 * 30, shape = 2), nrow = 80,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 40)
#' X[y == 1, 1:6] <- X[y == 1, 1:6] + 25
#' X[y == 0, 7:12] <- X[y == 0, 7:12] + 25
#' model <- fit_ot_slicedlrt(X, y)
#' score <- score_ot_slicedlrt(model, X)
#' }
#'
#' @references
#' Bonneel N, Rabin J, Peyre G, Pfister H. (2015) Sliced and Radon Wasserstein
#' barycenters of measures. \emph{Journal of Mathematical Imaging and Vision}
#' 51: 22-45.
#'
#' Neyman J, Pearson ES. (1933) On the problem of the most efficient tests of
#' statistical hypotheses. \emph{Philosophical Transactions of the Royal Society
#' of London. Series A} 231: 289-337.
#'
#' @export
fit_ot_slicedlrt <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_ot_slicedlrt", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_ot_slicedlrt", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_ot_slicedlrt")
  hp <- .ot_slicedlrt_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_ot_slicedlrt: X_train must contain at least hp$min_features features")
  }

  anchor_idx <- .ot_slicedlrt_anchor_indices(
    y = y,
    max_anchors_per_class = hp$max_anchors_per_class,
    seed = hp$seed
  )
  case_anchors <- X_train[anchor_idx$case, , drop = FALSE]
  control_anchors <- X_train[anchor_idx$control, , drop = FALSE]

  # Stored for inspection / tests. The score path re-derives the same helper
  # over the present-feature subset from frozen raw anchors and the frozen seed.
  repr <- .ot_slicedlrt_repr(
    case_anchors = case_anchors,
    control_anchors = control_anchors,
    hp = hp,
    present_dim = ncol(X_train)
  )

  model <- list(
    feature_universe = colnames(X_train),
    case_anchors = case_anchors,
    control_anchors = control_anchors,
    repr = repr,
    hp = hp
  )
  class(model) <- "ot_slicedlrt_model"
  model
}


#' @title Score sliced random-projection Gaussian LRT discriminator
#'
#' @description
#' Scores each row of \code{X} independently with a model learned by
#' \code{\link{fit_ot_slicedlrt}}. Each specimen is transformed to its own rCLR
#' vector over the features shared with the frozen training universe, projected
#' onto deterministic random unit directions, and assigned the average
#' per-direction Gaussian log-likelihood ratio of the case class versus the
#' control class. Larger values are more case-like.
#'
#' @param model An \code{ot_slicedlrt_model} object returned by
#'   \code{\link{fit_ot_slicedlrt}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. If fewer than
#'   \code{model$hp$min_features} model-universe features are present in
#'   \code{X}, returns a neutral vector of zeros.
#'
#' @details
#' At scoring time, the present feature set is
#' \code{intersect(model$feature_universe, colnames(X))}. The sliced Gaussian
#' representation is then re-derived exactly once from frozen case/control
#' anchors restricted to that present set plus the frozen seed. A scored row
#' contributes only its own rCLR vector; no statistic is estimated from any other
#' row of \code{X}. This row-local computation is compatible with singleton
#' deployment and is invariant to row order, duplication, and batch composition.
#'
#' The default \code{hp$min_features} is \code{3L}. Note that a present overlap of
#' exactly one feature is degenerate: the single-feature per-sample rCLR is
#' identically zero (\eqn{\log v - \log v = 0}), so every projection is zero and
#' the score is a constant neutral \code{0} with no discrimination. Keep
#' \code{min_features} \eqn{\ge 2} (the default \code{3L} already does) for a
#' discriminating configuration.
#'
#' For one row with projection coordinates \eqn{s_m}, case moments
#' \eqn{\mu_{1m}, \sigma^2_{1m}}, and control moments
#' \eqn{\mu_{0m}, \sigma^2_{0m}}, the score is
#' \deqn{\frac{1}{M} \sum_m
#'   \left[-\frac{(s_m-\mu_{1m})^2}{2\sigma^2_{1m}}
#'          -\frac{1}{2}\log\sigma^2_{1m}
#'          +\frac{(s_m-\mu_{0m})^2}{2\sigma^2_{0m}}
#'          +\frac{1}{2}\log\sigma^2_{0m}\right].}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(80 * 30, shape = 2), nrow = 80,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 40)
#' X[y == 1, 1:6] <- X[y == 1, 1:6] + 25
#' X[y == 0, 7:12] <- X[y == 0, 7:12] + 25
#' model <- fit_ot_slicedlrt(X, y)
#' score_ot_slicedlrt(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Bonneel N, Rabin J, Peyre G, Pfister H. (2015) Sliced and Radon Wasserstein
#' barycenters of measures. \emph{Journal of Mathematical Imaging and Vision}
#' 51: 22-45.
#'
#' Neyman J, Pearson ES. (1933) On the problem of the most efficient tests of
#' statistical hypotheses. \emph{Philosophical Transactions of the Royal Society
#' of London. Series A} 231: 289-337.
#'
#' @export
score_ot_slicedlrt <- function(model, X, meta = NULL) {
  if (!inherits(model, "ot_slicedlrt_model")) {
    stop("score_ot_slicedlrt: model must have class ot_slicedlrt_model")
  }
  X <- .reo_check_matrix(X, "score_ot_slicedlrt", "X")
  .reo_check_meta(meta, nrow(X), "score_ot_slicedlrt", "meta")

  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  repr <- .ot_slicedlrt_repr(
    case_anchors = model$case_anchors[, present, drop = FALSE],
    control_anchors = model$control_anchors[, present, drop = FALSE],
    hp = model$hp,
    present_dim = length(present)
  )
  X_use <- X[, present, drop = FALSE]

  out <- vapply(seq_len(nrow(X_use)), function(i) {
    x_i <- X_use[i, ]
    names(x_i) <- present
    z <- .ot_slicedlrt_rclr(x_i)
    s <- as.numeric(repr$Theta %*% z)
    .ot_slicedlrt_gaussian_llr(s, repr)
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_ot_slicedlrt: scorer produced non-finite or wrong-length output")
  }
  out
}

.ot_slicedlrt_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_ot_slicedlrt: hp must be a list")
  if (length(hp) > 0L) {
    if (is.null(names(hp)) || any(!nzchar(names(hp)))) {
      stop("fit_ot_slicedlrt: hp fields must be named")
    }
    if (any(duplicated(names(hp)))) {
      stop("fit_ot_slicedlrt: hp fields must be unique")
    }
  }
  allowed <- c("n_slices", "max_anchors_per_class", "min_features", "eps",
               "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_ot_slicedlrt: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  n_slices <- hp$n_slices
  if (is.null(n_slices)) n_slices <- 100L
  if (!is.numeric(n_slices) || length(n_slices) != 1L ||
      !is.finite(n_slices) || n_slices < 1L ||
      n_slices > .Machine$integer.max ||
      n_slices != as.integer(n_slices)) {
    stop("fit_ot_slicedlrt: hp$n_slices must be a positive integer")
  }

  max_anchors_per_class <- hp$max_anchors_per_class
  if (is.null(max_anchors_per_class)) max_anchors_per_class <- 200L
  if (!is.numeric(max_anchors_per_class) ||
      length(max_anchors_per_class) != 1L ||
      !is.finite(max_anchors_per_class) || max_anchors_per_class < 1L ||
      max_anchors_per_class > .Machine$integer.max ||
      max_anchors_per_class != as.integer(max_anchors_per_class)) {
    stop("fit_ot_slicedlrt: hp$max_anchors_per_class must be a positive integer")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_ot_slicedlrt: hp$min_features must be a positive integer")
  }

  eps <- hp$eps
  if (is.null(eps)) eps <- 1e-6
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_ot_slicedlrt: hp$eps must be a positive finite number")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_ot_slicedlrt: hp$seed must be a non-negative integer")
  }

  list(
    n_slices = as.integer(n_slices),
    max_anchors_per_class = as.integer(max_anchors_per_class),
    min_features = as.integer(min_features),
    eps = as.numeric(eps),
    seed = as.integer(seed)
  )
}

.ot_slicedlrt_anchor_indices <- function(y, max_anchors_per_class, seed) {
  case_idx <- which(y == 1L)
  control_idx <- which(y == 0L)
  needs_subsample <- length(case_idx) > max_anchors_per_class ||
    length(control_idx) > max_anchors_per_class
  if (!needs_subsample) {
    return(list(case = case_idx, control = control_idx))
  }

  .ot_slicedlrt_with_seed(seed, {
    if (length(case_idx) > max_anchors_per_class) {
      case_idx <- sort(sample(case_idx, max_anchors_per_class, replace = FALSE))
    }
    if (length(control_idx) > max_anchors_per_class) {
      control_idx <- sort(sample(control_idx, max_anchors_per_class,
                                 replace = FALSE))
    }
    list(case = case_idx, control = control_idx)
  })
}

.ot_slicedlrt_with_seed <- function(seed, expr) {
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
  force(expr)
}

.ot_slicedlrt_repr <- function(case_anchors, control_anchors, hp, present_dim) {
  if (!is.matrix(case_anchors) || !is.matrix(control_anchors)) {
    stop(".ot_slicedlrt_repr: anchors must be matrices")
  }
  if (!is.numeric(present_dim) || length(present_dim) != 1L ||
      !is.finite(present_dim) || present_dim < 1L ||
      present_dim > .Machine$integer.max ||
      present_dim != as.integer(present_dim)) {
    stop(".ot_slicedlrt_repr: present_dim must be a positive integer")
  }
  present_dim <- as.integer(present_dim)
  if (ncol(case_anchors) != present_dim ||
      ncol(control_anchors) != present_dim) {
    stop(".ot_slicedlrt_repr: present_dim must match anchor columns")
  }
  if (hp$n_slices > floor(.Machine$integer.max / present_dim)) {
    stop(".ot_slicedlrt_repr: hp$n_slices * present_dim is too large")
  }

  Theta <- .ot_slicedlrt_with_seed(hp$seed, {
    matrix(stats::rnorm(hp$n_slices * present_dim),
           nrow = hp$n_slices, ncol = present_dim)
  })
  theta_norm <- sqrt(rowSums(Theta^2))
  bad <- !is.finite(theta_norm) | theta_norm <= 0
  if (any(bad)) {
    Theta[bad, ] <- 0
    Theta[bad, 1L] <- 1
    theta_norm[bad] <- 1
  }
  Theta <- sweep(Theta, 1L, theta_norm, "/")

  Z_case <- .ot_slicedlrt_rclr_matrix(case_anchors)
  Z_control <- .ot_slicedlrt_rclr_matrix(control_anchors)
  P_case <- Z_case %*% t(Theta)
  P_control <- Z_control %*% t(Theta)

  list(
    Theta = Theta,
    mu_case = as.numeric(colMeans(P_case)),
    v_case = .ot_slicedlrt_col_var(P_case, hp$eps),
    mu_control = as.numeric(colMeans(P_control)),
    v_control = .ot_slicedlrt_col_var(P_control, hp$eps)
  )
}

.ot_slicedlrt_rclr <- function(x) {
  z <- numeric(length(x))
  pos <- x > 0
  if (any(pos)) {
    log_pos <- log(x[pos])
    z[pos] <- log_pos - mean(log_pos)
  }
  names(z) <- names(x)
  z
}

.ot_slicedlrt_rclr_matrix <- function(X) {
  Z <- do.call(rbind, lapply(seq_len(nrow(X)), function(i) {
    .ot_slicedlrt_rclr(X[i, ])
  }))
  dimnames(Z) <- dimnames(X)
  storage.mode(Z) <- "double"
  Z
}

.ot_slicedlrt_col_var <- function(P, eps) {
  if (nrow(P) < 2L) {
    return(rep(eps, ncol(P)))
  }
  v <- if (ncol(P) == 1L) {
    stats::var(P[, 1L])
  } else {
    apply(P, 2L, stats::var)
  }
  v <- as.numeric(v)
  v[!is.finite(v)] <- 0
  pmax(v, eps)
}

.ot_slicedlrt_gaussian_llr <- function(s, repr) {
  mean(
    -((s - repr$mu_case)^2) / (2 * repr$v_case) -
      0.5 * log(repr$v_case) +
      ((s - repr$mu_control)^2) / (2 * repr$v_control) +
      0.5 * log(repr$v_control)
  )
}
