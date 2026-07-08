#' @title Fit unconstrained least-squares importance-fitting (uLSIF) density-ratio LRT
#'
#' @description
#' Learns a frozen two-class discriminator from training data by \emph{directly}
#' estimating the class-conditional density ratio
#' \eqn{w(x) = p_{\mathrm{case}}(x) / p_{\mathrm{control}}(x)} -- without ever
#' estimating either density -- using unconstrained Least-Squares Importance
#' Fitting (uLSIF; Kanamori, Hido & Sugiyama 2009). The estimated ratio is itself
#' the Neyman-Pearson-optimal class-conditional likelihood ratio, so larger
#' \eqn{w(z)} means more case-like.
#'
#' Each specimen is first mapped to a self-contained, per-sample robust
#' centred-log-ratio (rCLR) representation \eqn{z}, where
#' \eqn{z_j = \log v_j - \mathrm{mean}_{k:\,v_k>0}\log v_k} on the nonzero parts
#' of the specimen and \eqn{z_j = 0} on its zero parts. Because the centring uses
#' only that specimen's own nonzero values, \eqn{z(c\,v) = z(v)} for any
#' \eqn{c>0}: the representation -- and therefore every downstream score -- is
#' exactly invariant to per-sample scaling (library size / input amount).
#'
#' The density ratio is modelled as a non-negative combination of Gaussian
#' basis functions \eqn{\phi_l(z) = \exp(-\lVert z - c_l\rVert^2 / (2\sigma^2))}
#' centred at \eqn{b} frozen kernel centres \eqn{c_l}. The centres are a
#' deterministic frozen subsample of the \emph{case} rCLR anchors (cap
#' \code{hp$n_centers}); the bandwidth \eqn{\sigma} is either supplied as
#' \code{hp$sigma} or learned by the median heuristic on the pairwise Euclidean
#' distances among the centres. Writing
#' \eqn{\Phi_{\mathrm{control}}} and \eqn{\Phi_{\mathrm{case}}} for the kernel
#' matrices of the control / case anchors against the centres, the closed-form
#' uLSIF coefficients solve the ridge-regularized normal equations
#' \deqn{(H + \lambda I)\,\alpha = h,\qquad
#'       H = \tfrac{1}{n_{\mathrm{control}}}
#'           \Phi_{\mathrm{control}}^\top \Phi_{\mathrm{control}},\qquad
#'       h = \mathrm{colMeans}(\Phi_{\mathrm{case}}),}
#' and the non-negative parts \eqn{\alpha \leftarrow \max(\alpha, 0)} are kept so
#' that \eqn{\hat w(z) = \sum_l \alpha_l \phi_l(z) \ge 0} is a valid density
#' ratio. There is no iteration: \eqn{H} is a Gram matrix (positive
#' semidefinite), so \eqn{H + \lambda I} is positive definite for \eqn{\lambda>0}
#' and the single linear solve is exact (an increase-to-positive-definite ridge
#' loop guards only pathological numerical rank deficiency).
#'
#' Fitting stores the raw case/control anchor rows and the raw centre rows (so the
#' kernel model can be re-derived consistently over any present-feature subset at
#' scoring), the frozen training feature universe, the resolved hyperparameters,
#' and -- for inspection and tests -- the full-universe representation
#' (\code{centers}, \code{sigma}, \code{alpha}). No test data and no test-time
#' batch statistic are used.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{sigma} Gaussian kernel
#'   bandwidth (default \code{NULL}, the median heuristic; otherwise a positive
#'   finite number); \code{lambda} positive uLSIF ridge (default \code{1e-3});
#'   \code{n_centers} positive-integer cap on the number of case-anchor kernel
#'   centres (default \code{100L}); \code{max_anchors_per_class} positive-integer
#'   anchor cap per class (default \code{200L}); \code{min_features}
#'   positive-integer feature-overlap floor at scoring (default \code{3L});
#'   \code{eps} positive ridge increment / bandwidth-and-score floor (default
#'   \code{1e-6}); and \code{seed} non-negative integer used only for
#'   deterministic anchor / centre subsampling while restoring the global RNG
#'   state (default \code{1L}).
#'
#' @return A plain list of class \code{dre_ulsif_model} containing
#'   \code{case_anchors}, \code{control_anchors}, \code{center_anchors},
#'   \code{feature_universe}, \code{repr} (the full-universe frozen
#'   centres/bandwidth/coefficients), and the resolved \code{hp}.
#'
#' @details
#' The score (see \code{\link{score_dre_ulsif}}) is the log density ratio
#' \eqn{S(z) = \log\max(\hat w(z), \varepsilon)} with
#' \eqn{\hat w(z) = \sum_l \alpha_l \exp(-\lVert z - c_l\rVert^2 / (2\sigma^2))}.
#' The log scale is the canonical likelihood-ratio-test statistic
#' \eqn{\log(p_{\mathrm{case}}/p_{\mathrm{control}})}; it is a strictly monotone
#' transform of \eqn{\hat w} (so it preserves the density-ratio ranking and
#' discrimination) but compresses the heavy right tail of the ratio, which keeps
#' the statistic numerically stable. The \eqn{\varepsilon} floor (which is only
#' ever active when the non-negative kernel model evaluates to (near) zero) keeps
#' the logarithm finite.
#'
#' Degenerate inputs stay finite without caller-side special-casing. A feature
#' that is constant across a class's anchors contributes a constant rCLR column
#' that simply shifts the kernel arguments; the Gram matrix \eqn{H} remains
#' positive semidefinite and the ridge keeps \eqn{H + \lambda I} invertible. If
#' the kernel centres collapse to a single point (or are all identical), the
#' median heuristic has no pairwise distance to use and the bandwidth falls back
#' to a within-centre root-mean-square length scale floored at \code{hp$eps}, so
#' \eqn{\sigma > 0} always.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120; p <- 30; k <- 10
#' L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2  # case density shift
#' X <- exp(L)
#' model <- fit_dre_ulsif(X, y)
#' score <- score_dre_ulsif(model, X)
#' }
#'
#' @references
#' Kanamori T, Hido S, Sugiyama M. (2009) A least-squares approach to direct
#' importance estimation. \emph{Journal of Machine Learning Research}
#' 10: 1391-1445.
#'
#' Sugiyama M, Suzuki T, Kanamori T. (2012) \emph{Density Ratio Estimation in
#' Machine Learning}. Cambridge University Press.
#'
#' @export
fit_dre_ulsif <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_dre_ulsif", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_dre_ulsif", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_dre_ulsif")
  hp <- .dre_ulsif_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_dre_ulsif: X_train must contain at least hp$min_features features")
  }

  sel <- .dre_ulsif_select_indices(
    y = y,
    max_anchors_per_class = hp$max_anchors_per_class,
    n_centers = hp$n_centers,
    seed = hp$seed
  )
  case_anchors <- X_train[sel$case, , drop = FALSE]
  control_anchors <- X_train[sel$control, , drop = FALSE]
  center_anchors <- case_anchors[sel$center_local, , drop = FALSE]

  # Full-universe representation: stored for inspection / tests. The score path
  # recomputes this over the present-feature subset through the SAME helper, so
  # there is one kernel/coefficient definition and no fit/score drift.
  repr <- .dre_ulsif_repr(case_anchors, control_anchors, center_anchors, hp)

  model <- list(
    case_anchors = case_anchors,
    control_anchors = control_anchors,
    center_anchors = center_anchors,
    feature_universe = colnames(X_train),
    repr = repr,
    hp = hp
  )
  class(model) <- "dre_ulsif_model"
  model
}


#' @title Score the uLSIF density-ratio LRT discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model learned by
#' \code{\link{fit_dre_ulsif}}. For one specimen the abundances are mapped to the
#' self-contained per-sample rCLR representation \eqn{z}, and the score is the log
#' estimated density ratio
#' \deqn{S(z) = \log\max\!\Bigl(\sum_l \alpha_l
#'             \exp\bigl(-\lVert z - c_l\rVert^2 / (2\sigma^2)\bigr),\,
#'             \varepsilon\Bigr),}
#' with larger values more case-like.
#'
#' At scoring time the present feature set is
#' \code{intersect(model$feature_universe, colnames(X))}. The kernel centres,
#' bandwidth and coefficients are re-derived over exactly that present set from
#' the frozen raw anchors / centres -- using only frozen training data, never a
#' scored-batch statistic -- which keeps partial-overlap scores consistent and
#' equal to the fit-time representation at full overlap. The score of a row
#' therefore depends only on that row and the frozen model, and is exactly
#' invariant to per-sample scaling. If fewer than \code{model$hp$min_features}
#' features are present, the documented neutral score \code{0} is returned for
#' every row.
#'
#' @param model A \code{dre_ulsif_model} object returned by
#'   \code{\link{fit_dre_ulsif}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Larger values
#'   are more case-like (the specimen lies where the case density dominates the
#'   control density). Scoring uses only each row's own values plus the frozen
#'   centres, bandwidth and coefficients.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120; p <- 30; k <- 10
#' L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
#' X <- exp(L)
#' model <- fit_dre_ulsif(X, y)
#' score_dre_ulsif(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Kanamori T, Hido S, Sugiyama M. (2009) A least-squares approach to direct
#' importance estimation. \emph{Journal of Machine Learning Research}
#' 10: 1391-1445.
#'
#' @export
score_dre_ulsif <- function(model, X, meta = NULL) {
  if (!inherits(model, "dre_ulsif_model")) {
    stop("score_dre_ulsif: model must have class dre_ulsif_model")
  }
  X <- .reo_check_matrix(X, "score_dre_ulsif", "X")
  .reo_check_meta(meta, nrow(X), "score_dre_ulsif", "meta")

  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  repr <- .dre_ulsif_repr(
    model$case_anchors[, present, drop = FALSE],
    model$control_anchors[, present, drop = FALSE],
    model$center_anchors[, present, drop = FALSE],
    model$hp
  )
  X_use <- X[, present, drop = FALSE]
  C <- repr$centers
  sigma <- repr$sigma
  alpha <- repr$alpha
  eps <- model$hp$eps

  out <- vapply(seq_len(nrow(X_use)), function(i) {
    z <- .dre_ulsif_rclr(X_use[i, ])
    k <- .dre_ulsif_kernel(matrix(z, nrow = 1L), C, sigma)
    w <- sum(alpha * as.numeric(k))
    log(max(w, eps))
  }, numeric(1))

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_dre_ulsif: scorer produced non-finite or wrong-length output")
  }
  out
}

.dre_ulsif_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_dre_ulsif: hp must be a list")
  allowed <- c("sigma", "lambda", "n_centers", "max_anchors_per_class",
               "min_features", "eps", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_dre_ulsif: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  sigma <- hp$sigma
  if (!is.null(sigma) &&
      (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) ||
       sigma <= 0)) {
    stop("fit_dre_ulsif: hp$sigma must be NULL or a positive finite number")
  }

  lambda <- hp$lambda
  if (is.null(lambda)) lambda <- 1e-3
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) ||
      lambda <= 0) {
    stop("fit_dre_ulsif: hp$lambda must be a positive finite number")
  }

  n_centers <- hp$n_centers
  if (is.null(n_centers)) n_centers <- 100L
  if (!is.numeric(n_centers) || length(n_centers) != 1L ||
      !is.finite(n_centers) || n_centers < 1L ||
      n_centers > .Machine$integer.max ||
      n_centers != as.integer(n_centers)) {
    stop("fit_dre_ulsif: hp$n_centers must be a positive integer")
  }

  max_anchors_per_class <- hp$max_anchors_per_class
  if (is.null(max_anchors_per_class)) max_anchors_per_class <- 200L
  if (!is.numeric(max_anchors_per_class) ||
      length(max_anchors_per_class) != 1L ||
      !is.finite(max_anchors_per_class) || max_anchors_per_class < 1L ||
      max_anchors_per_class > .Machine$integer.max ||
      max_anchors_per_class != as.integer(max_anchors_per_class)) {
    stop("fit_dre_ulsif: hp$max_anchors_per_class must be a positive integer")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_dre_ulsif: hp$min_features must be a positive integer")
  }

  eps <- hp$eps
  if (is.null(eps)) eps <- 1e-6
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_dre_ulsif: hp$eps must be a positive finite number")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_dre_ulsif: hp$seed must be a non-negative integer")
  }

  list(
    sigma = if (is.null(sigma)) NULL else as.numeric(sigma),
    lambda = as.numeric(lambda),
    n_centers = as.integer(n_centers),
    max_anchors_per_class = as.integer(max_anchors_per_class),
    min_features = as.integer(min_features),
    eps = as.numeric(eps),
    seed = as.integer(seed)
  )
}

# Deterministic per-class anchor subsampling AND kernel-centre selection, with a
# single global-RNG save/restore region (mirrors the lrt-copula / ig-fisherrao
# pattern). The centres are a subsample of the (already subsampled) CASE anchors.
# When no subsampling is needed -- neither class exceeds its anchor cap and the
# case anchors do not exceed the centre cap -- the data are used as-is and the
# RNG is never touched. The sample() calls always run in the fixed order
# case -> control -> centres, so the selection is reproducible from hp$seed.
.dre_ulsif_select_indices <- function(y, max_anchors_per_class, n_centers,
                                      seed) {
  case_idx <- which(y == 1L)
  control_idx <- which(y == 0L)
  n_case_anchor <- min(length(case_idx), max_anchors_per_class)
  needs_subsample <- length(case_idx) > max_anchors_per_class ||
    length(control_idx) > max_anchors_per_class ||
    n_case_anchor > n_centers
  if (!needs_subsample) {
    return(list(case = case_idx, control = control_idx,
                center_local = seq_along(case_idx)))
  }

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
  if (length(case_idx) > max_anchors_per_class) {
    case_idx <- sort(sample(case_idx, max_anchors_per_class, replace = FALSE))
  }
  if (length(control_idx) > max_anchors_per_class) {
    control_idx <- sort(sample(control_idx, max_anchors_per_class,
                               replace = FALSE))
  }
  n_case_anchor <- length(case_idx)
  if (n_case_anchor > n_centers) {
    center_local <- sort(sample(seq_len(n_case_anchor), n_centers,
                                replace = FALSE))
  } else {
    center_local <- seq_len(n_case_anchor)
  }
  list(case = case_idx, control = control_idx, center_local = center_local)
}

# Self-contained per-sample robust CLR of one specimen. The geometric-mean
# centring uses ONLY this specimen's own nonzero parts, so it is a within-sample
# transform: z(c*v) = z(v) for any c > 0 (exact per-sample scale-invariance), and
# zero parts map to 0. An all-zero (or single-nonzero) specimen maps to the
# all-zero vector. This deliberately does NOT reuse any package rclr helper that
# centres using cross-sample / reference statistics.
.dre_ulsif_rclr <- function(v) {
  z <- rep.int(0, length(v))
  pos <- which(v > 0)
  if (length(pos) >= 1L) {
    lv <- log(v[pos])
    z[pos] <- lv - mean(lv)
  }
  z
}

# Row-wise rCLR of an anchor matrix. rbind (not t(vapply)) keeps an n x p shape
# even when n == 1 or p == 1, where the collapsed forms would break alignment.
.dre_ulsif_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .dre_ulsif_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Squared Euclidean distances between every row of Z (n x p) and every centre
# row of C (b x p), returned as an n x b matrix via the identity
# ||z - c||^2 = ||z||^2 + ||c||^2 - 2 z.c. A small negative-clip absorbs the
# floating-point round-off that can make an exact zero distance slightly
# negative.
.dre_ulsif_sqdist <- function(Z, C) {
  zz <- rowSums(Z^2)
  cc <- rowSums(C^2)
  cross <- Z %*% t(C)
  d2 <- outer(zz, cc, "+") - 2 * cross
  d2[d2 < 0] <- 0
  d2
}

# Gaussian kernel matrix exp(-||z - c||^2 / (2 sigma^2)) of rows of Z against the
# centres C. Entries lie in (0, 1], so any non-negative combination is >= 0.
.dre_ulsif_kernel <- function(Z, C, sigma) {
  exp(-.dre_ulsif_sqdist(Z, C) / (2 * sigma^2))
}

# Median-heuristic Gaussian bandwidth: the median of the strictly positive
# pairwise Euclidean distances among the kernel centres, floored at eps so it is
# strictly positive. With a single centre (or all-identical centres) there is no
# pairwise distance; the bandwidth then falls back to a within-centre
# root-mean-square length scale, again floored at eps. Panel-frozen in the sense
# that it is computed once over whichever feature set is in play (full universe at
# fit, present subset at score) purely from frozen training anchors.
.dre_ulsif_sigma <- function(C, eps) {
  if (nrow(C) >= 2L) {
    d <- as.numeric(stats::dist(C, method = "euclidean"))
    d <- d[is.finite(d) & d > 0]
    if (length(d) > 0L) {
      return(max(stats::median(d), eps))
    }
  }
  max(sqrt(mean(C^2)), eps)
}

# Solve the ridge-regularized uLSIF normal equations (H + lambda I) alpha = h.
# H is a Gram matrix (positive semidefinite) so H + lambda I is positive definite
# for lambda > 0 and the first solve succeeds; the loop is the documented guard
# for pathological numerical rank deficiency, increasing the ridge by eps until a
# finite solution is obtained (mirrors the lrt-copula increase-to-PD ridge loop).
.dre_ulsif_solve_alpha <- function(H, h, lambda, eps) {
  b <- nrow(H)
  Imat <- diag(b)
  ridge <- lambda
  alpha <- NULL
  for (k in 0:50) {
    alpha <- tryCatch(solve(H + ridge * Imat, h), error = function(e) NULL)
    if (!is.null(alpha) && all(is.finite(alpha))) break
    alpha <- NULL
    ridge <- ridge + eps
  }
  if (is.null(alpha)) {
    stop("fit_dre_ulsif: failed to solve the uLSIF normal equations (singular H + ridge)")
  }
  alpha
}

# Derive the uLSIF kernel-model representation (centres in rCLR space, bandwidth,
# non-negative coefficients) from raw anchor / centre matrices over their current
# (full or present) feature set. Used by both fit (full universe) and score
# (present subset), so the model is defined once and there is no fit/score drift.
# H = (1/n_control) Phi_control^T Phi_control is the second-moment (Gram) matrix
# of the control kernel features; h = colMeans(Phi_case) is the case mean kernel
# feature; alpha solves (H + lambda I) alpha = h and is clipped to its
# non-negative part so w(z) = sum_l alpha_l phi_l(z) >= 0 is a valid density ratio.
.dre_ulsif_repr <- function(case_anchors, control_anchors, center_anchors, hp) {
  Zc <- .dre_ulsif_rclr_matrix(case_anchors)
  Z0 <- .dre_ulsif_rclr_matrix(control_anchors)
  C <- .dre_ulsif_rclr_matrix(center_anchors)

  sigma <- if (is.null(hp$sigma)) .dre_ulsif_sigma(C, hp$eps) else hp$sigma

  Phi_case <- .dre_ulsif_kernel(Zc, C, sigma)
  Phi_ctrl <- .dre_ulsif_kernel(Z0, C, sigma)

  H <- crossprod(Phi_ctrl) / nrow(Z0)
  h <- colMeans(Phi_case)
  alpha <- .dre_ulsif_solve_alpha(H, h, hp$lambda, hp$eps)
  alpha <- pmax(alpha, 0)

  list(centers = C, sigma = sigma, alpha = alpha)
}
