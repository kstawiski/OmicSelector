#' @title ECOD + COPOD novelty discriminator (single-sample, frozen-ECDF)
#'
#' @description
#' Single-sample-deployable novelty/outlier-detection-as-discriminator scorer
#' (family NC) that bundles the two empirical-cumulative-distribution outlier
#' detectors ECOD (Li et al., 2022) and COPOD (Li et al., 2020) from the Python
#' \pkg{pyod} library. Both detectors flag a specimen as anomalous when its
#' per-feature values fall in the extreme tails of an in-distribution reference,
#' which here is the training CONTROL cohort (mirroring the healthy-reference
#' convention of \code{\link{fit_conformal_anomaly}}).
#'
#' \strong{Why the per-feature ECDF state is frozen rather than re-fitted at
#' score time.} \pkg{pyod}'s \code{ECOD$decision_function()} /
#' \code{COPOD$decision_function()} are \emph{transductive}: they re-fit each
#' per-feature empirical CDF on whatever batch is passed, so the score of one
#' specimen depends on which other specimens are scored alongside it (measured
#' batch-vs-row-by-row max difference \eqn{\approx 0.68} for ECOD). That is a
#' batch-coupled scorer and would fail the single-sample deployability gate
#' (\code{\link{singlesample_assert_row_equivariant}}). This implementation
#' therefore NEVER calls \pkg{pyod} at score time. At fit it freezes the
#' training-control per-feature order statistics and the per-feature skewness
#' signs that ECOD uses; at score it evaluates each specimen's ECOD and COPOD
#' aggregates against those frozen reference ECDFs in pure R, one specimen at a
#' time. On the training-control rows the frozen reference ECDF equals the
#' fit-time ECDF, so the pure-R aggregates reproduce \pkg{pyod}'s
#' \code{decision_scores_} to machine precision (asserted by the test suite).
#'
#' \strong{Single-sample transform.} Each specimen is mapped to the
#' self-contained per-sample robust CLR over its OWN strictly-positive support
#' (geometric-mean centring on \code{v > 0}); this is exactly invariant to
#' per-specimen scaling and uses no cross-row statistic. Missing universe
#' features are absent from a specimen's rCLR support and contribute no tail
#' term (documented single-sample-safe rule); the per-feature reference ECDF is
#' simply skipped for an absent feature, exactly as if that coordinate carried
#' the neutral aggregate \code{0}.
#'
#' \strong{Combining the two detectors and orientation.} ECOD and COPOD live on
#' different (feature-count-dependent) scales, so each per-specimen aggregate is
#' z-standardized against the frozen control-aggregate mean / sd captured at fit,
#' and the two z-scores are averaged. A single frozen orientation scalar
#' (\eqn{+1} or \eqn{-1}), chosen at fit so that larger = more case-like
#' (training AUC of the combined anomaly score against \code{y_train}; flipped
#' if anomaly is anti-correlated with case), is applied at score. Orientation
#' and standardization are frozen constants, so the score of a specimen depends
#' only on that specimen and the frozen model.
#'
#' Degenerate inputs (no present features, an all-zero specimen, or an empty
#' control reference) return the neutral score \code{0}.
#'
#' @references
#' Li Z, Zhao Y, Botta N, Ionescu C, Hu X. (2020) COPOD: Copula-Based Outlier
#' Detection. \emph{ICDM}.
#'
#' Li Z, Zhao Y, Hu X, Botta N, Ionescu C, Chen G. (2022) ECOD: Unsupervised
#' Outlier Detection Using Empirical Cumulative Distribution Functions.
#' \emph{IEEE TKDE}.
#'
#' Vovk V, Gammerman A, Shafer G. (2005) Algorithmic Learning in a Random World.
#' Springer.
#'
#' @name singlesample-ecod-copod
NULL


# ----------------------------------------------------------------------------
# Single-sample transform and pure-R frozen-ECDF tail aggregates
# ----------------------------------------------------------------------------

# Self-contained per-sample robust CLR of one specimen over its OWN positive
# support. z(c * v) = z(v) for any c > 0 (exact per-sample scale-invariance);
# zero / non-positive parts map to 0 (they are absent from the support). This
# deliberately does NOT reuse a package rclr helper that centres on a fixed
# pseudocount or on cross-sample reference statistics (either would break exact
# scale-invariance / single-sample deployability).
.ecod_copod_rclr <- function(v) {
  z <- rep.int(0, length(v))
  pos <- which(v > 0)
  if (length(pos) >= 1L) {
    vp <- v[pos]
    # log(v / max(v)) - mean(log(v / max(v))) is algebraically identical to
    # log(v) - mean(log(v)), but removes the per-row scale before taking logs.
    lv <- log(vp / max(vp))
    z[pos] <- lv - mean(lv)
  }
  z
}

# Row-wise rCLR of a matrix. rbind (not t(vapply)) keeps the n x p shape even
# when n == 1 or p == 1.
.ecod_copod_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .ecod_copod_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Population skewness m3 / m2^1.5 of a numeric vector, matching the sign of
# scipy.stats.skew (the biased / population estimator pyod's ECOD uses; nan -> 0
# for a constant column). Verified to share the sign of scipy.stats.skew.
.ecod_copod_skew <- function(x) {
  m <- mean(x)
  d <- x - m
  m2 <- mean(d^2)
  if (!is.finite(m2) || m2 <= 0) return(0)
  m3 <- mean(d^3)
  m3 / m2^1.5
}

# Tie tolerance for snapping a query value back to an exact frozen reference
# point after harmless floating-point round-off. The ECDF is a STEP function of
# x, so a round-off-sized perturbation across a reference value can flip one
# rank count. The tolerance is therefore sized to expected arithmetic round-off,
# not to biological/platform noise:
#   tol(x) = .ECOD_COPOD_TIE_ULP * .Machine$double.eps * (abs(x) + 1)
# with ULP = 2048. Crucially, this tolerance is NOT applied by widening the ECDF
# inequalities. PyOD's column_ecdf uses exact maximum-rank ties, and large real
# references can contain genuinely distinct values inside this ULP window. The
# helper below only snaps the query to its nearest frozen reference value; the
# subsequent ECDF count remains exact, preserving pyod parity and avoiding
# near-neighbour rank merging.
.ECOD_COPOD_TIE_ULP <- 2048L

# Snap only the query value, not the reference ECDF, when per-sample scaling
# round-off lands within the ULP window of a frozen reference point. PyOD's
# column_ecdf uses exact maximum-rank ties; widening the inequality itself would
# merge genuinely distinct near-neighbour ranks in large references.
.ecod_copod_snap_to_ref <- function(ref, x) {
  tol <- .ECOD_COPOD_TIE_ULP * .Machine$double.eps * (abs(x) + 1)
  n <- length(ref)
  if (n < 1L) return(x)
  k <- findInterval(x, ref)
  cand <- unique(pmin(pmax(c(k, k + 1L), 1L), n))
  d <- abs(ref[cand] - x)
  best <- which.min(d)
  if (length(best) == 1L && is.finite(d[[best]]) && d[[best]] <= tol) {
    ref[[cand[[best]]]]
  } else {
    x
  }
}

.ecod_copod_snap_to_ref_vec <- function(ref, x) {
  n <- length(ref)
  if (n < 1L || length(x) < 1L) return(x)
  tol <- .ECOD_COPOD_TIE_ULP * .Machine$double.eps * (abs(x) + 1)
  k <- findInterval(x, ref)
  lo <- pmin(pmax(k, 1L), n)
  hi <- pmin(pmax(k + 1L, 1L), n)
  d_lo <- abs(ref[lo] - x)
  d_hi <- abs(ref[hi] - x)
  use_hi <- d_hi < d_lo
  nearest <- ref[lo]
  nearest[use_hi] <- ref[hi[use_hi]]
  d <- d_lo
  d[use_hi] <- d_hi[use_hi]
  snap <- is.finite(d) & d <= tol
  x[snap] <- nearest[snap]
  x
}

# Left-tail empirical CDF P(ref <= x) = (#{ref_i <= x}) / n_ref, evaluated for a
# scalar x against a single frozen reference column `ref`. Matches pyod's
# column_ecdf: tied values all receive the maximum rank probability, i.e. the
# exact count of values <= x over n. For scaled single-row scoring, `x` is first
# snapped to the nearest frozen reference value only when it is within the
# measured floating-point round-off window; the subsequent ECDF count remains an
# exact pyod-style rank count and does not merge neighbouring distinct values.
# The count is floored at 0.5 (strictly below the smallest training-row count of
# 1) BEFORE dividing by n, so: (a) on a training-reference row the true count is
# >= 1 and the floor never bites -> exact match to pyod's decision_scores_; (b) a
# new specimen more extreme than every reference value (count 0) yields a finite
# -log(0.5/n) instead of -log(0) = +Inf. The symmetric right tail P(ref >= x) is
# the same construction on the exact count of values >= snapped x.
.ecod_copod_ecdf_left <- function(ref, x) {
  x <- .ecod_copod_snap_to_ref(ref, x)
  max(findInterval(x, ref), 0.5) / length(ref)
}
.ecod_copod_ecdf_right <- function(ref, x) {
  x <- .ecod_copod_snap_to_ref(ref, x)
  max(length(ref) - findInterval(x, ref, left.open = TRUE), 0.5) / length(ref)
}

# ECOD and COPOD per-specimen aggregates against the frozen per-feature control
# reference columns `ref_cols` (a list, one numeric vector per universe feature
# in `feat_order`) and frozen per-feature skew signs `skew_sign`. `present` are
# the names of the features the specimen actually carries (intersect of universe
# and the scored matrix columns); absent features contribute no tail term (their
# aggregate coordinate is the neutral 0), realising the documented missing-
# feature rule. `z` is the specimen's rCLR aligned to `feat_order`.
#
# ECOD:  O_j = max(U_l, U_r, U_skew); score = sum_j O_j
# COPOD: O_j = max(U_skew, (U_l + U_r) / 2); score = sum_j O_j
# with U_l = -log P(ref <= x), U_r = -log P(ref >= x), and the skew-driven
# U_skew = U_r (skew>0) / U_l (skew<0) / U_l + U_r (skew==0). This is the exact
# pyod ECOD/COPOD aggregation, evaluated against the FROZEN reference.
.ecod_copod_aggregates <- function(z, feat_order, present_idx, ref_cols,
                                   skew_sign) {
  ecod <- 0
  copod <- 0
  for (j in present_idx) {
    ref <- ref_cols[[j]]
    x <- z[j]
    u_l <- -log(.ecod_copod_ecdf_left(ref, x))
    u_r <- -log(.ecod_copod_ecdf_right(ref, x))
    s <- skew_sign[[j]]
    u_skew <- if (s > 0) u_r else if (s < 0) u_l else (u_l + u_r)
    ecod <- ecod + max(u_l, u_r, u_skew)
    copod <- copod + max(u_skew, (u_l + u_r) / 2)
  }
  c(ecod = ecod, copod = copod)
}

.ecod_copod_aggregates_matrix <- function(Z, present_idx, ref_cols,
                                          skew_sign) {
  ecod <- numeric(nrow(Z))
  copod <- numeric(nrow(Z))
  for (j in present_idx) {
    ref <- ref_cols[[j]]
    x <- .ecod_copod_snap_to_ref_vec(ref, Z[, j])
    u_l <- -log(pmax(findInterval(x, ref), 0.5) / length(ref))
    u_r <- -log(pmax(length(ref) - findInterval(x, ref, left.open = TRUE),
                     0.5) / length(ref))
    s <- skew_sign[[j]]
    u_skew <- if (s > 0) u_r else if (s < 0) u_l else (u_l + u_r)
    ecod <- ecod + pmax(u_l, u_r, u_skew)
    copod <- copod + pmax(u_skew, (u_l + u_r) / 2)
  }
  cbind(ecod = ecod, copod = copod)
}

# Fit-time control aggregates against the full control reference. This is the
# vectorized equivalent of `.ecod_copod_aggregates()` for rows that are already
# in `Z_ctrl`: pyod's exact maximum-rank ECDF is just rank(..., ties="max") / n.
.ecod_copod_control_aggregates <- function(Z_ctrl, skew_sign) {
  n <- nrow(Z_ctrl)
  p <- ncol(Z_ctrl)
  # matrix(..., nrow = n) keeps an n x p shape even when n == 1, where vapply()
  # otherwise simplifies the per-feature numeric(1) results to a length-p vector
  # and the matrix-indexed u_skew[, pos] assignment below would error. For n == 1
  # every rank is 1, so u_l = u_r = -log(1) = 0 (matching pyod's single-row ref).
  u_l <- matrix(vapply(seq_len(p), function(j) {
    -log(rank(Z_ctrl[, j], ties.method = "max") / n)
  }, numeric(n)), nrow = n)
  u_r <- matrix(vapply(seq_len(p), function(j) {
    -log(rank(-Z_ctrl[, j], ties.method = "max") / n)
  }, numeric(n)), nrow = n)

  u_skew <- u_l
  pos <- skew_sign > 0
  zero <- skew_sign == 0
  if (any(pos)) u_skew[, pos] <- u_r[, pos, drop = FALSE]
  if (any(zero)) u_skew[, zero] <- u_l[, zero, drop = FALSE] +
    u_r[, zero, drop = FALSE]

  cbind(
    ecod = rowSums(pmax(u_l, u_r, u_skew)),
    copod = rowSums(pmax(u_skew, (u_l + u_r) / 2))
  )
}

# Combined, oriented per-specimen score from the two raw aggregates and the
# frozen standardization / orientation constants. Each aggregate is centred and
# scaled by the frozen control mean/sd (sd floored at eps so a constant control
# aggregate does not divide by zero), the two z-scores averaged, and the frozen
# orientation applied so larger = more case-like.
.ecod_copod_combine <- function(agg, mu, sigma, orientation, eps) {
  z_ecod <- (agg[["ecod"]] - mu[["ecod"]]) / max(sigma[["ecod"]], eps)
  z_copod <- (agg[["copod"]] - mu[["copod"]]) / max(sigma[["copod"]], eps)
  orientation * ((z_ecod + z_copod) / 2)
}


# ----------------------------------------------------------------------------
# Hyperparameter resolver
# ----------------------------------------------------------------------------

.ecod_copod_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_ecod_copod: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_ecod_copod: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_ecod_copod: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("min_features", "eps", "verify_pyod")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_ecod_copod: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_ecod_copod: hp$min_features must be a positive integer")
  }

  eps <- hp$eps
  if (is.null(eps)) eps <- 1e-8
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_ecod_copod: hp$eps must be a positive finite number")
  }

  verify_pyod <- hp$verify_pyod
  if (is.null(verify_pyod)) verify_pyod <- TRUE
  if (!is.logical(verify_pyod) || length(verify_pyod) != 1L ||
      is.na(verify_pyod)) {
    stop("fit_ecod_copod: hp$verify_pyod must be a single logical")
  }

  list(
    min_features = as.integer(min_features),
    eps = as.numeric(eps),
    verify_pyod = verify_pyod
  )
}


# ----------------------------------------------------------------------------
# Python availability guard (mirrors .tabdl_require_python_module)
# ----------------------------------------------------------------------------

.ecod_copod_require_pyod <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_ecod_copod: package 'reticulate' is required to fit the ECOD/COPOD ",
         "detectors. Install with install.packages('reticulate').",
         call. = FALSE)
  }
  if (!reticulate::py_module_available("pyod")) {
    stop("fit_ecod_copod: Python module 'pyod' is not available in the active ",
         "reticulate environment. Install with: ",
         "reticulate::py_install('pyod') (or pip install pyod).",
         call. = FALSE)
  }
  invisible(TRUE)
}


# ----------------------------------------------------------------------------
# fit_ecod_copod
# ----------------------------------------------------------------------------

#' @title Fit the ECOD + COPOD novelty discriminator
#'
#' @description
#' Fits the bundled ECOD + COPOD novelty discriminator. The training matrix is
#' mapped to the per-sample robust CLR (over each row's own positive support),
#' and the in-distribution reference is taken to be the training CONTROL rows
#' (\code{y_train == 0}). \pkg{pyod}'s \code{ECOD} and \code{COPOD} are fitted on
#' that control rCLR (leakage-safe: training only) and used to (a) source /
#' verify the per-feature skewness signs and (b) provide the
#' \code{decision_scores_} anchor that the frozen pure-R reimplementation must
#' reproduce. The frozen model stores only R-side state -- the per-feature
#' control reference columns, the skew signs, the per-detector control-aggregate
#' mean / sd, and the orientation scalar -- so scoring is pure R and never calls
#' Python.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and otherwise ignored.
#' @param hp Optional list of hyperparameters. Allowed fields:
#'   \code{min_features} (positive integer feature-overlap floor at scoring;
#'   default \code{3L}), \code{eps} (positive finite floor for the
#'   standardization sd and a numeric guard; default \code{1e-8}), and
#'   \code{verify_pyod} (single logical; when \code{TRUE} (default) the frozen
#'   pure-R reimplementation is checked against \pkg{pyod}'s
#'   \code{decision_scores_} on the control rows at fit and the fit fails if the
#'   maximum difference exceeds \code{1e-6}).
#'
#' @return Object of class \code{ecod_copod_model}: a list with
#'   \code{feature_universe}, \code{ref_cols} (frozen per-feature control
#'   reference columns), \code{skew_sign}, \code{score_mean} / \code{score_sd}
#'   (per-detector control-aggregate standardization constants),
#'   \code{orientation}, \code{n_control}, and \code{hp}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120; p <- 30; k <- 8
#' L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
#' X <- exp(L)
#' model <- fit_ecod_copod(X, y)
#' score_ecod_copod(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Li Z, Zhao Y, Hu X, Botta N, Ionescu C, Chen G. (2022) ECOD. \emph{IEEE TKDE}.
#'
#' Li Z, Zhao Y, Botta N, Ionescu C, Hu X. (2020) COPOD. \emph{ICDM}.
#'
#' @export
fit_ecod_copod <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_ecod_copod", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_ecod_copod", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_ecod_copod")
  hp <- .ecod_copod_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_ecod_copod: X_train must contain at least hp$min_features features")
  }
  .ecod_copod_require_pyod()

  feat_order <- colnames(X_train)
  Z <- .ecod_copod_rclr_matrix(X_train)            # n x p, full universe
  control_idx <- which(y == 0L)
  Z_ctrl <- Z[control_idx, , drop = FALSE]         # in-distribution reference

  # Per-feature skew sign from scipy via the fitted pyod ECOD object, recomputed
  # in pure R and cross-checked so the frozen signs match pyod's exactly.
  ecod_mod <- reticulate::import("pyod.models.ecod", delay_load = FALSE)
  copod_mod <- reticulate::import("pyod.models.copod", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)
  Z_ctrl_py <- np$asarray(Z_ctrl, dtype = "float64")
  ecod_fit <- ecod_mod$ECOD()
  ecod_fit$fit(Z_ctrl_py)
  copod_fit <- copod_mod$COPOD()
  copod_fit$fit(Z_ctrl_py)

  skew_sign <- vapply(seq_len(ncol(Z_ctrl)), function(j) {
    sign(.ecod_copod_skew(Z_ctrl[, j]))
  }, numeric(1L))

  # Frozen sorted per-feature reference columns (one numeric vector per universe
  # feature). Stored as a list so the score path indexes by feature position.
  ref_cols <- lapply(seq_len(ncol(Z_ctrl)), function(j) sort(Z_ctrl[, j]))
  names(ref_cols) <- feat_order
  names(skew_sign) <- feat_order

  # Frozen pure-R control aggregates -> standardization constants. Every universe
  # feature is present for a control row, so present_idx is the full set here.
  full_idx <- seq_len(ncol(Z_ctrl))
  ctrl_agg <- .ecod_copod_control_aggregates(Z_ctrl, skew_sign)
  score_mean <- c(ecod = mean(ctrl_agg[, "ecod"]),
                  copod = mean(ctrl_agg[, "copod"]))
  score_sd <- c(ecod = stats::sd(ctrl_agg[, "ecod"]),
                copod = stats::sd(ctrl_agg[, "copod"]))
  if (!is.finite(score_sd[["ecod"]])) score_sd["ecod"] <- 0
  if (!is.finite(score_sd[["copod"]])) score_sd["copod"] <- 0

  # Optional §7 fit-time anchor: the frozen pure-R aggregates on the control
  # rows must equal pyod's decision_scores_ (same reference, same convention).
  if (hp$verify_pyod) {
    ds_ecod <- as.numeric(ecod_fit$decision_scores_)
    ds_copod <- as.numeric(copod_fit$decision_scores_)
    d_ecod <- max(abs(ctrl_agg[, "ecod"] - ds_ecod))
    d_copod <- max(abs(ctrl_agg[, "copod"] - ds_copod))
    if (d_ecod > 1e-6 || d_copod > 1e-6) {
      stop(sprintf(paste0("fit_ecod_copod: frozen pure-R reimplementation does ",
                          "not match pyod decision_scores_ (ECOD maxdiff %.3g, ",
                          "COPOD maxdiff %.3g); the ECDF/skew convention is wrong."),
                   d_ecod, d_copod), call. = FALSE)
    }
  }

  # Orientation: combined anomaly z-score on the TRAINING rows, oriented so the
  # combined score is positively associated with case. Computed with a +1
  # orientation, then flipped if anomaly is anti-correlated with case (training
  # AUC < 0.5). Uses the frozen standardization constants.
  train_agg <- .ecod_copod_aggregates_matrix(Z, full_idx, ref_cols, skew_sign)
  train_combined <- ((train_agg[, "ecod"] - score_mean[["ecod"]]) /
                       max(score_sd[["ecod"]], hp$eps) +
                       (train_agg[, "copod"] - score_mean[["copod"]]) /
                       max(score_sd[["copod"]], hp$eps)) / 2
  auc_pos <- .ecod_copod_auc(y, train_combined)
  orientation <- if (is.finite(auc_pos) && auc_pos < 0.5) -1 else 1

  model <- list(
    feature_universe = feat_order,
    ref_cols = ref_cols,
    skew_sign = skew_sign,
    score_mean = score_mean,
    score_sd = score_sd,
    orientation = orientation,
    n_control = nrow(Z_ctrl),
    hp = hp
  )
  class(model) <- "ecod_copod_model"
  model
}

# Mann-Whitney AUC of `score` against 0/1 labels `y` (P(score_case > score_ctrl)
# with ties at 0.5). Self-contained so the fit does not depend on pROC.
.ecod_copod_auc <- function(y, score) {
  case <- score[y == 1L]
  ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}


# ----------------------------------------------------------------------------
# score_ecod_copod
# ----------------------------------------------------------------------------

#' @title Score the ECOD + COPOD novelty discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_ecod_copod}}. For one specimen the abundances are mapped to
#' the per-sample robust CLR over its own positive support; its ECOD and COPOD
#' tail aggregates are computed in pure R against the frozen per-feature control
#' reference ECDFs (never \pkg{pyod}); each aggregate is z-standardized with the
#' frozen control mean / sd, the two z-scores are averaged, and the frozen
#' orientation is applied so larger = more case-like.
#'
#' The present feature set is \code{intersect(model$feature_universe,
#' colnames(X))}; features missing from a specimen contribute no tail term. If
#' fewer than \code{model$hp$min_features} universe features are present, the
#' neutral score \code{0} is returned for every row. The score of a row depends
#' only on that row and the frozen model, and is exactly invariant to
#' per-specimen scaling.
#'
#' @param model An \code{ecod_copod_model} object from
#'   \code{\link{fit_ecod_copod}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with named feature columns.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}; larger values
#'   are more case-like.
#'
#' @examples
#' \dontrun{
#' model <- fit_ecod_copod(X, y)
#' score_ecod_copod(model, X[1, , drop = FALSE])
#' }
#'
#' @export
score_ecod_copod <- function(model, X, meta = NULL) {
  if (!inherits(model, "ecod_copod_model")) {
    stop("score_ecod_copod: model must have class ecod_copod_model")
  }
  X <- .reo_check_matrix(X, "score_ecod_copod", "X")
  .reo_check_meta(meta, nrow(X), "score_ecod_copod", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))
  if (length(present) < model$hp$min_features || model$n_control < 1L) {
    return(rep(0, nrow(X)))
  }
  present_idx <- match(present, feat_order)

  # Align X to the frozen universe order; absent features stay 0 (no support).
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  X_use[, present] <- X[, present, drop = FALSE]

  mu <- model$score_mean
  sd <- model$score_sd
  orientation <- model$orientation
  eps <- model$hp$eps
  ref_cols <- model$ref_cols
  skew_sign <- model$skew_sign

  Z <- .ecod_copod_rclr_matrix(X_use)
  empty <- rowSums(Z != 0) == 0                      # all-zero specimen: neutral
  agg <- .ecod_copod_aggregates_matrix(Z, present_idx, ref_cols, skew_sign)
  out <- orientation * (((agg[, "ecod"] - mu[["ecod"]]) /
                           max(sd[["ecod"]], eps) +
                           (agg[, "copod"] - mu[["copod"]]) /
                           max(sd[["copod"]], eps)) / 2)
  out[empty] <- 0

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_ecod_copod: scorer produced non-finite or wrong-length output")
  }
  out
}


# ----------------------------------------------------------------------------
# Stable model digest for the equivariance harness
# ----------------------------------------------------------------------------

# Hash ONLY the frozen R-side state of the model (feature order, reference
# columns, skew signs, standardization and orientation constants, n_control,
# hp). The model holds no live Python / external-pointer object, so this is a
# deterministic snapshot suitable as the `model_digest` for
# singlesample_assert_row_equivariant() (clause (d): the bytes must be identical
# before and after scoring).
.ecod_copod_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe,
    ref_cols = model$ref_cols,
    skew_sign = model$skew_sign,
    score_mean = model$score_mean,
    score_sd = model$score_sd,
    orientation = model$orientation,
    n_control = model$n_control,
    hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
