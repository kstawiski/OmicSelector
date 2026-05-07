#' @title Paper 3 outlier detection and conformal claim-gating (Module D)
#'
#' @description
#' Outlier detection and conformal anomaly-scoring primitives introduced in
#' Paper 3 (Module D; Stawiski et al., in preparation). These methods flag
#' out-of-distribution samples before a biomarker claim is reported, providing
#' either distribution-free FPR guarantees (conformal approach) or
#' robust compositional distance metrics (Mahalanobis / isolation forest).
#'
#' Methods provided:
#' \itemize{
#'   \item \code{\link{fit_compositional_mahalanobis}} /
#'     \code{\link{apply_compositional_mahalanobis}}: MCD-robust Mahalanobis
#'     distance on ILR or rCLR log-ratio coordinates.
#'   \item \code{\link{fit_conformal_anomaly}} /
#'     \code{\link{os_conformal_anomaly}}: conformal p-value via k-NN distance
#'     against a held-out calibration partition.
#'   \item \code{\link{fit_isolation_forest_logratio}} /
#'     \code{\link{apply_isolation_forest_logratio}}: pure-R isolation forest
#'     on rCLR-transformed inputs.
#' }
#'
#' @references
#' Filzmoser P, Hron K, Reimann C. (2009) Principal component analysis for
#' compositional data with outliers. \emph{Environmetrics} 20: 621-632.
#'
#' Rousseeuw PJ, Van Driessen K. (1999) A Fast Algorithm for the Minimum
#' Covariance Determinant Estimator. \emph{Technometrics} 41(3): 212-223.
#'
#' Vovk V, Gammerman A, Shafer G. (2005) Algorithmic Learning in a Random
#' World. Springer.
#'
#' Liu FT, Ting KM, Zhou Z-H. (2008) Isolation Forest. \emph{ICDM}.
#'
#' Stawiski K. (in preparation) Paper 3 of the OmicSelector programme.
#'
#' @name paper3-outlier-detection
NULL


# ============================================================================
# Compositional Mahalanobis (MCD)
# ============================================================================

# --- internal transform helpers ---------------------------------------------

#' @noRd
.paper3_compose_transform <- function(x, transform = "rclr", pseudocount = 0.5) {
  x <- x + pseudocount
  log_x <- log(x)
  if (transform == "rclr") {
    centering <- rowMeans(log_x, na.rm = TRUE)
    return(sweep(log_x, 1L, centering, FUN = "-"))
  } else if (transform == "ilr") {
    D <- ncol(log_x)
    H <- .paper3_helmert_basis(D)
    return(log_x %*% H)
  }
  stop("unknown transform: ", transform)
}

#' @noRd
.paper3_helmert_basis <- function(D) {
  H <- matrix(0, D, D - 1L)
  for (k in seq_len(D - 1L)) {
    s <- sqrt(k / (k + 1))
    H[seq_len(k), k] <- s / k
    H[k + 1L, k] <- -s
  }
  H
}


# ----------------------------------------------------------------------------
# fit_compositional_mahalanobis
# ----------------------------------------------------------------------------

#' @title Fit robust Mahalanobis detector on compositional log-ratio coordinates
#'
#' @description
#' Fits a Minimum Covariance Determinant (MCD) robust Mahalanobis distance
#' model on ILR (default) or rCLR log-ratio coordinates of a training
#' compositional matrix. The MCD estimator is provided by the
#' \code{robustbase} package. The function fails closed if \code{robustbase}
#' is unavailable and \code{require_robust = TRUE} (the default).
#'
#' @param x_train Numeric matrix (samples \eqn{\times} features). All values
#'   must be non-negative.
#' @param transform \code{"ilr"} (default, \eqn{D \to D-1} dimensional, full
#'   rank under healthy assumptions) or \code{"rclr"} (sum-to-zero constraint,
#'   effective rank \eqn{D-1}).
#' @param alpha Numeric in (0.5, 1). MCD coverage fraction. Default 0.75.
#' @param pseudocount Numeric \eqn{> 0}. Added before log transform. Default
#'   0.5.
#' @param require_robust Logical. If \code{TRUE} (default), stops when
#'   \code{robustbase} is unavailable rather than silently downgrading to
#'   classical covariance.
#'
#' @return Object of class \code{compositional_mahalanobis_fit}.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x_train <- matrix(abs(rnorm(80 * 20, 100, 30)), nrow = 80,
#'                   dimnames = list(NULL, paste0("miR-", seq_len(20))))
#' fit <- fit_compositional_mahalanobis(x_train, transform = "rclr",
#'                                      require_robust = FALSE)
#' }
#'
#' @references
#' Rousseeuw PJ, Van Driessen K. (1999) \emph{Technometrics} 41: 212-223.
#'
#' @export
fit_compositional_mahalanobis <- function(x_train,
                                           transform = c("ilr", "rclr"),
                                           alpha = 0.75,
                                           pseudocount = 0.5,
                                           require_robust = TRUE) {
  transform <- match.arg(transform)
  if (!is.matrix(x_train)) x_train <- as.matrix(x_train)
  if (any(x_train < 0, na.rm = TRUE)) {
    stop("fit_compositional_mahalanobis: negative values not allowed (compositional input)")
  }
  if (pseudocount <= 0) stop("fit_compositional_mahalanobis: pseudocount must be > 0")

  z_train <- .paper3_compose_transform(x_train, transform = transform,
                                        pseudocount = pseudocount)
  D <- ncol(z_train)
  effective_df <- if (transform == "ilr") D else max(1L, D - 1L)

  rank_cov <- min(nrow(z_train) - 1L, effective_df)
  if (rank_cov < effective_df) {
    warning("fit_compositional_mahalanobis: covariance is rank-deficient (n=",
             nrow(z_train), ", df=", effective_df,
             "). Chi-square p-values will be conservative.")
  }

  has_robustbase <- requireNamespace("robustbase", quietly = TRUE)
  has_MASS <- requireNamespace("MASS", quietly = TRUE)

  if (has_robustbase) {
    mcd <- robustbase::covMcd(z_train, alpha = alpha)
    center <- mcd$center
    cov_mat <- mcd$cov
    method <- "MCD"
  } else {
    if (require_robust) {
      stop("fit_compositional_mahalanobis: robustbase unavailable and ",
            "require_robust=TRUE. Install robustbase or pass require_robust=FALSE.")
    }
    warning("robustbase unavailable - falling back to classical covariance.")
    center <- colMeans(z_train, na.rm = TRUE)
    cov_mat <- stats::cov(z_train, use = "pairwise.complete.obs")
    method <- "classical"
  }

  if (!has_MASS) stop("fit_compositional_mahalanobis requires the MASS package for ginv().")
  cov_inv <- MASS::ginv(cov_mat)

  fit <- list(
    transform = transform,
    pseudocount = pseudocount,
    feature_names = colnames(x_train),
    center = center,
    cov = cov_mat,
    cov_inv = cov_inv,
    method = method,
    df = effective_df,
    n_train = nrow(z_train)
  )
  class(fit) <- "compositional_mahalanobis_fit"
  fit
}


#' @title Apply compositional Mahalanobis distance to new samples
#'
#' @description
#' Projects new samples into the same log-ratio space and computes their
#' squared Mahalanobis distance using the frozen center and inverse-covariance
#' from \code{\link{fit_compositional_mahalanobis}}.
#'
#' @param fit A \code{compositional_mahalanobis_fit} object.
#' @param x_test Numeric matrix (samples \eqn{\times} features). All values
#'   must be non-negative. All training features must be present.
#' @param return_pvalue Logical. If \code{TRUE} (default), also returns
#'   chi-square reference p-values.
#'
#' @return List with \code{distance_sq}, \code{distance}, and (if
#'   \code{return_pvalue}) \code{pvalue_chisq}.
#'
#' @export
apply_compositional_mahalanobis <- function(fit, x_test, return_pvalue = TRUE) {
  if (!inherits(fit, "compositional_mahalanobis_fit")) {
    stop("fit must be class compositional_mahalanobis_fit")
  }
  if (!is.matrix(x_test)) x_test <- as.matrix(x_test)
  if (any(x_test < 0, na.rm = TRUE)) {
    stop("apply_compositional_mahalanobis: negative values not allowed")
  }
  missing <- setdiff(fit$feature_names, colnames(x_test))
  if (length(missing) > 0L) {
    stop("apply_compositional_mahalanobis: ", length(missing),
          " training features missing: ",
          paste(utils::head(missing, 5L), collapse = ", "), " ...")
  }
  x_test <- x_test[, fit$feature_names, drop = FALSE]
  z_test <- .paper3_compose_transform(x_test, transform = fit$transform,
                                       pseudocount = fit$pseudocount)
  delta <- sweep(z_test, 2L, fit$center, FUN = "-")
  d2 <- rowSums((delta %*% fit$cov_inv) * delta)

  out <- list(distance_sq = d2, distance = sqrt(pmax(0, d2)))
  if (return_pvalue) {
    out$pvalue_chisq <- stats::pchisq(d2, df = fit$df, lower.tail = FALSE)
  }
  out
}


# ============================================================================
# Conformal anomaly detector
# ============================================================================

#' @title Fit a conformal anomaly detector from a Tier R healthy reference cohort
#'
#' @description
#' Fits a conformal anomaly detector using k-nearest-neighbour distance as the
#' non-conformity measure. A held-out calibration fraction of the healthy
#' reference cohort is used to calibrate the conformal p-value, providing a
#' distribution-free, finite-sample FPR guarantee at level \code{alpha}
#' (Vovk et al., 2005; Romano et al., 2019).
#'
#' @param X_ref Numeric matrix (samples \eqn{\times} features) of healthy
#'   reference samples.
#' @param representation Function to project samples to a latent space. Default
#'   \code{identity} (use raw features).
#' @param k Integer - number of nearest neighbours. Default 10.
#' @param calibration_split Fraction of \code{X_ref} held out for calibration.
#'   Default 0.30.
#' @param seed Integer - RNG seed for the train/calibration split. Default 42.
#'
#' @return Object of class \code{conformal_anomaly_fit}.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X_ref <- matrix(rnorm(200 * 30), nrow = 200)
#' fit <- fit_conformal_anomaly(X_ref, k = 10L)
#' fit$n_cal
#' }
#'
#' @references
#' Vovk V, Gammerman A, Shafer G. (2005) Algorithmic Learning in a Random
#' World. Springer.
#'
#' @export
fit_conformal_anomaly <- function(X_ref,
                                   representation = identity,
                                   k = 10L,
                                   calibration_split = 0.30,
                                   seed = 42L) {
  X <- as.matrix(X_ref)
  X_repr <- representation(X)
  if (is.vector(X_repr)) X_repr <- matrix(X_repr, nrow = nrow(X))
  X_repr <- as.matrix(X_repr)

  set.seed(seed)
  n <- nrow(X_repr)
  cal_idx <- sample(seq_len(n), size = floor(calibration_split * n))
  ref_idx <- setdiff(seq_len(n), cal_idx)
  X_ref_manifold <- X_repr[ref_idx, , drop = FALSE]
  X_calibration  <- X_repr[cal_idx, , drop = FALSE]

  cal_scores <- apply(X_calibration, 1L, function(s) {
    dists <- sqrt(rowSums(sweep(X_ref_manifold, 2L, s, "-")^2))
    sort(dists)[k]
  })

  fit <- list(
    representation = representation,
    X_ref_manifold = X_ref_manifold,
    cal_scores = cal_scores,
    k = k,
    n_ref = nrow(X_ref_manifold),
    n_cal = length(cal_scores)
  )
  class(fit) <- "conformal_anomaly_fit"
  fit
}


#' @title Compute conformal anomaly p-value for new samples
#'
#' @description
#' Scores one or more new samples against the calibrated conformal anomaly
#' model from \code{\link{fit_conformal_anomaly}}. The conformal p-value is
#' the fraction of calibration scores at least as extreme as the test score:
#' \deqn{p = \frac{1 + \#\{s_{\mathrm{cal}} \geq s_{\mathrm{test}}\}}{1 + n_{\mathrm{cal}}}.}
#'
#' @param x Numeric vector or matrix (samples \eqn{\times} features).
#' @param fit A \code{conformal_anomaly_fit} object from
#'   \code{\link{fit_conformal_anomaly}}.
#' @param alpha Significance level for the binary anomaly call. Default 0.05.
#'
#' @return List with \code{score} (raw k-NN distance), \code{p_value}
#'   (conformal p in \eqn{[0,1]}), \code{is_anomaly} (logical at level
#'   \code{alpha}), \code{alpha}, and \code{calibration_summary}.
#'
#' @export
os_conformal_anomaly <- function(x, fit, alpha = 0.05) {
  stopifnot(inherits(fit, "conformal_anomaly_fit"))
  if (is.vector(x)) {
    x_mat <- matrix(x, nrow = 1L, dimnames = list(NULL, names(x)))
  } else {
    x_mat <- as.matrix(x)
  }
  x_repr <- fit$representation(x_mat)
  if (is.vector(x_repr)) x_repr <- matrix(x_repr, nrow = nrow(x_mat))
  x_repr <- as.matrix(x_repr)

  scores <- apply(x_repr, 1L, function(s) {
    dists <- sqrt(rowSums(sweep(fit$X_ref_manifold, 2L, s, "-")^2))
    sort(dists)[fit$k]
  })

  p_values <- vapply(scores, function(s) {
    (1 + sum(fit$cal_scores >= s)) / (1 + fit$n_cal)
  }, numeric(1L))

  is_anomaly <- p_values < alpha

  list(
    score = scores,
    p_value = p_values,
    is_anomaly = is_anomaly,
    alpha = alpha,
    calibration_summary = sprintf(
      "conformal_anomaly (k=%d, n_ref=%d, n_cal=%d, alpha=%.3f)",
      fit$k, fit$n_ref, fit$n_cal, alpha)
  )
}


# ============================================================================
# Isolation Forest on rCLR log-ratio inputs
# ============================================================================

# --- pure-R isolation tree helpers ------------------------------------------

#' @noRd
.paper3_itree_split <- function(x, max_depth = 8L, depth = 0L) {
  n <- nrow(x)
  if (n <= 1L || depth >= max_depth) {
    return(list(type = "leaf", n = n, depth = depth))
  }
  feat <- sample(ncol(x), 1L)
  vals <- x[, feat]
  if (length(unique(vals)) < 2L) {
    return(list(type = "leaf", n = n, depth = depth))
  }
  split_val <- stats::runif(1L, min(vals), max(vals))
  left  <- x[vals <  split_val, , drop = FALSE]
  right <- x[vals >= split_val, , drop = FALSE]
  list(
    type  = "split",
    feat  = feat,
    split = split_val,
    left  = .paper3_itree_split(left,  max_depth, depth + 1L),
    right = .paper3_itree_split(right, max_depth, depth + 1L)
  )
}

#' @noRd
.paper3_itree_path_length <- function(tree, x) {
  if (tree$type == "leaf") return(tree$depth + .paper3_c_factor(tree$n))
  if (x[tree$feat] < tree$split) .paper3_itree_path_length(tree$left,  x)
  else                            .paper3_itree_path_length(tree$right, x)
}

#' @noRd
.paper3_c_factor <- function(n) {
  if (n <= 1L) return(0)
  2 * (log(n - 1L) + 0.5772156649) - 2 * (n - 1L) / n
}

#' @noRd
.paper3_rclr <- function(x, pseudocount = 0.5) {
  x <- x + pseudocount
  log_x <- log(x)
  centering <- rowMeans(log_x, na.rm = TRUE)
  sweep(log_x, 1L, centering, FUN = "-")
}


# ----------------------------------------------------------------------------
# fit_isolation_forest_logratio
# ----------------------------------------------------------------------------

#' @title Fit a pure-R isolation forest on rCLR log-ratio inputs
#'
#' @description
#' Builds a pure-R isolation forest (Liu et al., 2008) on rCLR-transformed
#' compositional inputs. Anomalies have shorter expected path lengths to
#' isolation-tree leaves. The decision threshold is calibrated on training
#' scores to control the empirical false-positive rate at \code{fpr_target}
#' (applied at scoring time in \code{\link{apply_isolation_forest_logratio}}).
#'
#' This implementation uses no external dependencies and is intended for
#' robustness and auditability. For production use on large panels, a C++
#' backend (e.g., \code{isotree}) would be faster.
#'
#' @param x_train Numeric matrix (samples \eqn{\times} features).
#' @param n_trees Integer - number of isolation trees. Default 100.
#' @param sample_size Integer - subsample size per tree. Default 256.
#' @param max_depth Integer or \code{NULL}. Maximum tree depth. Default
#'   \code{ceiling(log2(min(sample_size, n)))}.
#' @param pseudocount Numeric. rCLR pseudocount. Default 0.5.
#' @param seed Integer - RNG seed. Default 42.
#'
#' @return Object of class \code{isolation_forest_logratio_fit}.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x_train <- matrix(abs(rnorm(80 * 20, 100, 30)), nrow = 80,
#'                   dimnames = list(NULL, paste0("miR-", seq_len(20))))
#' fit <- fit_isolation_forest_logratio(x_train, n_trees = 20L, sample_size = 32L)
#' }
#'
#' @references
#' Liu FT, Ting KM, Zhou Z-H. (2008) Isolation Forest. \emph{ICDM}.
#'
#' @export
fit_isolation_forest_logratio <- function(x_train,
                                           n_trees = 100L,
                                           sample_size = 256L,
                                           max_depth = NULL,
                                           pseudocount = 0.5,
                                           seed = 42L) {
  if (!is.matrix(x_train)) x_train <- as.matrix(x_train)
  z_train <- .paper3_rclr(x_train, pseudocount = pseudocount)
  n <- nrow(z_train)
  set.seed(seed)
  if (is.null(max_depth)) max_depth <- ceiling(log2(min(sample_size, n)))
  c_n <- .paper3_c_factor(min(sample_size, n))

  forest <- vector("list", n_trees)
  for (t in seq_len(n_trees)) {
    idx <- sample(n, min(sample_size, n), replace = FALSE)
    forest[[t]] <- .paper3_itree_split(z_train[idx, , drop = FALSE],
                                        max_depth = max_depth)
  }

  train_scores <- vapply(seq_len(n), function(i) {
    paths <- vapply(forest,
                    function(tr) .paper3_itree_path_length(tr, z_train[i, ]),
                    numeric(1L))
    2 ^ (-mean(paths) / c_n)
  }, numeric(1L))

  fit <- list(
    forest = forest,
    feature_names = colnames(x_train),
    pseudocount = pseudocount,
    c_n = c_n,
    n_trees = n_trees,
    sample_size = sample_size,
    train_scores = train_scores
  )
  class(fit) <- "isolation_forest_logratio_fit"
  fit
}


#' @title Apply isolation forest to new samples and flag anomalies
#'
#' @description
#' Scores new samples using a fitted isolation forest
#' (\code{\link{fit_isolation_forest_logratio}}), computes empirical p-values
#' relative to the training-score distribution, and flags samples exceeding
#' the training \eqn{(1 - \texttt{fpr\_target})} quantile.
#'
#' @param fit An \code{isolation_forest_logratio_fit} object.
#' @param x_test Numeric matrix (samples \eqn{\times} features). All training
#'   features must be present.
#' @param fpr_target Numeric in (0, 1). Empirical FPR target for the flagging
#'   threshold. Default 0.05.
#'
#' @return List with \code{score} (anomaly score in \eqn{[0,1]}; higher =
#'   more anomalous), \code{pvalue_emp}, \code{threshold}, \code{flagged}.
#'
#' @export
apply_isolation_forest_logratio <- function(fit, x_test, fpr_target = 0.05) {
  if (!inherits(fit, "isolation_forest_logratio_fit")) {
    stop("fit must be class isolation_forest_logratio_fit")
  }
  if (!is.matrix(x_test)) x_test <- as.matrix(x_test)
  shared <- intersect(fit$feature_names, colnames(x_test))
  if (length(shared) == 0L) {
    stop("apply_isolation_forest_logratio: no shared features between fit and x_test")
  }
  x_test <- x_test[, fit$feature_names, drop = FALSE]
  z_test <- .paper3_rclr(x_test, pseudocount = fit$pseudocount)

  test_scores <- vapply(seq_len(nrow(z_test)), function(i) {
    paths <- vapply(fit$forest,
                    function(tr) .paper3_itree_path_length(tr, z_test[i, ]),
                    numeric(1L))
    2 ^ (-mean(paths) / fit$c_n)
  }, numeric(1L))

  p_emp <- vapply(test_scores, function(s) {
    (1 + sum(fit$train_scores >= s)) / (1 + length(fit$train_scores))
  }, numeric(1L))

  threshold <- stats::quantile(fit$train_scores, probs = 1 - fpr_target)
  flagged <- test_scores > threshold

  list(score = test_scores, pvalue_emp = p_emp,
       threshold = threshold, flagged = flagged)
}
