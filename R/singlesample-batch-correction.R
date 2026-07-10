#' @title Single-sample batch-correction methods (Module C)
#'
#' @description
#' Frozen batch-correction primitives introduced in the single-sample scoring
#' bank (Module C). All methods follow a fit-then-freeze
#' paradigm: factor loadings or basis vectors are estimated once on training
#' data, then applied to test samples without re-estimation, making them
#' suitable for clinical deployment where a population-matched training cohort
#' may not be available at prediction time.
#'
#' Methods provided:
#' \itemize{
#'   \item \code{\link{fit_frozen_ruv}} / \code{\link{apply_frozen_ruv}}:
#'     frozen RUV (Remove Unwanted Variation) using empirical negative-control
#'     features. Reference: Risso et al. (Nat Biotechnol 2014).
#'   \item \code{\link{fit_robust_pca_residual}} /
#'     \code{\link{apply_robust_pca_residual}}: MAD-scaled SVD residual filter.
#'     NOT the Candès RPCA (ALM sparse-low-rank) decomposition.
#' }
#'
#' @references
#' Risso D, Ngai J, Speed TP, Dudoit S. (2014) Normalization of RNA-seq data
#' using factor analysis of control genes or samples.
#' \emph{Nature Biotechnology} 32(9): 896–902.
#'
#' Hubert M, Rousseeuw PJ, Vanden Branden K. (2005) ROBPCA: A New Approach
#' to Robust Principal Component Analysis.
#' \emph{Technometrics} 47(1): 64–79.
#'
#' @name singlesample-batch-correction
NULL


# ----------------------------------------------------------------------------
# fit_frozen_ruv / apply_frozen_ruv
# ----------------------------------------------------------------------------

#' @title Fit a frozen RUV factor model from training data
#'
#' @description
#' Estimates latent factors of unwanted variation from empirical
#' negative-control features on training data, then freezes the factor
#' loadings for deployment. This is the frozen-mode extension of RUVSeq
#' (Risso et al., Nat Biotechnol 2014): for a new sample \eqn{x}, the
#' unwanted variation component is estimated via the frozen pseudo-inverse of
#' training latent factors, then subtracted.
#'
#' If fewer than 5 named control features are found in \code{colnames(X_train)},
#' falls back to the top-30\% lowest-variance features as proxy controls, with a
#' warning.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features). Training
#'   cohort.
#' @param control_features Character vector of feature names expected to be
#'   non-cancer-relevant (empirical negative controls). Default = a curated set
#'   of stably-expressed serum/plasma miRNAs.
#' @param k Integer or integer vector (range). Number of unwanted factors to
#'   estimate. Default \code{2:6} — BIC selects within this range.
#' @param select Method for selecting \code{k}. \code{"BIC"} (default),
#'   \code{"fixed"} (uses the smallest value in \code{k}), or \code{"explicit"}
#'   (if \code{select} is a single numeric).
#'
#' @return Object of class \code{frozen_ruv_fit} with components:
#'   \code{feat_mean}, \code{full_loadings}, \code{latent_train_inv},
#'   \code{k_selected}, \code{bics}, \code{n_train}, \code{p_full},
#'   \code{n_controls}.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X <- matrix(rnorm(200 * 50), nrow = 200,
#'             dimnames = list(NULL, paste0("mir-", seq_len(50))))
#' fit <- fit_frozen_ruv(X, k = 2:4)
#' fit$k_selected
#' }
#'
#' @references
#' Risso D, Ngai J, Speed TP, Dudoit S. (2014) Normalization of RNA-seq data
#' using factor analysis of control genes or samples.
#' \emph{Nature Biotechnology} 32(9): 896–902.
#'
#' @export
fit_frozen_ruv <- function(X_train,
                            control_features = NULL,
                            k = 2:6,
                            select = c("BIC", "fixed", "explicit")) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("fit_frozen_ruv requires the MASS package.")
  }
  if (is.numeric(select) && length(select) == 1L) {
    k_target <- as.integer(select)
    select <- "explicit"
  } else {
    select <- match.arg(select)
  }
  X <- as.matrix(X_train)

  if (is.null(control_features)) {
    control_features <- c(
      "hsa-miR-93-5p", "hsa-miR-26a-5p", "hsa-miR-191-5p",
      "hsa-miR-92a-3p", "hsa-miR-30c-5p", "hsa-let-7g-5p",
      "hsa-miR-103a-3p", "hsa-miR-23a-3p", "hsa-miR-21-5p"
    )
  }
  ctrl_idx <- which(colnames(X) %in% control_features)
  if (length(ctrl_idx) < 5L) {
    warning("fit_frozen_ruv: only ", length(ctrl_idx),
            " control features found; using top-30% lowest-variance features as proxy.")
    feat_var <- apply(X, 2L, var, na.rm = TRUE)
    ctrl_idx <- order(feat_var)[seq_len(floor(0.30 * ncol(X)))]
  }

  feat_mean <- colMeans(X, na.rm = TRUE)
  X_centered <- sweep(X, 2L, feat_mean, "-")
  X_ctrl <- X_centered[, ctrl_idx, drop = FALSE]
  svd_res <- svd(X_ctrl)

  n <- nrow(X)
  p_ctrl <- ncol(X_ctrl)

  k_candidates <- if (select == "BIC") k
                  else if (select == "fixed") min(k)
                  else k_target

  log_lik <- function(k_val) {
    var_total <- sum(svd_res$d^2)
    var_top   <- sum(svd_res$d[seq_len(k_val)]^2)
    var_resid <- max(var_total - var_top, .Machine$double.eps)
    -0.5 * (n * p_ctrl) * log(var_resid / (n * p_ctrl)) -
      0.5 * k_val * (n + p_ctrl) * log(n)
  }

  bics <- vapply(k_candidates,
                  function(k_val) -2 * log_lik(k_val) + k_val * (n + p_ctrl) * log(n),
                  numeric(1L))
  k_selected <- k_candidates[which.min(bics)]

  latent_train <- svd_res$u[, seq_len(k_selected), drop = FALSE] %*%
    diag(svd_res$d[seq_len(k_selected)], nrow = k_selected, ncol = k_selected)

  full_loadings <- matrix(NA_real_, nrow = ncol(X_centered), ncol = k_selected)
  for (j in seq_len(ncol(X_centered))) {
    full_loadings[j, ] <- coef(.lm.fit(latent_train, X_centered[, j]))
  }

  fit <- list(
    feat_mean = feat_mean,
    full_loadings = full_loadings,
    latent_train_inv = MASS::ginv(latent_train),
    k_selected = k_selected,
    bics = stats::setNames(bics, k_candidates),
    n_train = n,
    p_full = ncol(X),
    n_controls = length(ctrl_idx)
  )
  class(fit) <- "frozen_ruv_fit"
  fit
}


#' @title Apply frozen RUV correction to new samples
#'
#' @description
#' Removes the unwanted-variation component estimated by
#' \code{\link{fit_frozen_ruv}} from one or more new samples. The latent
#' factors are approximated from the frozen pseudo-inverse learned during
#' training.
#'
#' @param x Numeric vector (single sample) or matrix (samples \eqn{\times}
#'   features).
#' @param fit A \code{frozen_ruv_fit} object from \code{\link{fit_frozen_ruv}}.
#'
#' @return Same shape as \code{x}, with the unwanted-variation component
#'   removed and training feature means re-added. The attribute
#'   \code{fit_summary} summarises the fit parameters.
#'
#' @export
apply_frozen_ruv <- function(x, fit) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("apply_frozen_ruv requires the MASS package.")
  }
  stopifnot(inherits(fit, "frozen_ruv_fit"))
  if (is.vector(x)) {
    return(apply_frozen_ruv(
      matrix(x, nrow = 1L, dimnames = list(NULL, names(x))), fit)[1L, ])
  }
  X <- as.matrix(x)
  X_centered <- sweep(X, 2L, fit$feat_mean, "-")
  latent_new <- X_centered %*% fit$full_loadings %*%
    MASS::ginv(t(fit$full_loadings) %*% fit$full_loadings)
  unwanted <- latent_new %*% t(fit$full_loadings)
  X_clean <- X_centered - unwanted
  X_out <- sweep(X_clean, 2L, fit$feat_mean, "+")
  attr(X_out, "fit_summary") <- sprintf(
    "frozen_ruv (k=%d, n_train=%d, n_controls=%d)",
    fit$k_selected, fit$n_train, fit$n_controls)
  X_out
}


# ----------------------------------------------------------------------------
# fit_robust_pca_residual / apply_robust_pca_residual
# ----------------------------------------------------------------------------

#' @title Fit a MAD-scaled SVD residual filter for batch correction
#'
#' @description
#' Fits a frozen MAD-scaled PCA residual filter. This is NOT the Candès RPCA
#' (ALM sparse-low-rank decomposition) — it is classical SVD on a robustified
#' pre-processing step (per-feature median centering + MAD scaling), which is
#' more resistant to outlier samples than vanilla CLR-then-PCA but is not the
#' L+S sparse-low-rank framework.
#'
#' Training pipeline: per-feature median + MAD scaling, then classical SVD;
#' retain top-\code{n_components} right singular vectors as a frozen PC basis.
#'
#' @param x_train Numeric matrix (samples \eqn{\times} features). Must have
#'   column names. All values must be non-negative.
#' @param n_components Integer — number of principal components to retain
#'   (and remove). Default 5.
#' @param pseudocount Numeric \eqn{> 0}. Added before log transform. Default
#'   0.5.
#'
#' @return Object of class \code{robust_pca_residual_fit} with components:
#'   \code{feature_names}, \code{feat_median}, \code{feat_mad}, \code{v}
#'   (right singular vectors, features \eqn{\times} \code{n_components}),
#'   \code{n_components}, \code{pseudocount}, \code{method}.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x_train <- matrix(abs(rnorm(100 * 30, mean = 100, sd = 30)), nrow = 100,
#'                   dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' fit <- fit_robust_pca_residual(x_train, n_components = 3L)
#' }
#'
#' @references
#' Hubert M, Rousseeuw PJ, Vanden Branden K. (2005) ROBPCA.
#' \emph{Technometrics} 47(1): 64–79.
#'
#' @export
fit_robust_pca_residual <- function(x_train,
                                     n_components = 5L,
                                     pseudocount = 0.5) {
  if (!is.matrix(x_train)) x_train <- as.matrix(x_train)
  if (any(x_train < 0, na.rm = TRUE)) {
    stop("fit_robust_pca_residual: negative values not allowed (compositional input)")
  }
  if (pseudocount <= 0) stop("fit_robust_pca_residual: pseudocount must be > 0")
  if (n_components < 1L) stop("fit_robust_pca_residual: n_components must be >= 1")
  feat_names <- colnames(x_train)
  if (is.null(feat_names)) {
    stop("fit_robust_pca_residual: x_train must have feature names (colnames)")
  }
  if (n_components >= min(nrow(x_train) - 1L, ncol(x_train))) {
    stop("fit_robust_pca_residual: n_components (", n_components,
          ") must be < min(n_train - 1, n_features) = ",
          min(nrow(x_train) - 1L, ncol(x_train)))
  }

  log_x <- log(x_train + pseudocount)
  if (any(!is.finite(log_x))) {
    stop("fit_robust_pca_residual: non-finite log values after pseudocount addition")
  }
  feat_median <- apply(log_x, 2L, stats::median, na.rm = TRUE)
  feat_mad <- apply(log_x, 2L, stats::mad, na.rm = TRUE)
  feat_mad[is.na(feat_mad) | feat_mad < 1e-6] <- 1

  z_train <- sweep(log_x, 2L, feat_median, FUN = "-")
  z_train <- sweep(z_train, 2L, feat_mad, FUN = "/")

  svd_fit <- svd(z_train, nu = 0L, nv = n_components)
  v <- svd_fit$v

  fit <- list(
    feature_names = feat_names,
    feat_median = feat_median,
    feat_mad = feat_mad,
    v = v,
    n_components = n_components,
    pseudocount = pseudocount,
    method = "MAD-scaled SVD (not Candes RPCA)"
  )
  class(fit) <- "robust_pca_residual_fit"
  fit
}


#' @title Apply MAD-scaled SVD residual filter to new samples
#'
#' @description
#' Projects test samples into the frozen PC space estimated by
#' \code{\link{fit_robust_pca_residual}}, then returns the residual (test
#' sample minus reconstruction from top-\code{n_components} PCs) as a
#' low-rank-removed representation.
#'
#' All training features must be present in \code{x_test}.
#'
#' @param fit A \code{robust_pca_residual_fit} object from
#'   \code{\link{fit_robust_pca_residual}}.
#' @param x_test Numeric matrix (samples \eqn{\times} features). All values
#'   must be non-negative.
#' @param return_reconstructed Logical. If \code{TRUE}, returns a list with
#'   components \code{scores}, \code{reconstructed}, and \code{residual}.
#'   If \code{FALSE} (default), returns only the residual matrix.
#'
#' @return The residual matrix (samples \eqn{\times} features) on the
#'   MAD-scaled log scale, or a list when \code{return_reconstructed = TRUE}.
#'
#' @export
apply_robust_pca_residual <- function(fit, x_test, return_reconstructed = FALSE) {
  if (!inherits(fit, "robust_pca_residual_fit")) {
    stop("fit must be class robust_pca_residual_fit")
  }
  if (!is.matrix(x_test)) x_test <- as.matrix(x_test)
  if (any(x_test < 0, na.rm = TRUE)) {
    stop("apply_robust_pca_residual: negative values not allowed")
  }
  missing <- setdiff(fit$feature_names, colnames(x_test))
  if (length(missing) > 0L) {
    stop("apply_robust_pca_residual: ", length(missing),
          " training features missing from test set: ",
          paste(utils::head(missing, 5L), collapse = ", "), " ...")
  }
  x_test <- x_test[, fit$feature_names, drop = FALSE]
  log_x <- log(x_test + fit$pseudocount)
  z_test <- sweep(log_x, 2L, fit$feat_median, FUN = "-")
  z_test <- sweep(z_test, 2L, fit$feat_mad, FUN = "/")

  scores <- z_test %*% fit$v
  recon  <- scores %*% t(fit$v)
  resid  <- z_test - recon

  if (return_reconstructed) {
    list(scores = scores, reconstructed = recon, residual = resid)
  } else {
    resid
  }
}
