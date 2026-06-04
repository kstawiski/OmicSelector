#' @title Paper 3 additional within-sample methods (Module A, P2)
#'
#' @description
#' Two additional Module-A within-sample methods introduced in Paper 3
#' (Stawiski et al., in preparation), complementing the five core methods in
#' \code{singlesample-within-sample.R}.
#'
#' Methods provided:
#' \itemize{
#'   \item \code{\link{fit_logistic_normal_eb}} /
#'     \code{\link{apply_logistic_normal_eb}}: frozen-reference empirical-Bayes
#'     denoiser using logistic-normal posterior shrinkage (Aitchison & Shen
#'     1980; Efron 2010). Per-feature depth-weighted CLR shrinkage toward a
#'     training-cohort prior.
#'   \item \code{\link{fit_frozen_quantile}} /
#'     \code{\link{apply_frozen_quantile}}: monotone quantile calibrator for
#'     cross-platform mapping. Frozen empirical-CDF mapping from test-sample
#'     feature values to training-distribution quantile values (Bolstad et al.
#'     2003; Hicks et al. 2018).
#' }
#'
#' @references
#' Aitchison J, Shen SM. (1980) Logistic-normal distributions: some properties
#' and uses. \emph{Biometrika} 67(2): 261–272.
#'
#' Efron B. (2010) Large-Scale Inference. Cambridge University Press.
#'
#' Bolstad BM, Irizarry RA, Astrand M, Speed TP. (2003) A comparison of
#' normalization methods for high density oligonucleotide array data based on
#' variance and bias. \emph{Bioinformatics} 19(2): 185–193.
#'
#' Hicks SC, Okrah K, Paulson JN, Quackenbush J, Irizarry RA, Bravo HC. (2018)
#' Smooth quantile normalization. \emph{Biostatistics} 19(2): 185–198.
#'
#' Stawiski K. (in preparation) Paper 3 of the OmicSelector programme.
#'
#' @name singlesample-additional-within-sample
NULL


# ============================================================================
# Logistic-normal Empirical-Bayes within-sample denoiser
# ============================================================================

#' @title Fit frozen logistic-normal empirical-Bayes prior from training data
#'
#' @description
#' Estimates per-feature CLR mean and variance on a training cohort. At apply
#' time (\code{\link{apply_logistic_normal_eb}}), each test-sample CLR
#' observation is shrunk toward the per-feature prior using a
#' precision-weighted Bayesian update:
#' \deqn{\hat{z}_{i,f} = w_f z_{i,f} + (1 - w_f) \mu_f,}
#' where \eqn{w_f = \sigma^2_f / (\sigma^2_f + \sigma^2_{\mathrm{obs},i})}.
#' The observation variance \eqn{\sigma^2_{\mathrm{obs},i}} scales
#' inversely with sample total counts so shallow-coverage samples shrink
#' harder. Per-feature variance is truly load-bearing (not a scalar heuristic).
#'
#' @param x_train Numeric matrix (samples \eqn{\times} features). Must have
#'   column names. All values must be non-negative.
#' @param pseudocount Numeric \eqn{> 0}. Added before log transform. Default
#'   0.5.
#'
#' @return Object of class \code{logistic_normal_eb_fit} with components:
#'   \code{feature_names}, \code{pseudocount}, \code{prior_mean},
#'   \code{prior_var}, \code{median_total_train}.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x_train <- matrix(abs(rnorm(50 * 20, 100, 30)), nrow = 50,
#'                   dimnames = list(NULL, paste0("miR-", seq_len(20))))
#' fit <- fit_logistic_normal_eb(x_train)
#' length(fit$prior_var)
#' }
#'
#' @references
#' Aitchison J, Shen SM. (1980) \emph{Biometrika} 67: 261–272.
#' Efron B. (2010) Large-Scale Inference. Cambridge UP.
#'
#' @export
fit_logistic_normal_eb <- function(x_train, pseudocount = 0.5) {
  if (!is.matrix(x_train)) x_train <- as.matrix(x_train)
  if (any(x_train < 0, na.rm = TRUE)) {
    stop("fit_logistic_normal_eb: negative values not allowed (compositional input)")
  }
  if (pseudocount <= 0) stop("fit_logistic_normal_eb: pseudocount must be > 0")
  feat_names <- colnames(x_train)
  if (is.null(feat_names)) {
    stop("fit_logistic_normal_eb: x_train must have feature names (colnames)")
  }

  log_x <- log(x_train + pseudocount)
  centering <- rowMeans(log_x, na.rm = TRUE)
  z_train <- sweep(log_x, 1L, centering, FUN = "-")

  prior_mean <- colMeans(z_train, na.rm = TRUE)
  prior_var <- apply(z_train, 2L, stats::var, na.rm = TRUE)
  prior_var[prior_var < 1e-6 | is.na(prior_var)] <- 1e-6

  total_counts_train <- rowSums(x_train, na.rm = TRUE)
  median_total_train <- stats::median(total_counts_train, na.rm = TRUE)

  fit <- list(
    feature_names = feat_names,
    pseudocount = pseudocount,
    prior_mean = prior_mean,
    prior_var = prior_var,
    median_total_train = median_total_train
  )
  class(fit) <- "logistic_normal_eb_fit"
  fit
}


#' @title Apply frozen logistic-normal EB shrinkage to new samples
#'
#' @description
#' Shrinks CLR-transformed test samples toward the per-feature training prior
#' using the model fitted by \code{\link{fit_logistic_normal_eb}}.
#'
#' Missing test features are handled by \code{missing_features}:
#' \code{"error"} (default, fail-closed) or \code{"fill_with_prior"} (pads
#' missing features with \code{NA} in the test matrix, then the posterior
#' collapses fully to the prior mean for those features).
#'
#' @param fit A \code{logistic_normal_eb_fit} object from
#'   \code{\link{fit_logistic_normal_eb}}.
#' @param x_test Numeric matrix (samples \eqn{\times} features). All values
#'   must be non-negative.
#' @param missing_features \code{"error"} (default) or
#'   \code{"fill_with_prior"}.
#'
#' @return Numeric matrix (same row/col structure as ordered training features)
#'   of posterior-shrunk CLR coordinates.
#'
#' @export
apply_logistic_normal_eb <- function(fit, x_test,
                                      missing_features = c("error", "fill_with_prior")) {
  missing_features <- match.arg(missing_features)
  if (!inherits(fit, "logistic_normal_eb_fit")) {
    stop("fit must be class logistic_normal_eb_fit")
  }
  if (!is.matrix(x_test)) x_test <- as.matrix(x_test)
  if (any(x_test < 0, na.rm = TRUE)) {
    stop("apply_logistic_normal_eb: negative values not allowed (compositional input)")
  }
  if (is.null(colnames(x_test))) {
    stop("apply_logistic_normal_eb: x_test must have feature names")
  }

  missing <- setdiff(fit$feature_names, colnames(x_test))
  if (length(missing) > 0L) {
    if (missing_features == "error") {
      stop("apply_logistic_normal_eb: ", length(missing),
            " training features missing from test set: ",
            paste(utils::head(missing, 5L), collapse = ", "))
    }
    pad <- matrix(NA_real_, nrow = nrow(x_test), ncol = length(missing),
                   dimnames = list(rownames(x_test), missing))
    x_test <- cbind(x_test, pad)
  }
  x_test <- x_test[, fit$feature_names, drop = FALSE]

  log_x <- log(x_test + fit$pseudocount)
  centering <- rowMeans(log_x, na.rm = TRUE)
  z_test <- sweep(log_x, 1L, centering, FUN = "-")

  total_counts_test <- rowSums(x_test, na.rm = TRUE)
  total_counts_test[total_counts_test <= 0] <- 1
  median_prior_var <- stats::median(fit$prior_var)
  obs_var_per_sample <- (fit$median_total_train / total_counts_test) * median_prior_var

  z_post <- matrix(NA_real_, nrow = nrow(z_test), ncol = ncol(z_test),
                    dimnames = dimnames(z_test))
  for (i in seq_len(nrow(z_test))) {
    sigma2_obs_i <- obs_var_per_sample[i]
    w_obs <- fit$prior_var / (fit$prior_var + sigma2_obs_i)
    z_obs <- z_test[i, ]
    z_obs[is.na(z_obs)] <- fit$prior_mean[is.na(z_obs)]
    z_post[i, ] <- w_obs * z_obs + (1 - w_obs) * fit$prior_mean
  }
  z_post
}


# ============================================================================
# Frozen monotone quantile calibrator
# ============================================================================

#' @title Fit a frozen monotone quantile calibrator from training data
#'
#' @description
#' Estimates per-feature empirical quantile grids on a training cohort and
#' freezes them for deployment. At apply time
#' (\code{\link{apply_frozen_quantile}}), test-sample feature values are
#' mapped to their empirical quantile under the training distribution and then
#' transformed back to the corresponding training-scale value. Tail values
#' beyond the training range are clamped (intentional winsorisation, documented
#' limitation).
#'
#' @param x_train Numeric matrix (samples \eqn{\times} features). Must have
#'   column names. All values must be non-negative.
#' @param n_quantiles Integer — number of quantile grid levels (resolution).
#'   Default 1000.
#' @param na_action \code{"median"} (default, replaces NA test values with
#'   per-feature training median) or \code{"leave_na"}.
#'
#' @return Object of class \code{frozen_quantile_fit} with components:
#'   \code{quantile_ref} (matrix of quantile values), \code{quantile_levels},
#'   \code{feature_names}, \code{median_per_feat}, \code{n_quantiles},
#'   \code{na_action}.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x_train <- matrix(abs(rnorm(50 * 20, 100, 30)), nrow = 50,
#'                   dimnames = list(NULL, paste0("miR-", seq_len(20))))
#' fit <- fit_frozen_quantile(x_train, n_quantiles = 100L)
#' }
#'
#' @references
#' Bolstad BM, et al. (2003) \emph{Bioinformatics} 19: 185–193.
#'
#' @export
fit_frozen_quantile <- function(x_train,
                                 n_quantiles = 1000L,
                                 na_action = c("median", "leave_na")) {
  na_action <- match.arg(na_action)
  if (!is.matrix(x_train)) x_train <- as.matrix(x_train)
  if (any(x_train < 0, na.rm = TRUE)) {
    stop("fit_frozen_quantile: negative values not allowed (compositional input)")
  }
  feat_names <- colnames(x_train)
  if (is.null(feat_names)) {
    stop("fit_frozen_quantile: x_train must have feature names (colnames)")
  }

  qs <- seq(0, 1, length.out = n_quantiles + 1L)
  ref <- apply(x_train, 2L, function(col) {
    if (all(is.na(col))) return(rep(NA_real_, length(qs)))
    stats::quantile(col, probs = qs, na.rm = TRUE, names = FALSE)
  })

  median_per_feat <- apply(x_train, 2L, stats::median, na.rm = TRUE)

  fit <- list(
    quantile_ref = ref,
    quantile_levels = qs,
    feature_names = feat_names,
    median_per_feat = median_per_feat,
    n_quantiles = n_quantiles,
    na_action = na_action
  )
  class(fit) <- "frozen_quantile_fit"
  fit
}


#' @title Apply frozen quantile calibration to new samples
#'
#' @description
#' Maps test-sample feature values to the training distribution using the
#' frozen empirical-CDF grid from \code{\link{fit_frozen_quantile}}. Each
#' test value is located in the training CDF, and the corresponding training
#' reference value is returned. This aligns the test distribution to the
#' training distribution in a monotone, feature-wise manner.
#'
#' @param fit A \code{frozen_quantile_fit} object from
#'   \code{\link{fit_frozen_quantile}}.
#' @param x_test Numeric matrix (samples \eqn{\times} features). All values
#'   must be non-negative.
#' @param missing_features \code{"error"} (default, fail-closed when training
#'   features are absent) or \code{"skip"} (process only shared features and
#'   return a reduced-column matrix).
#'
#' @return Numeric matrix with columns for the shared features, values on the
#'   training quantile scale.
#'
#' @export
apply_frozen_quantile <- function(fit, x_test,
                                   missing_features = c("error", "skip")) {
  missing_features <- match.arg(missing_features)
  if (!is.matrix(x_test)) x_test <- as.matrix(x_test)
  if (any(x_test < 0, na.rm = TRUE)) {
    stop("apply_frozen_quantile: negative values not allowed")
  }
  if (!inherits(fit, "frozen_quantile_fit")) {
    stop("fit must be class frozen_quantile_fit")
  }

  shared <- intersect(fit$feature_names, colnames(x_test))
  missing <- setdiff(fit$feature_names, colnames(x_test))
  if (length(missing) > 0L && missing_features == "error") {
    stop("apply_frozen_quantile: ", length(missing),
          " training features missing from test set: ",
          paste(utils::head(missing, 5L), collapse = ", "), " ...")
  }
  if (length(shared) == 0L) {
    stop("apply_frozen_quantile: no shared features between fit and x_test")
  }

  out <- matrix(NA_real_, nrow = nrow(x_test), ncol = length(shared),
                dimnames = list(rownames(x_test), shared))

  for (f in shared) {
    ref_col <- fit$quantile_ref[, f]
    if (all(is.na(ref_col))) next
    test_col <- x_test[, f]
    nonna <- !is.na(test_col)
    if (sum(nonna) == 0L) next

    train_cdf <- function(v) {
      idx <- findInterval(v, ref_col, all.inside = TRUE)
      idx / length(ref_col)
    }
    test_q <- train_cdf(test_col[nonna])
    out_idx <- pmax(1L, pmin(length(ref_col),
                              round(test_q * (length(ref_col) - 1L)) + 1L))
    out[nonna, f] <- ref_col[out_idx]

    if (fit$na_action == "median" && sum(!nonna) > 0L) {
      out[!nonna, f] <- fit$median_per_feat[f]
    }
  }
  out
}
