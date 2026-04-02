#' @title Safe Evaluation Utilities for Biomarker Validation
#'
#' @description
#' Functions that enforce correct ROC computation and non-inferiority testing.
#' These utilities prevent the most common evaluation bugs:
#' hardcoded ROC direction (inverting AUC), wrong non-inferiority references,
#' and unstable percentile estimates from insufficient random draws.
#'
#' @name safe-evaluation
NULL


#' @title Safe ROC Computation (Direction-Guarded)
#'
#' @description
#' Wrapper around \code{pROC::roc()} that enforces \code{direction = "auto"}
#' by default and warns loudly when direction is hardcoded.
#'
#' **Common bug this prevents**: Hardcoding \code{direction = "<"} when the
#' predictor's score direction depends on class weights, calibration, or
#' model type. With extreme class weighting (e.g., 17x), the positive class
#' can have LOWER predicted probabilities, inverting the AUC to 1 - true_AUC.
#'
#' @param response Binary outcome vector (0/1 or factor)
#' @param predictor Numeric prediction scores
#' @param direction ROC direction. Default \code{"auto"} (strongly recommended).
#'   If set to \code{"<"} or \code{">"}, a warning is issued.
#' @param levels Two-element vector specifying the levels of the response.
#'   Default \code{c(0, 1)}.
#' @param quiet Logical; suppress pROC messages. Default TRUE.
#' @param ... Additional arguments passed to \code{pROC::roc()}
#'
#' @return A pROC roc object
#'
#' @examples
#' \dontrun{
#' roc_obj <- safe_roc(y_true, y_pred)
#' auc_val <- as.numeric(pROC::auc(roc_obj))
#' }
#'
#' @export
safe_roc <- function(response, predictor, direction = "auto",
                     levels = c(0, 1), quiet = TRUE, ...) {

  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required. Install with: install.packages('pROC')", call. = FALSE)
  }

  if (direction != "auto") {
    warning(
      sprintf(
        "Hardcoded direction='%s' in ROC computation. ",
        direction
      ),
      "This can INVERT the AUC when class weights or calibration change the ",
      "score direction. Use direction='auto' unless you have verified the ",
      "predictor's score polarity. See: OmicSelector methodology audit, ",
      "pROC direction bug (miRPOC 2026-04-02).",
      call. = FALSE
    )
  }

  pROC::roc(
    response  = response,
    predictor = predictor,
    direction = direction,
    levels    = levels,
    quiet     = quiet,
    ...
  )
}


#' @title Safe AUC with Confidence Interval
#'
#' @description
#' Computes AUC with DeLong confidence interval using \code{safe_roc()}.
#' Returns a named list with AUC, CI, and the ROC object.
#'
#' @param response Binary outcome vector (0/1)
#' @param predictor Numeric prediction scores
#' @param ci_method CI method: "delong" (default) or "bootstrap"
#' @param ... Additional arguments passed to \code{safe_roc()}
#'
#' @return A list with:
#'   \item{auc}{Numeric AUC value}
#'   \item{ci_lower}{Lower bound of 95\% CI}
#'   \item{ci_upper}{Upper bound of 95\% CI}
#'   \item{roc}{The pROC roc object}
#'   \item{direction}{Direction used by pROC}
#'
#' @export
safe_auc <- function(response, predictor, ci_method = "delong", ...) {

  roc_obj <- safe_roc(response, predictor, ...)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_obj  <- pROC::ci.auc(roc_obj, method = ci_method)

  list(
    auc      = auc_val,
    ci_lower = as.numeric(ci_obj[1]),
    ci_upper = as.numeric(ci_obj[3]),
    roc      = roc_obj,
    direction = roc_obj$direction
  )
}


#' @title Paired DeLong Non-Inferiority Test
#'
#' @description
#' Tests whether a candidate model's AUC is non-inferior to a reference AUC,
#' using a paired DeLong test on the same predictions.
#'
#' **Common bugs this prevents**:
#' \enumerate{
#'   \item Using the wrong reference AUC (e.g., superseded Phase 2.1 instead of
#'         current Phase 2.1b)
#'   \item Using an unpaired CI check instead of a formal paired DeLong test
#'   \item Comparing AUCs from different sample denominators
#' }
#'
#' @param response Binary outcome vector (0/1), same for both models
#' @param predictor_candidate Candidate model predictions
#' @param predictor_reference Reference model predictions. Must be on the
#'   SAME samples as predictor_candidate.
#' @param delta Non-inferiority margin (positive value, e.g., 0.03).
#'   Candidate is non-inferior if its AUC is within delta of the reference.
#' @param alpha Significance level (default 0.05 for one-sided test)
#'
#' @return A list with:
#'   \item{non_inferior}{Logical: TRUE if candidate is non-inferior}
#'   \item{auc_candidate}{Candidate AUC}
#'   \item{auc_reference}{Reference AUC}
#'   \item{delta_auc}{Candidate - Reference AUC}
#'   \item{delta}{Non-inferiority margin used}
#'   \item{se_diff}{Standard error of AUC difference (paired DeLong)}
#'   \item{z_stat}{Z-statistic for non-inferiority}
#'   \item{p_value}{One-sided p-value}
#'   \item{ci_lower_diff}{Lower bound of 95\% CI for AUC difference}
#'
#' @examples
#' \dontrun{
#' result <- test_noninferiority(
#'   response = y_test,
#'   predictor_candidate = pred_reduced_panel,
#'   predictor_reference = pred_full_panel,
#'   delta = 0.03
#' )
#' cat("Non-inferior:", result$non_inferior, "\n")
#' }
#'
#' @export
test_noninferiority <- function(response, predictor_candidate, predictor_reference,
                                delta = 0.03, alpha = 0.05) {

  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required.", call. = FALSE)
  }

  if (length(response) != length(predictor_candidate) ||
      length(response) != length(predictor_reference)) {
    stop("response, predictor_candidate, and predictor_reference must have equal length",
         call. = FALSE)
  }

  if (delta <= 0) {
    stop("delta must be positive (e.g., 0.03)", call. = FALSE)
  }

  # Compute ROC objects (direction = auto for both)
  roc_cand <- safe_roc(response, predictor_candidate)
  roc_ref  <- safe_roc(response, predictor_reference)

  auc_cand <- as.numeric(pROC::auc(roc_cand))
  auc_ref  <- as.numeric(pROC::auc(roc_ref))
  delta_auc <- auc_cand - auc_ref

  # Paired DeLong variance for AUC difference
  # Use pROC::roc.test for the paired comparison
  delong_test <- pROC::roc.test(roc_ref, roc_cand, method = "delong", paired = TRUE)

  # Extract standard error from the DeLong test
  # The test statistic is (AUC1 - AUC2) / SE
  # We need SE for the non-inferiority formulation
  se_diff <- abs(delta_auc / delong_test$statistic)
  if (!is.finite(se_diff) || se_diff == 0) {
    # Fallback: estimate from the DeLong CI
    se_diff <- abs(delta_auc) / abs(qnorm(delong_test$p.value / 2))
    if (!is.finite(se_diff)) se_diff <- NA_real_
  }

  # Non-inferiority test: H0: AUC_cand - AUC_ref < -delta
  #                       H1: AUC_cand - AUC_ref >= -delta
  # Z = (delta_auc + delta) / SE
  if (is.finite(se_diff) && se_diff > 0) {
    z_stat <- (delta_auc + delta) / se_diff
    p_value <- pnorm(z_stat, lower.tail = FALSE)  # One-sided
    # Actually for non-inferiority: we want to reject H0: diff < -delta
    # So p = P(Z > z) under H0
    p_value <- 1 - pnorm(z_stat)
    ci_lower_diff <- delta_auc - qnorm(1 - alpha) * se_diff
  } else {
    z_stat <- NA_real_
    p_value <- NA_real_
    ci_lower_diff <- NA_real_
  }

  non_inferior <- !is.na(p_value) && p_value < alpha

  list(
    non_inferior   = non_inferior,
    auc_candidate  = auc_cand,
    auc_reference  = auc_ref,
    delta_auc      = delta_auc,
    delta          = delta,
    se_diff        = as.numeric(se_diff),
    z_stat         = as.numeric(z_stat),
    p_value        = as.numeric(p_value),
    ci_lower_diff  = as.numeric(ci_lower_diff)
  )
}


#' @title Validate Random-Panel Null Benchmark
#'
#' @description
#' Checks whether the number of random draws is sufficient for a stable
#' percentile estimate. Warns when fewer than 500 draws are used.
#'
#' @param n_draws Number of random panel draws performed
#' @param percentile The percentile being estimated (e.g., 0.95)
#'
#' @return Invisible TRUE if adequate, FALSE with warning otherwise
#'
#' @export
check_null_benchmark_draws <- function(n_draws, percentile = 0.95) {
  # Order statistic SE approximation
  min_recommended <- ceiling(20 / (1 - percentile))  # ~400 for 95th pctile

  if (n_draws < min_recommended) {
    warning(
      sprintf(
        "Only %d draws for %.0fth percentile estimate (recommend >= %d). ",
        n_draws, percentile * 100, min_recommended
      ),
      "The Monte Carlo error on this order statistic is non-trivial. ",
      "Consider rerunning with more draws or reporting a bootstrap CI ",
      "on the percentile.",
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  invisible(TRUE)
}
