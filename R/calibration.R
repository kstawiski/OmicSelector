#' @title Calibration Metrics for OmicSelector 2.0
#'
#' @description
#' Functions for assessing and improving probability calibration of
#' classification models. Well-calibrated probabilities are essential
#' for clinical decision-making.
#'
#' @details
#' ## Why Calibration Matters
#'
#' A model with high AUC can still produce poorly calibrated probabilities.
#' For clinical use, we need:
#' - When model says 80% probability, ~80% of such cases should be positive
#' - Brier score decomposes into discrimination + calibration + uncertainty
#' - ECE (Expected Calibration Error) measures mean deviation from perfect calibration
#'
#' ## Calibration Repair
#'
#' Post-hoc calibration can improve probability estimates:
#' - **Platt scaling**: Fits logistic regression on predictions
#' - **Isotonic regression**: Non-parametric monotonic recalibration
#' - **Temperature scaling**: Single parameter softmax recalibration
#'
#' @references
#' Niculescu-Mizil, A., & Caruana, R. (2005). Predicting good probabilities
#' with supervised learning. ICML.
#'
#' @name calibration
NULL


#' @title Compute Expected Calibration Error (ECE)
#'
#' @description
#' Calculates the Expected Calibration Error, which measures how well
#' predicted probabilities match observed frequencies. Lower is better.
#'
#' @param probs Numeric vector of predicted probabilities (0-1)
#' @param labels Binary labels (0/1 or logical)
#' @param n_bins Number of bins for grouping probabilities (default: 10)
#' @param weighting How to weight bins: "samples" (default) or "uniform"
#'
#' @return A list containing:
#'   - ece: Expected Calibration Error (0-1, lower is better)
#'   - mce: Maximum Calibration Error
#'   - bin_data: Data frame with per-bin statistics
#'
#' @details
#' ECE is computed as:
#' \deqn{ECE = \sum_{b=1}^{B} \frac{n_b}{N} |accuracy_b - confidence_b|}
#'
#' where B is the number of bins, n_b is samples in bin b, N is total samples,
#' accuracy_b is fraction of positives in bin, and confidence_b is mean
#' predicted probability in bin.
#'
#' @examples
#' \dontrun{
#' # Simulated predictions
#' probs <- runif(100)
#' labels <- rbinom(100, 1, probs)  # Well-calibrated by construction
#'
#' result <- compute_ece(probs, labels)
#' print(result$ece)  # Should be low
#' }
#'
#' @export
compute_ece <- function(probs, labels, n_bins = 10, weighting = c("samples", "uniform")) {
  weighting <- match.arg(weighting)

  # Convert factor labels to numeric (0/1) before validation
  # Factors are common in mlr3 classification workflows
  if (is.factor(labels)) {
    labels <- as.integer(labels) - 1L  # Factor levels to 0/1
  }

  # Validate inputs
  checkmate::assert_numeric(probs, lower = 0, upper = 1, any.missing = FALSE)
  checkmate::assert_integerish(labels, lower = 0, upper = 1, any.missing = FALSE)
  checkmate::assert_true(length(probs) == length(labels))

  labels <- as.numeric(labels)
  n <- length(probs)

  # Create bins
  bin_edges <- seq(0, 1, length.out = n_bins + 1)
  bin_indices <- findInterval(probs, bin_edges, all.inside = TRUE)

  # Compute per-bin statistics
  bin_data <- data.frame(
    bin = 1:n_bins,
    lower = bin_edges[1:n_bins],
    upper = bin_edges[2:(n_bins + 1)],
    n_samples = 0,
    mean_prob = NA_real_,
    mean_label = NA_real_,
    calibration_error = NA_real_
  )

  for (b in 1:n_bins) {
    idx <- bin_indices == b
    n_b <- sum(idx)

    if (n_b > 0) {
      bin_data$n_samples[b] <- n_b
      bin_data$mean_prob[b] <- mean(probs[idx])
      bin_data$mean_label[b] <- mean(labels[idx])
      bin_data$calibration_error[b] <- abs(bin_data$mean_label[b] - bin_data$mean_prob[b])
    }
  }

  # Filter to non-empty bins
  non_empty <- bin_data$n_samples > 0

  if (sum(non_empty) == 0) {
    return(list(ece = NA, mce = NA, bin_data = bin_data))
  }

  # Compute ECE
  if (weighting == "samples") {
    weights <- bin_data$n_samples[non_empty] / n
  } else {
    weights <- rep(1 / sum(non_empty), sum(non_empty))
  }

  ece <- sum(weights * bin_data$calibration_error[non_empty])
  mce <- max(bin_data$calibration_error[non_empty])

  list(
    ece = ece,
    mce = mce,
    n_bins = n_bins,
    n_samples = n,
    bin_data = bin_data
  )
}


#' @title Decompose Brier Score
#'
#' @description
#' Decomposes Brier score into reliability (calibration), resolution
#' (discrimination), and uncertainty components.
#'
#' @param probs Numeric vector of predicted probabilities
#' @param labels Binary labels (0/1)
#'
#' @return A list containing:
#'   - brier: Total Brier score
#'   - reliability: Calibration component (lower is better)
#'   - resolution: Discrimination component (higher is better)
#'   - uncertainty: Inherent uncertainty (fixed for given class balance)
#'
#' @details
#' Brier score decomposes as:
#' \deqn{Brier = Reliability - Resolution + Uncertainty}
#'
#' - **Reliability** (calibration): measures deviation from perfect calibration
#' - **Resolution**: measures how much predictions differ across groups
#' - **Uncertainty**: base rate of positive class, p(1-p)
#'
#' A well-calibrated model has low reliability and high resolution.
#'
#' @examples
#' \dontrun{
#' probs <- c(0.1, 0.4, 0.6, 0.9)
#' labels <- c(0, 0, 1, 1)
#' decompose_brier(probs, labels)
#' }
#'
#' @export
decompose_brier <- function(probs, labels) {
  # Convert factor labels to numeric (0/1) before validation
  if (is.factor(labels)) {
    labels <- as.integer(labels) - 1L
  }

  checkmate::assert_numeric(probs, lower = 0, upper = 1, any.missing = FALSE)
  checkmate::assert_integerish(labels, lower = 0, upper = 1, any.missing = FALSE)

  labels <- as.numeric(labels)
  n <- length(probs)

  # Total Brier score
  brier <- mean((probs - labels)^2)

  # Base rate (uncertainty)
  p_bar <- mean(labels)
  uncertainty <- p_bar * (1 - p_bar)

  # Use binned estimates for decomposition
  n_bins <- min(10, floor(sqrt(n)))
  bin_edges <- seq(0, 1, length.out = n_bins + 1)
  bin_indices <- findInterval(probs, bin_edges, all.inside = TRUE)

  reliability <- 0
  resolution <- 0

  for (b in 1:n_bins) {
    idx <- bin_indices == b
    n_b <- sum(idx)

    if (n_b > 0) {
      p_b <- mean(probs[idx])  # Mean predicted probability
      o_b <- mean(labels[idx])  # Observed frequency

      # Reliability: how far observed is from predicted
      reliability <- reliability + n_b * (o_b - p_b)^2

      # Resolution: how far observed is from overall base rate
      resolution <- resolution + n_b * (o_b - p_bar)^2
    }
  }

  reliability <- reliability / n
  resolution <- resolution / n

  list(
    brier = brier,
    reliability = reliability,
    resolution = resolution,
    uncertainty = uncertainty,
    # Sanity check: should approximately equal brier
    reconstructed = reliability - resolution + uncertainty
  )
}


#' @title Platt Scaling (Logistic Calibration)
#'
#' @description
#' Fits a logistic regression model to recalibrate probabilities.
#' Parameters should be estimated on a held-out calibration set.
#'
#' @param probs Calibration set predicted probabilities
#' @param labels Calibration set labels
#'
#' @return A function that recalibrates probabilities
#'
#' @details
#' Platt scaling fits: P(y=1|f) = 1 / (1 + exp(A*f + B))
#' where f is the original score/probability and A, B are fitted parameters.
#'
#' @examples
#' \dontrun{
#' # Fit on calibration set
#' cal_probs <- runif(100)
#' cal_labels <- rbinom(100, 1, cal_probs)
#' calibrator <- fit_platt_scaling(cal_probs, cal_labels)
#'
#' # Apply to new predictions
#' new_probs <- c(0.2, 0.5, 0.8)
#' calibrated <- calibrator(new_probs)
#' }
#'
#' @export
fit_platt_scaling <- function(probs, labels) {
  # Convert factor labels to numeric (0/1) before validation
  if (is.factor(labels)) {
    labels <- as.integer(labels) - 1L
  }

  checkmate::assert_numeric(probs, lower = 0, upper = 1, any.missing = FALSE)
  checkmate::assert_integerish(labels, lower = 0, upper = 1, any.missing = FALSE)

  labels <- as.numeric(labels)

  # Transform probabilities to log-odds for fitting
  # Clip to avoid -Inf/Inf
  probs_clipped <- pmax(pmin(probs, 1 - 1e-10), 1e-10)

  # Fit logistic regression
  df <- data.frame(y = labels, score = probs_clipped)
  fit <- suppressWarnings(glm(y ~ score, data = df, family = binomial()))

  # Return calibration function
  function(new_probs) {
    new_probs_clipped <- pmax(pmin(new_probs, 1 - 1e-10), 1e-10)
    new_df <- data.frame(score = new_probs_clipped)
    as.numeric(predict(fit, newdata = new_df, type = "response"))
  }
}


#' @title Isotonic Regression Calibration
#'
#' @description
#' Fits isotonic (monotonic) regression for probability calibration.
#' Non-parametric method that produces a monotonically increasing
#' calibration mapping.
#'
#' @param probs Calibration set predicted probabilities
#' @param labels Calibration set labels
#'
#' @return A function that recalibrates probabilities
#'
#' @details
#' Isotonic regression finds the best monotonically increasing function
#' that maps scores to calibrated probabilities. It's more flexible than
#' Platt scaling but requires more calibration data.
#'
#' IMPORTANT: Must be fit on a held-out calibration set, not training data,
#' to avoid leakage.
#'
#' @examples
#' \dontrun{
#' # Fit on calibration set
#' calibrator <- fit_isotonic_calibration(cal_probs, cal_labels)
#'
#' # Apply to new predictions
#' calibrated <- calibrator(new_probs)
#' }
#'
#' @export
fit_isotonic_calibration <- function(probs, labels) {
  # Convert factor labels to numeric (0/1) before validation
  if (is.factor(labels)) {
    labels <- as.integer(labels) - 1L
  }

  checkmate::assert_numeric(probs, lower = 0, upper = 1, any.missing = FALSE)
  checkmate::assert_integerish(labels, lower = 0, upper = 1, any.missing = FALSE)

  labels <- as.numeric(labels)

  # Sort by predicted probability
  ord <- order(probs)
  probs_sorted <- probs[ord]
  labels_sorted <- labels[ord]

  # Fit isotonic regression
  iso_fit <- isoreg(probs_sorted, labels_sorted)

  # Store for prediction
  # isoreg returns stepfun via stepfun(iso_fit)

  function(new_probs) {
    # Use stepfun for prediction
    sf <- as.stepfun(iso_fit)
    calibrated <- sf(new_probs)
    # Clip to [0, 1]
    pmax(pmin(calibrated, 1), 0)
  }
}


#' @title Temperature Scaling
#'
#' @description
#' Fits temperature scaling for probability calibration.
#' Uses a single temperature parameter T to soften/sharpen predictions.
#'
#' @param probs Calibration set predicted probabilities
#' @param labels Calibration set labels
#'
#' @return A function that recalibrates probabilities
#'
#' @details
#' Temperature scaling modifies the logits by dividing by T:
#' \deqn{p_{calibrated} = \sigma(logit(p) / T)}
#'
#' - T > 1: softens predictions (moves toward 0.5)
#' - T < 1: sharpens predictions (moves toward 0 or 1)
#' - T = 1: no change
#'
#' Commonly used for neural network calibration.
#'
#' @export
fit_temperature_scaling <- function(probs, labels) {
  # Convert factor labels to numeric (0/1) before validation
  if (is.factor(labels)) {
    labels <- as.integer(labels) - 1L
  }

  checkmate::assert_numeric(probs, lower = 0, upper = 1, any.missing = FALSE)
  checkmate::assert_integerish(labels, lower = 0, upper = 1, any.missing = FALSE)

  labels <- as.numeric(labels)

  # Clip probabilities
  probs_clipped <- pmax(pmin(probs, 1 - 1e-10), 1e-10)

  # Convert to logits
  logits <- log(probs_clipped / (1 - probs_clipped))

  # Objective: minimize negative log-likelihood
  nll <- function(T) {
    scaled_probs <- 1 / (1 + exp(-logits / T))
    -mean(labels * log(scaled_probs + 1e-10) + (1 - labels) * log(1 - scaled_probs + 1e-10))
  }

  # Optimize T
  opt <- optimize(nll, interval = c(0.1, 10))
  T_opt <- opt$minimum

  function(new_probs) {
    new_probs_clipped <- pmax(pmin(new_probs, 1 - 1e-10), 1e-10)
    new_logits <- log(new_probs_clipped / (1 - new_probs_clipped))
    1 / (1 + exp(-new_logits / T_opt))
  }
}


#' @title Create Reliability Diagram Data
#'
#' @description
#' Generates data for plotting a reliability (calibration) diagram.
#'
#' @param probs Predicted probabilities
#' @param labels True labels
#' @param n_bins Number of bins
#'
#' @return A data frame suitable for ggplot2
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' diagram_data <- reliability_diagram_data(probs, labels)
#'
#' ggplot(diagram_data, aes(x = mean_predicted, y = fraction_positive)) +
#'   geom_point(aes(size = n_samples)) +
#'   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#'   xlim(0, 1) + ylim(0, 1) +
#'   labs(x = "Mean Predicted Probability", y = "Fraction of Positives")
#' }
#'
#' @export
reliability_diagram_data <- function(probs, labels, n_bins = 10) {
  ece_result <- compute_ece(probs, labels, n_bins = n_bins)

  data.frame(
    bin = ece_result$bin_data$bin,
    bin_lower = ece_result$bin_data$lower,
    bin_upper = ece_result$bin_data$upper,
    mean_predicted = ece_result$bin_data$mean_prob,
    fraction_positive = ece_result$bin_data$mean_label,
    n_samples = ece_result$bin_data$n_samples,
    calibration_error = ece_result$bin_data$calibration_error
  )
}


#' @title Calibration Summary for Model Results
#'
#' @description
#' Computes comprehensive calibration metrics from benchmark results.
#' Intended for use with BenchmarkService outputs.
#'
#' This function is metrics-only and does NOT fit calibration repair models.
#' To fit calibrators, use the separate functions \code{fit_platt_scaling()},
#' \code{fit_isotonic_calibration()}, or \code{fit_temperature_scaling()}
#' on a held-out calibration set (NOT the evaluation/test set).
#'
#' @param probs Vector or list of predicted probabilities
#' @param labels Vector or list of true labels
#'
#' @return A CalibrationResult object with calibration metrics
#'
#' @details
#' To avoid data leakage, calibration repair models should be fit on a
#' separate calibration set, not on the evaluation data used for metrics.
#' This function only computes metrics; fitting calibrators is a separate step.
#'
#' @seealso \code{\link{fit_platt_scaling}}, \code{\link{fit_isotonic_calibration}},
#'   \code{\link{fit_temperature_scaling}}
#'
#' @export
calibration_summary <- function(probs, labels) {
  # Handle list inputs (from CV folds)
  if (is.list(probs) && !is.data.frame(probs)) {
    probs <- unlist(probs)
    labels <- unlist(labels)
  }

  probs <- as.numeric(probs)
  labels <- as.numeric(labels)

  # Compute metrics only - no fitting to prevent leakage
  ece_result <- compute_ece(probs, labels)
  brier_result <- decompose_brier(probs, labels)

  result <- list(
    ece = ece_result$ece,
    mce = ece_result$mce,
    brier = brier_result$brier,
    reliability = brier_result$reliability,
    resolution = brier_result$resolution,
    uncertainty = brier_result$uncertainty,
    n_samples = length(probs),
    bin_data = ece_result$bin_data
  )

  class(result) <- c("CalibrationResult", "list")
  result
}


#' @title Print method for CalibrationResult
#'
#' @param x CalibrationResult object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.CalibrationResult <- function(x, ...) {
  cat("=== Calibration Summary ===\n\n")
  cat(sprintf("Samples: %d\n\n", x$n_samples))

  cat("Calibration Metrics:\n")
  cat(sprintf("  ECE (Expected Calibration Error): %.4f\n", x$ece))
  cat(sprintf("  MCE (Maximum Calibration Error):  %.4f\n", x$mce))

  cat("\nBrier Score Decomposition:\n")
  cat(sprintf("  Total Brier Score:  %.4f\n", x$brier))
  cat(sprintf("  - Reliability:      %.4f (calibration, lower is better)\n", x$reliability))
  cat(sprintf("  - Resolution:       %.4f (discrimination, higher is better)\n", x$resolution))
  cat(sprintf("  - Uncertainty:      %.4f (inherent, fixed)\n", x$uncertainty))

  cat("\nTo fit calibrators, use fit_platt_scaling() or fit_isotonic_calibration()\n")
  cat("on a held-out calibration set (not the evaluation data).\n")

  invisible(x)
}
