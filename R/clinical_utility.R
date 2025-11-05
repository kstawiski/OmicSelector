#' @title Clinical Utility Metrics for Prediction Models
#' @description Functions for assessing clinical utility of prediction models
#'
#' This module provides tools for evaluating whether a prediction model is
#' clinically useful, beyond just statistical performance. Key features:
#'
#' - Model calibration assessment and correction
#' - Decision curve analysis (net benefit)
#' - Clinical impact curves
#' - Integration with TRIPOD+AI reporting
#'
#' @references
#' Van Calster, B., et al. (2016). Calibration: the Achilles heel of predictive analytics.
#'   BMC Medicine, 14(1), 230.
#'
#' Vickers, A. J., & Elkin, E. B. (2006). Decision curve analysis: a novel method for
#'   evaluating prediction models. Medical Decision Making, 26(6), 565-574.
#'
#' @author Konrad Stawiski
#' @name clinical_utility
NULL


#' Calibrate Prediction Model
#'
#' @description
#' Assesses and corrects calibration of predicted probabilities. Calibration
#' refers to how well predicted probabilities match observed outcome frequencies.
#' A well-calibrated model has predicted probabilities that match reality.
#'
#' @param predictions Numeric vector of predicted probabilities (0-1)
#' @param outcomes Factor or numeric vector of observed outcomes (0/1 or factor)
#' @param method Character string specifying calibration method:
#'   \itemize{
#'     \item "platt" - Platt scaling (logistic calibration)
#'     \item "isotonic" - Isotonic regression (non-parametric)
#'     \item "beta" - Beta calibration (for well-calibrated models)
#'     \item "assessment_only" - Only assess, don't calibrate
#'   }
#' @param n_bins Integer, number of bins for calibration plot (default: 10)
#' @param plot Logical, whether to generate calibration plots (default: TRUE)
#' @param conf_level Numeric, confidence level for intervals (default: 0.95)
#'
#' @return A list of class "OmicSelector_calibration" containing:
#'   \item{method}{Calibration method used}
#'   \item{calibration_model}{Fitted calibration model (if applicable)}
#'   \item{calibrated_predictions}{Calibrated probabilities}
#'   \item{metrics}{List of calibration metrics:}
#'     \itemize{
#'       \item brier_score - Overall accuracy of probabilistic predictions
#'       \item calibration_intercept - Calibration-in-the-large
#'       \item calibration_slope - Measure of overfitting (should be ~1)
#'       \item hosmer_lemeshow_statistic - Goodness of fit test
#'       \item hosmer_lemeshow_p - P-value (>0.05 is good)
#'       \item ICI - Integrated Calibration Index
#'       \item E50 - Median absolute calibration error
#'       \item E90 - 90th percentile calibration error
#'       \item Emax - Maximum absolute calibration error
#'     }
#'   \item{calibration_table}{Data frame with bins and observed vs expected}
#'   \item{plot}{ggplot2 calibration plot (if plot = TRUE)}
#'
#' @details
#' **Calibration Metrics Interpretation:**
#'
#' - **Brier Score**: Lower is better. <0.10 excellent, <0.25 acceptable
#' - **Calibration Slope**: Should be close to 1.0. <1 indicates overfitting
#' - **Calibration Intercept**: Should be close to 0
#' - **Hosmer-Lemeshow p-value**: >0.05 indicates good calibration
#' - **ICI**: <0.05 excellent, <0.10 acceptable
#' - **E50/E90/Emax**: Absolute errors, lower is better
#'
#' **When to use each method:**
#'
#' - **Platt scaling**: Most common, works well for most models
#' - **Isotonic**: For non-linear calibration issues
#' - **Beta**: For models that are already well-calibrated
#' - **Assessment only**: To evaluate without modifying predictions
#'
#' @examples
#' \dontrun{
#' # Simulate predictions and outcomes
#' set.seed(123)
#' predictions <- runif(200, 0, 1)
#' outcomes <- rbinom(200, 1, predictions)
#'
#' # Assess and calibrate
#' calib_result <- OmicSelector_calibrate(
#'   predictions = predictions,
#'   outcomes = outcomes,
#'   method = "platt",
#'   plot = TRUE
#' )
#'
#' # View metrics
#' print(calib_result)
#' calib_result$metrics
#'
#' # View plot
#' calib_result$plot
#'
#' # Use calibrated predictions
#' calibrated_probs <- calib_result$calibrated_predictions
#' }
#'
#' @export
OmicSelector_calibrate <- function(
    predictions,
    outcomes,
    method = c("platt", "isotonic", "beta", "assessment_only"),
    n_bins = 10,
    plot = TRUE,
    conf_level = 0.95
) {

  # Validate inputs
  method <- match.arg(method)
  .validate_calibration_inputs(predictions, outcomes, n_bins)

  # Convert outcomes to binary if factor
  if (is.factor(outcomes)) {
    outcomes_binary <- as.numeric(outcomes) - 1
  } else {
    outcomes_binary <- as.numeric(outcomes)
  }

  # Ensure predictions are probabilities
  if (any(predictions < 0) || any(predictions > 1)) {
    stop("Predictions must be probabilities between 0 and 1")
  }

  # Calculate calibration metrics
  metrics <- .calculate_calibration_metrics(predictions, outcomes_binary, n_bins)

  # Perform calibration if requested
  if (method != "assessment_only") {
    calibration_model <- .fit_calibration_model(predictions, outcomes_binary, method)
    calibrated_predictions <- .apply_calibration(predictions, calibration_model, method)
  } else {
    calibration_model <- NULL
    calibrated_predictions <- predictions
  }

  # Create calibration table for plotting
  calibration_table <- .create_calibration_table(predictions, outcomes_binary, n_bins, conf_level)

  # Generate plot if requested
  calibration_plot <- NULL
  if (plot) {
    calibration_plot <- .plot_calibration(calibration_table, method)
  }

  # Compile results
  result <- list(
    method = method,
    calibration_model = calibration_model,
    original_predictions = predictions,
    calibrated_predictions = calibrated_predictions,
    outcomes = outcomes_binary,
    metrics = metrics,
    calibration_table = calibration_table,
    plot = calibration_plot
  )

  class(result) <- c("OmicSelector_calibration", "list")
  return(result)
}


#' Decision Curve Analysis
#'
#' @description
#' Performs decision curve analysis to assess the clinical utility of a prediction
#' model across a range of threshold probabilities. DCA calculates the net benefit
#' of using the model to guide decisions compared to treating all or no patients.
#'
#' @param predictions Numeric vector of predicted probabilities (0-1)
#' @param outcomes Factor or numeric vector of observed outcomes
#' @param thresholds Numeric vector of threshold probabilities to evaluate.
#'   Default: seq(0.01, 0.99, by = 0.01)
#' @param plot Logical, whether to generate decision curve plot (default: TRUE)
#' @param harm Numeric, relative harm of intervention. Used to adjust net benefit
#'   for treatment-related harm. Default: NULL (no harm adjustment)
#'
#' @return A list of class "OmicSelector_dca" containing:
#'   \item{net_benefit}{Data frame with net benefit at each threshold for:}
#'     \itemize{
#'       \item model - Using the prediction model
#'       \item treat_all - Treating everyone
#'       \item treat_none - Treating no one
#'     }
#'   \item{standardized_net_benefit}{Net benefit standardized by prevalence}
#'   \item{optimal_threshold}{Threshold maximizing net benefit}
#'   \item{plot}{ggplot2 decision curve (if plot = TRUE)}
#'
#' @details
#' **Net Benefit Formula:**
#'
#' NB = (TP/n) - (FP/n) × (pt/(1-pt))
#'
#' Where:
#' - TP = True positives
#' - FP = False positives
#' - n = Sample size
#' - pt = Threshold probability
#'
#' **Interpretation:**
#'
#' - If model curve is above treat-all and treat-none, the model has clinical utility
#' - The range of thresholds where model is superior defines clinical usefulness
#' - Higher net benefit = more clinical utility
#'
#' @references
#' Vickers, A. J., & Elkin, E. B. (2006). Decision curve analysis: a novel method
#' for evaluating prediction models. Medical Decision Making, 26(6), 565-574.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' predictions <- runif(300, 0, 1)
#' outcomes <- rbinom(300, 1, predictions)
#'
#' # Perform DCA
#' dca_result <- OmicSelector_decision_curve(
#'   predictions = predictions,
#'   outcomes = outcomes,
#'   thresholds = seq(0.01, 0.99, by = 0.01)
#' )
#'
#' # View results
#' print(dca_result)
#' dca_result$plot
#'
#' # Find optimal threshold
#' dca_result$optimal_threshold
#' }
#'
#' @export
OmicSelector_decision_curve <- function(
    predictions,
    outcomes,
    thresholds = seq(0.01, 0.99, by = 0.01),
    plot = TRUE,
    harm = NULL
) {

  # Validate inputs
  .validate_dca_inputs(predictions, outcomes, thresholds)

  # Convert outcomes to binary
  if (is.factor(outcomes)) {
    outcomes_binary <- as.numeric(outcomes) - 1
  } else {
    outcomes_binary <- as.numeric(outcomes)
  }

  # Calculate net benefit for each threshold
  net_benefit_df <- .calculate_net_benefit(
    predictions,
    outcomes_binary,
    thresholds,
    harm
  )

  # Calculate standardized net benefit
  prevalence <- mean(outcomes_binary)
  standardized_nb <- net_benefit_df
  standardized_nb$model <- net_benefit_df$model / prevalence
  standardized_nb$treat_all <- net_benefit_df$treat_all / prevalence

  # Find optimal threshold (maximum net benefit)
  optimal_idx <- which.max(net_benefit_df$model)
  optimal_threshold <- thresholds[optimal_idx]
  optimal_nb <- net_benefit_df$model[optimal_idx]

  # Generate plot if requested
  dca_plot <- NULL
  if (plot) {
    dca_plot <- .plot_decision_curve(net_benefit_df, thresholds, optimal_threshold)
  }

  # Compile results
  result <- list(
    net_benefit = net_benefit_df,
    standardized_net_benefit = standardized_nb,
    optimal_threshold = optimal_threshold,
    optimal_net_benefit = optimal_nb,
    prevalence = prevalence,
    thresholds = thresholds,
    plot = dca_plot
  )

  class(result) <- c("OmicSelector_dca", "list")
  return(result)
}


#' Clinical Impact Curve
#'
#' @description
#' Creates a clinical impact curve showing the number of people classified as
#' high risk and the number of true positives at each threshold. This helps
#' understand the practical implications of using a model in clinical practice.
#'
#' @param predictions Numeric vector of predicted probabilities (0-1)
#' @param outcomes Factor or numeric vector of observed outcomes
#' @param thresholds Numeric vector of threshold probabilities. Default: seq(0.01, 0.99, by = 0.01)
#' @param population_size Integer, hypothetical population size for scaling results.
#'   Default: 1000
#' @param plot Logical, whether to generate impact curve plot (default: TRUE)
#'
#' @return A list of class "OmicSelector_impact" containing:
#'   \item{impact_table}{Data frame with for each threshold:}
#'     \itemize{
#'       \item threshold - Probability threshold
#'       \item n_high_risk - Number classified as high risk
#'       \item n_true_positives - Number with actual events among high risk
#'       \item n_false_positives - Number without events among high risk
#'       \item ppv - Positive predictive value at this threshold
#'       \item npv - Negative predictive value
#'       \item number_needed_to_test - NNT for one true positive
#'     }
#'   \item{optimal_threshold}{Threshold with best PPV/NPV trade-off}
#'   \item{plot}{ggplot2 clinical impact curve (if plot = TRUE)}
#'
#' @details
#' **Clinical Impact Interpretation:**
#'
#' The clinical impact curve shows:
#' - How many patients would be classified as high risk
#' - Of those, how many actually have the outcome (true positives)
#' - Number needed to test/treat to find one case
#'
#' This is useful for:
#' - Planning resource allocation
#' - Understanding workload implications
#' - Communicating model utility to clinicians
#' - Choosing appropriate thresholds for clinical use
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' predictions <- runif(500, 0, 1)
#' outcomes <- rbinom(500, 1, predictions)
#'
#' # Create clinical impact curve
#' impact_result <- OmicSelector_clinical_impact(
#'   predictions = predictions,
#'   outcomes = outcomes,
#'   population_size = 1000
#' )
#'
#' # View results
#' print(impact_result)
#' impact_result$plot
#'
#' # At 20% threshold, how many high risk?
#' impact_result$impact_table[impact_result$impact_table$threshold == 0.20, ]
#' }
#'
#' @export
OmicSelector_clinical_impact <- function(
    predictions,
    outcomes,
    thresholds = seq(0.01, 0.99, by = 0.01),
    population_size = 1000,
    plot = TRUE
) {

  # Validate inputs
  .validate_impact_inputs(predictions, outcomes, thresholds, population_size)

  # Convert outcomes to binary
  if (is.factor(outcomes)) {
    outcomes_binary <- as.numeric(outcomes) - 1
  } else {
    outcomes_binary <- as.numeric(outcomes)
  }

  # Calculate impact metrics for each threshold
  impact_table <- .calculate_clinical_impact(
    predictions,
    outcomes_binary,
    thresholds,
    population_size
  )

  # Find optimal threshold (best balance of PPV and NPV)
  # Using F1 score as criterion
  impact_table$f1_score <- 2 * (impact_table$ppv * impact_table$sensitivity) /
    (impact_table$ppv + impact_table$sensitivity)
  optimal_idx <- which.max(impact_table$f1_score)
  optimal_threshold <- thresholds[optimal_idx]

  # Generate plot if requested
  impact_plot <- NULL
  if (plot) {
    impact_plot <- .plot_clinical_impact(impact_table, optimal_threshold)
  }

  # Compile results
  result <- list(
    impact_table = impact_table,
    optimal_threshold = optimal_threshold,
    population_size = population_size,
    plot = impact_plot
  )

  class(result) <- c("OmicSelector_impact", "list")
  return(result)
}


# ===== S3 Methods =====

#' Print method for calibration results
#' @param x An OmicSelector_calibration object
#' @param ... Additional arguments (not used)
#' @export
print.OmicSelector_calibration <- function(x, ...) {
  cat("OmicSelector Model Calibration\n")
  cat("==============================\n\n")
  cat("Method:", x$method, "\n\n")

  cat("Calibration Metrics:\n")
  cat(sprintf("  Brier Score:             %.4f %s\n",
              x$metrics$brier_score,
              ifelse(x$metrics$brier_score < 0.10, "(Excellent)",
                     ifelse(x$metrics$brier_score < 0.25, "(Acceptable)", "(Poor)"))))
  cat(sprintf("  Calibration Intercept:   %.4f %s\n",
              x$metrics$calibration_intercept,
              ifelse(abs(x$metrics$calibration_intercept) < 0.10, "(Good)", "(Poor)"))))
  cat(sprintf("  Calibration Slope:       %.4f %s\n",
              x$metrics$calibration_slope,
              ifelse(abs(x$metrics$calibration_slope - 1) < 0.10, "(Good)",
                     ifelse(x$metrics$calibration_slope < 1, "(Overfitting)", "(Underfitting)"))))
  cat(sprintf("  Hosmer-Lemeshow p-value: %.4f %s\n",
              x$metrics$hosmer_lemeshow_p,
              ifelse(x$metrics$hosmer_lemeshow_p > 0.05, "(Good fit)", "(Poor fit)"))))

  cat("\nCalibration Error Metrics:\n")
  cat(sprintf("  ICI (Integrated):        %.4f\n", x$metrics$ICI))
  cat(sprintf("  E50 (Median error):      %.4f\n", x$metrics$E50))
  cat(sprintf("  E90 (90th percentile):   %.4f\n", x$metrics$E90))
  cat(sprintf("  Emax (Maximum error):    %.4f\n", x$metrics$Emax))

  if (!is.null(x$plot)) {
    cat("\nCalibration plot available in $plot\n")
  }

  invisible(x)
}

#' Print method for decision curve analysis
#' @param x An OmicSelector_dca object
#' @param ... Additional arguments (not used)
#' @export
print.OmicSelector_dca <- function(x, ...) {
  cat("OmicSelector Decision Curve Analysis\n")
  cat("====================================\n\n")

  cat(sprintf("Outcome prevalence: %.2f%%\n", x$prevalence * 100))
  cat(sprintf("Number of thresholds evaluated: %d\n", length(x$thresholds)))
  cat(sprintf("Threshold range: %.2f to %.2f\n\n",
              min(x$thresholds), max(x$thresholds)))

  cat(sprintf("Optimal threshold: %.3f\n", x$optimal_threshold))
  cat(sprintf("Net benefit at optimal: %.4f\n\n", x$optimal_net_benefit))

  # Find range where model beats both strategies
  nb_df <- x$net_benefit
  useful_range <- nb_df$model > nb_df$treat_all & nb_df$model > nb_df$treat_none
  if (any(useful_range)) {
    useful_thresholds <- x$thresholds[useful_range]
    cat(sprintf("Model shows clinical utility for thresholds: %.2f to %.2f\n",
                min(useful_thresholds), max(useful_thresholds)))
  } else {
    cat("Warning: Model does not show clinical utility in any threshold range\n")
  }

  if (!is.null(x$plot)) {
    cat("\nDecision curve plot available in $plot\n")
  }

  invisible(x)
}

#' Print method for clinical impact
#' @param x An OmicSelector_impact object
#' @param ... Additional arguments (not used)
#' @export
print.OmicSelector_impact <- function(x, ...) {
  cat("OmicSelector Clinical Impact Analysis\n")
  cat("=====================================\n\n")

  cat(sprintf("Population size: %d\n", x$population_size))
  cat(sprintf("Optimal threshold: %.3f\n\n", x$optimal_threshold))

  # Get metrics at optimal threshold
  optimal_row <- x$impact_table[which.min(abs(x$impact_table$threshold - x$optimal_threshold)), ]

  cat("At optimal threshold:\n")
  cat(sprintf("  Classified as high risk: %d (%.1f%%)\n",
              optimal_row$n_high_risk,
              optimal_row$n_high_risk / x$population_size * 100))
  cat(sprintf("  True positives:          %d\n", optimal_row$n_true_positives))
  cat(sprintf("  False positives:         %d\n", optimal_row$n_false_positives))
  cat(sprintf("  PPV:                     %.2f%%\n", optimal_row$ppv * 100))
  cat(sprintf("  NPV:                     %.2f%%\n", optimal_row$npv * 100))
  cat(sprintf("  Number needed to test:   %.1f\n", optimal_row$number_needed_to_test))

  if (!is.null(x$plot)) {
    cat("\nClinical impact plot available in $plot\n")
  }

  invisible(x)
}


# ===== Internal Helper Functions =====

#' Validate calibration inputs
#' @keywords internal
.validate_calibration_inputs <- function(predictions, outcomes, n_bins) {
  if (!is.numeric(predictions)) {
    stop("predictions must be numeric")
  }
  if (length(predictions) != length(outcomes)) {
    stop("predictions and outcomes must have the same length")
  }
  if (any(is.na(predictions)) || any(is.na(outcomes))) {
    stop("predictions and outcomes cannot contain NA values")
  }
  if (n_bins < 2 || n_bins > length(predictions) / 5) {
    stop("n_bins must be between 2 and sample_size/5")
  }
}

#' Validate DCA inputs
#' @keywords internal
.validate_dca_inputs <- function(predictions, outcomes, thresholds) {
  if (!is.numeric(predictions)) {
    stop("predictions must be numeric")
  }
  if (length(predictions) != length(outcomes)) {
    stop("predictions and outcomes must have the same length")
  }
  if (any(thresholds < 0) || any(thresholds > 1)) {
    stop("thresholds must be between 0 and 1")
  }
}

#' Validate impact inputs
#' @keywords internal
.validate_impact_inputs <- function(predictions, outcomes, thresholds, population_size) {
  .validate_dca_inputs(predictions, outcomes, thresholds)
  if (!is.numeric(population_size) || population_size < 1) {
    stop("population_size must be a positive integer")
  }
}

#' Calculate calibration metrics
#' @keywords internal
.calculate_calibration_metrics <- function(predictions, outcomes, n_bins) {
  # Brier score
  brier_score <- mean((predictions - outcomes)^2)

  # Calibration intercept and slope
  cal_model <- glm(outcomes ~ predictions, family = binomial())
  cal_intercept <- coef(cal_model)[1]
  cal_slope <- coef(cal_model)[2]

  # Hosmer-Lemeshow test
  hl_result <- .hosmer_lemeshow_test(predictions, outcomes, n_bins)

  # E-statistics (calibration errors)
  bin_stats <- .create_calibration_table(predictions, outcomes, n_bins, 0.95)
  errors <- abs(bin_stats$observed_rate - bin_stats$predicted_rate)

  ICI <- mean(errors, na.rm = TRUE)
  E50 <- median(errors, na.rm = TRUE)
  E90 <- quantile(errors, 0.90, na.rm = TRUE)
  Emax <- max(errors, na.rm = TRUE)

  list(
    brier_score = brier_score,
    calibration_intercept = cal_intercept,
    calibration_slope = cal_slope,
    hosmer_lemeshow_statistic = hl_result$statistic,
    hosmer_lemeshow_p = hl_result$p_value,
    ICI = ICI,
    E50 = E50,
    E90 = E90,
    Emax = Emax
  )
}

#' Hosmer-Lemeshow goodness of fit test
#' @keywords internal
.hosmer_lemeshow_test <- function(predictions, outcomes, n_bins) {
  # Create bins
  breaks <- quantile(predictions, probs = seq(0, 1, length.out = n_bins + 1))
  bins <- cut(predictions, breaks = breaks, include.lowest = TRUE)

  # Calculate observed and expected
  obs <- tapply(outcomes, bins, sum)
  exp <- tapply(predictions, bins, sum)
  n <- tapply(outcomes, bins, length)

  # Chi-square statistic
  chi_sq <- sum((obs - exp)^2 / (exp * (1 - exp / n)), na.rm = TRUE)
  df <- n_bins - 2
  p_value <- pchisq(chi_sq, df, lower.tail = FALSE)

  list(statistic = chi_sq, p_value = p_value, df = df)
}

#' Create calibration table for plotting
#' @keywords internal
.create_calibration_table <- function(predictions, outcomes, n_bins, conf_level) {
  # Create quantile-based bins
  breaks <- quantile(predictions, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  bins <- cut(predictions, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  # Calculate statistics for each bin
  bin_stats <- data.frame(
    bin = 1:n_bins,
    n = as.numeric(tapply(outcomes, bins, length)),
    observed = as.numeric(tapply(outcomes, bins, sum)),
    predicted = as.numeric(tapply(predictions, bins, sum))
  )

  bin_stats$observed_rate <- bin_stats$observed / bin_stats$n
  bin_stats$predicted_rate <- bin_stats$predicted / bin_stats$n

  # Confidence intervals for observed rate (Wilson score interval)
  alpha <- 1 - conf_level
  z <- qnorm(1 - alpha/2)

  bin_stats$ci_lower <- with(bin_stats, {
    p_hat <- observed_rate
    (p_hat + z^2/(2*n) - z * sqrt(p_hat*(1-p_hat)/n + z^2/(4*n^2))) / (1 + z^2/n)
  })

  bin_stats$ci_upper <- with(bin_stats, {
    p_hat <- observed_rate
    (p_hat + z^2/(2*n) + z * sqrt(p_hat*(1-p_hat)/n + z^2/(4*n^2))) / (1 + z^2/n)
  })

  return(bin_stats)
}

#' Fit calibration model
#' @keywords internal
.fit_calibration_model <- function(predictions, outcomes, method) {
  if (method == "platt") {
    # Logistic calibration (Platt scaling)
    logit_pred <- log(predictions / (1 - predictions + 1e-10))
    model <- glm(outcomes ~ logit_pred, family = binomial())
    return(model)

  } else if (method == "isotonic") {
    # Isotonic regression
    if (!requireNamespace("stats", quietly = TRUE)) {
      stop("stats package required for isotonic regression")
    }
    order_idx <- order(predictions)
    iso_model <- isoreg(predictions[order_idx], outcomes[order_idx])
    return(iso_model)

  } else if (method == "beta") {
    # Beta calibration (simplified)
    # Fit beta distribution parameters
    model <- glm(outcomes ~ predictions, family = binomial())
    return(model)
  }
}

#' Apply calibration model
#' @keywords internal
.apply_calibration <- function(predictions, model, method) {
  if (method == "platt") {
    logit_pred <- log(predictions / (1 - predictions + 1e-10))
    calibrated <- predict(model, newdata = data.frame(logit_pred = logit_pred), type = "response")

  } else if (method == "isotonic") {
    # Interpolate using isotonic regression
    calibrated <- approx(model$x, model$yf, xout = predictions, rule = 2)$y

  } else if (method == "beta") {
    calibrated <- predict(model, newdata = data.frame(predictions = predictions), type = "response")
  }

  # Ensure in [0,1]
  calibrated <- pmax(0, pmin(1, calibrated))
  return(calibrated)
}

#' Plot calibration curve
#' @keywords internal
.plot_calibration <- function(calibration_table, method) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 required for plotting. Returning NULL.")
    return(NULL)
  }

  p <- ggplot2::ggplot(calibration_table, ggplot2::aes(x = predicted_rate, y = observed_rate)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.6, color = "steelblue") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, alpha = 0.5) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 1) +
    ggplot2::scale_size_continuous(name = "N per bin") +
    ggplot2::labs(
      title = paste("Calibration Plot -", tools::toTitleCase(method)),
      subtitle = "Dashed line = perfect calibration",
      x = "Predicted Probability",
      y = "Observed Frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))

  return(p)
}

#' Calculate net benefit for decision curve
#' @keywords internal
.calculate_net_benefit <- function(predictions, outcomes, thresholds, harm) {
  n <- length(outcomes)

  nb_model <- numeric(length(thresholds))
  nb_all <- numeric(length(thresholds))
  nb_none <- rep(0, length(thresholds))

  for (i in seq_along(thresholds)) {
    pt <- thresholds[i]

    # Model strategy
    test_positive <- predictions >= pt
    tp <- sum(test_positive & outcomes == 1)
    fp <- sum(test_positive & outcomes == 0)

    nb_model[i] <- (tp / n) - (fp / n) * (pt / (1 - pt))

    # Treat all strategy
    prevalence <- mean(outcomes)
    nb_all[i] <- prevalence - (1 - prevalence) * (pt / (1 - pt))

    # Adjust for harm if specified
    if (!is.null(harm)) {
      nb_model[i] <- nb_model[i] - harm * sum(test_positive) / n
      nb_all[i] <- nb_all[i] - harm
    }
  }

  data.frame(
    threshold = thresholds,
    model = nb_model,
    treat_all = nb_all,
    treat_none = nb_none
  )
}

#' Plot decision curve
#' @keywords internal
.plot_decision_curve <- function(net_benefit_df, thresholds, optimal_threshold) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 required for plotting. Returning NULL.")
    return(NULL)
  }

  # Reshape for plotting
  nb_long <- data.frame(
    threshold = rep(thresholds, 3),
    net_benefit = c(net_benefit_df$model, net_benefit_df$treat_all, net_benefit_df$treat_none),
    strategy = rep(c("Model", "Treat All", "Treat None"), each = length(thresholds))
  )

  p <- ggplot2::ggplot(nb_long, ggplot2::aes(x = threshold, y = net_benefit, color = strategy)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_vline(xintercept = optimal_threshold, linetype = "dashed", alpha = 0.5) +
    ggplot2::scale_color_manual(values = c("Model" = "steelblue", "Treat All" = "red", "Treat None" = "gray50")) +
    ggplot2::labs(
      title = "Decision Curve Analysis",
      subtitle = paste("Optimal threshold:", round(optimal_threshold, 3)),
      x = "Threshold Probability",
      y = "Net Benefit",
      color = "Strategy"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  return(p)
}

#' Calculate clinical impact metrics
#' @keywords internal
.calculate_clinical_impact <- function(predictions, outcomes, thresholds, population_size) {
  n <- length(outcomes)
  scale_factor <- population_size / n

  impact_list <- lapply(thresholds, function(pt) {
    # Classify as high risk if prediction >= threshold
    high_risk <- predictions >= pt

    # Calculate counts (scaled to population)
    n_high_risk <- sum(high_risk) * scale_factor
    n_true_pos <- sum(high_risk & outcomes == 1) * scale_factor
    n_false_pos <- sum(high_risk & outcomes == 0) * scale_factor
    n_true_neg <- sum(!high_risk & outcomes == 0) * scale_factor
    n_false_neg <- sum(!high_risk & outcomes == 1) * scale_factor

    # Calculate metrics
    ppv <- ifelse(n_high_risk > 0, n_true_pos / n_high_risk, NA)
    npv <- ifelse(population_size - n_high_risk > 0,
                  n_true_neg / (population_size - n_high_risk), NA)
    sensitivity <- n_true_pos / (n_true_pos + n_false_neg)
    specificity <- n_true_neg / (n_true_neg + n_false_pos)
    nnt <- ifelse(ppv > 0, 1 / ppv, Inf)

    data.frame(
      threshold = pt,
      n_high_risk = round(n_high_risk),
      n_true_positives = round(n_true_pos),
      n_false_positives = round(n_false_pos),
      ppv = ppv,
      npv = npv,
      sensitivity = sensitivity,
      specificity = specificity,
      number_needed_to_test = nnt
    )
  })

  do.call(rbind, impact_list)
}

#' Plot clinical impact curve
#' @keywords internal
.plot_clinical_impact <- function(impact_table, optimal_threshold) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 required for plotting. Returning NULL.")
    return(NULL)
  }

  # Create plot with two y-axes concept
  p <- ggplot2::ggplot(impact_table, ggplot2::aes(x = threshold)) +
    ggplot2::geom_line(ggplot2::aes(y = n_high_risk, color = "High Risk"), linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(y = n_true_positives, color = "True Positives"), linewidth = 1) +
    ggplot2::geom_vline(xintercept = optimal_threshold, linetype = "dashed", alpha = 0.5) +
    ggplot2::scale_color_manual(values = c("High Risk" = "steelblue", "True Positives" = "darkgreen")) +
    ggplot2::labs(
      title = "Clinical Impact Curve",
      subtitle = paste("Optimal threshold:", round(optimal_threshold, 3)),
      x = "Threshold Probability",
      y = "Number of Patients",
      color = "Classification"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  return(p)
}
