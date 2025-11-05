# Test suite for clinical_utility.R
# Tests calibration, decision curve analysis, and clinical impact functions

library(testthat)

# ===== Test Data Setup =====

# Create test data with known calibration properties
create_calibration_test_data <- function(n = 300, calibrated = TRUE, seed = 123) {
  set.seed(seed)

  if (calibrated) {
    # Well-calibrated predictions
    true_probs <- runif(n, 0, 1)
    predictions <- true_probs + rnorm(n, 0, 0.1)
    predictions <- pmax(0.01, pmin(0.99, predictions))
    outcomes <- rbinom(n, 1, true_probs)
  } else {
    # Poorly calibrated (overconfident)
    true_probs <- runif(n, 0.3, 0.7)
    predictions <- ifelse(true_probs > 0.5,
                          runif(n, 0.7, 0.95),
                          runif(n, 0.05, 0.3))
    outcomes <- rbinom(n, 1, true_probs)
  }

  list(predictions = predictions, outcomes = outcomes)
}

calibrated_data <- create_calibration_test_data(calibrated = TRUE)
uncalibrated_data <- create_calibration_test_data(calibrated = FALSE)

# ===== Test 1: Calibration - Input Validation =====

test_that("OmicSelector_calibrate validates inputs correctly", {

  # Test missing predictions
  expect_error(
    OmicSelector_calibrate(
      predictions = NULL,
      outcomes = calibrated_data$outcomes
    ),
    "predictions"
  )

  # Test missing outcomes
  expect_error(
    OmicSelector_calibrate(
      predictions = calibrated_data$predictions,
      outcomes = NULL
    ),
    "outcomes"
  )

  # Test length mismatch
  expect_error(
    OmicSelector_calibrate(
      predictions = calibrated_data$predictions[1:100],
      outcomes = calibrated_data$outcomes
    ),
    "same length"
  )

  # Test predictions out of range
  bad_predictions <- c(calibrated_data$predictions[1:100], 1.5, -0.5)
  bad_outcomes <- c(calibrated_data$outcomes[1:100], 0, 1)

  expect_error(
    OmicSelector_calibrate(
      predictions = bad_predictions,
      outcomes = bad_outcomes
    ),
    "between 0 and 1"
  )

  # Test invalid method
  expect_error(
    OmicSelector_calibrate(
      predictions = calibrated_data$predictions,
      outcomes = calibrated_data$outcomes,
      method = "invalid_method"
    ),
    "should be one of"
  )
})

# ===== Test 2: Calibration - Platt Scaling =====

test_that("Platt scaling calibration works correctly", {

  result <- OmicSelector_calibrate(
    predictions = uncalibrated_data$predictions,
    outcomes = uncalibrated_data$outcomes,
    method = "platt",
    plot = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_calibration")

  # Check required elements
  expect_true("method" %in% names(result))
  expect_true("calibrated_predictions" %in% names(result))
  expect_true("metrics" %in% names(result))
  expect_true("calibration_table" %in% names(result))

  # Method should be recorded
  expect_equal(result$method, "platt")

  # Calibration model should exist
  expect_true(!is.null(result$calibration_model))

  # Calibrated predictions should be probabilities
  expect_true(all(result$calibrated_predictions >= 0))
  expect_true(all(result$calibrated_predictions <= 1))

  # Should have same length as input
  expect_equal(length(result$calibrated_predictions),
               length(uncalibrated_data$predictions))

  # Check metrics exist
  metrics <- result$metrics
  expect_true("brier_score" %in% names(metrics))
  expect_true("calibration_slope" %in% names(metrics))
  expect_true("calibration_intercept" %in% names(metrics))
  expect_true("ICI" %in% names(metrics))
})

# ===== Test 3: Calibration - Isotonic Regression =====

test_that("Isotonic regression calibration works", {

  result <- OmicSelector_calibrate(
    predictions = uncalibrated_data$predictions,
    outcomes = uncalibrated_data$outcomes,
    method = "isotonic",
    plot = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_calibration")

  # Calibrated predictions should be monotonic
  # (isotonic regression preserves ordering)
  order_original <- order(uncalibrated_data$predictions)
  calibrated_ordered <- result$calibrated_predictions[order_original]

  # Check if mostly monotonic (allowing small violations due to ties)
  diffs <- diff(calibrated_ordered)
  expect_true(sum(diffs >= -0.01) / length(diffs) > 0.95)
})

# ===== Test 4: Calibration - Beta Method =====

test_that("Beta calibration works", {

  result <- OmicSelector_calibrate(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    method = "beta",
    plot = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_calibration")

  # Should produce valid probabilities
  expect_true(all(result$calibrated_predictions >= 0))
  expect_true(all(result$calibrated_predictions <= 1))
})

# ===== Test 5: Calibration - Assessment Only =====

test_that("Assessment-only mode works", {

  result <- OmicSelector_calibrate(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    method = "assessment_only",
    plot = FALSE
  )

  # Calibrated predictions should equal original
  expect_equal(result$calibrated_predictions, result$original_predictions)

  # Should not have calibration model
  expect_null(result$calibration_model)

  # But should still have metrics
  expect_true(!is.null(result$metrics))
})

# ===== Test 6: Calibration Metrics =====

test_that("Calibration metrics are calculated correctly", {

  result <- OmicSelector_calibrate(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    method = "assessment_only",
    plot = FALSE
  )

  metrics <- result$metrics

  # Brier score: 0 to 1
  expect_true(metrics$brier_score >= 0 && metrics$brier_score <= 1)

  # Calibration slope: typically 0.5 to 2
  expect_true(metrics$calibration_slope > 0)

  # Hosmer-Lemeshow p-value: 0 to 1
  expect_true(metrics$hosmer_lemeshow_p >= 0 && metrics$hosmer_lemeshow_p <= 1)

  # E-statistics: 0 to 1
  expect_true(metrics$ICI >= 0 && metrics$ICI <= 1)
  expect_true(metrics$E50 >= 0 && metrics$E50 <= 1)
  expect_true(metrics$E90 >= 0 && metrics$E90 <= 1)
  expect_true(metrics$Emax >= 0 && metrics$Emax <= 1)

  # E-statistics ordering
  expect_true(metrics$E50 <= metrics$E90)
  expect_true(metrics$E90 <= metrics$Emax)
})

# ===== Test 7: Calibration Improves Predictions =====

test_that("Calibration improves poorly calibrated predictions", {

  # Original uncalibrated result
  result_before <- OmicSelector_calibrate(
    predictions = uncalibrated_data$predictions,
    outcomes = uncalibrated_data$outcomes,
    method = "assessment_only",
    plot = FALSE
  )

  # After Platt calibration
  result_after <- OmicSelector_calibrate(
    predictions = uncalibrated_data$predictions,
    outcomes = uncalibrated_data$outcomes,
    method = "platt",
    plot = FALSE
  )

  # Recalculate metrics for calibrated predictions
  result_improved <- OmicSelector_calibrate(
    predictions = result_after$calibrated_predictions,
    outcomes = uncalibrated_data$outcomes,
    method = "assessment_only",
    plot = FALSE
  )

  # Calibrated should have better or similar Brier score
  expect_true(result_improved$metrics$brier_score <= result_before$metrics$brier_score + 0.05)

  # Calibration slope should be closer to 1
  expect_true(abs(result_improved$metrics$calibration_slope - 1) <=
              abs(result_before$metrics$calibration_slope - 1))
})

# ===== Test 8: Decision Curve Analysis - Input Validation =====

test_that("OmicSelector_decision_curve validates inputs", {

  # Test invalid predictions
  expect_error(
    OmicSelector_decision_curve(
      predictions = NULL,
      outcomes = calibrated_data$outcomes
    ),
    "predictions"
  )

  # Test length mismatch
  expect_error(
    OmicSelector_decision_curve(
      predictions = calibrated_data$predictions[1:100],
      outcomes = calibrated_data$outcomes
    ),
    "same length"
  )

  # Test invalid thresholds
  expect_error(
    OmicSelector_decision_curve(
      predictions = calibrated_data$predictions,
      outcomes = calibrated_data$outcomes,
      thresholds = c(-0.1, 0.5, 1.5)
    ),
    "between 0 and 1"
  )
})

# ===== Test 9: Decision Curve Analysis - Basic Functionality =====

test_that("Decision curve analysis works correctly", {

  result <- OmicSelector_decision_curve(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    thresholds = seq(0.1, 0.9, by = 0.1),
    plot = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_dca")

  # Check required elements
  expect_true("net_benefit" %in% names(result))
  expect_true("optimal_threshold" %in% names(result))
  expect_true("prevalence" %in% names(result))

  # Net benefit data frame should have correct structure
  nb_df <- result$net_benefit
  expect_true("threshold" %in% names(nb_df))
  expect_true("model" %in% names(nb_df))
  expect_true("treat_all" %in% names(nb_df))
  expect_true("treat_none" %in% names(nb_df))

  # Treat none should always be 0
  expect_true(all(nb_df$treat_none == 0))

  # Optimal threshold should be within range
  expect_true(result$optimal_threshold >= min(seq(0.1, 0.9, by = 0.1)))
  expect_true(result$optimal_threshold <= max(seq(0.1, 0.9, by = 0.1)))
})

# ===== Test 10: DCA Net Benefit Properties =====

test_that("Net benefit has expected properties", {

  result <- OmicSelector_decision_curve(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    thresholds = seq(0.01, 0.99, by = 0.01),
    plot = FALSE
  )

  nb_df <- result$net_benefit

  # At very low thresholds, treat_all should be best
  low_threshold_rows <- nb_df[nb_df$threshold < 0.1, ]
  # (Not always true, but generally)

  # At very high thresholds, treat_none should be best
  high_threshold_rows <- nb_df[nb_df$threshold > 0.9, ]
  # Treat none net benefit is 0, treat all should be negative
  expect_true(mean(high_threshold_rows$treat_all) < 0)

  # Net benefit should be bounded
  expect_true(all(nb_df$model >= -1))  # Can't be too negative
  expect_true(all(nb_df$treat_all >= -1))
})

# ===== Test 11: DCA with Harm Adjustment =====

test_that("DCA works with harm adjustment", {

  result_no_harm <- OmicSelector_decision_curve(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    thresholds = seq(0.1, 0.9, by = 0.1),
    harm = NULL,
    plot = FALSE
  )

  result_with_harm <- OmicSelector_decision_curve(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    thresholds = seq(0.1, 0.9, by = 0.1),
    harm = 0.05,  # 5% harm
    plot = FALSE
  )

  # With harm, net benefits should be lower
  expect_true(mean(result_with_harm$net_benefit$model) <
              mean(result_no_harm$net_benefit$model))
})

# ===== Test 12: Clinical Impact - Input Validation =====

test_that("OmicSelector_clinical_impact validates inputs", {

  # Test invalid population size
  expect_error(
    OmicSelector_clinical_impact(
      predictions = calibrated_data$predictions,
      outcomes = calibrated_data$outcomes,
      population_size = -100
    ),
    "population_size"
  )

  # Test invalid thresholds
  expect_error(
    OmicSelector_clinical_impact(
      predictions = calibrated_data$predictions,
      outcomes = calibrated_data$outcomes,
      thresholds = c(-0.1, 0.5)
    ),
    "between 0 and 1"
  )
})

# ===== Test 13: Clinical Impact - Basic Functionality =====

test_that("Clinical impact analysis works correctly", {

  result <- OmicSelector_clinical_impact(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    thresholds = seq(0.1, 0.9, by = 0.1),
    population_size = 1000,
    plot = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_impact")

  # Check required elements
  expect_true("impact_table" %in% names(result))
  expect_true("optimal_threshold" %in% names(result))
  expect_true("population_size" %in% names(result))

  # Impact table should have correct structure
  impact_df <- result$impact_table
  expect_true("threshold" %in% names(impact_df))
  expect_true("n_high_risk" %in% names(impact_df))
  expect_true("n_true_positives" %in% names(impact_df))
  expect_true("n_false_positives" %in% names(impact_df))
  expect_true("ppv" %in% names(impact_df))
  expect_true("npv" %in% names(impact_df))
  expect_true("number_needed_to_test" %in% names(impact_df))
})

# ===== Test 14: Clinical Impact - Scaling =====

test_that("Clinical impact scales correctly with population size", {

  result_1000 <- OmicSelector_clinical_impact(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    thresholds = c(0.5),
    population_size = 1000,
    plot = FALSE
  )

  result_2000 <- OmicSelector_clinical_impact(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    thresholds = c(0.5),
    population_size = 2000,
    plot = FALSE
  )

  # Counts should scale proportionally
  expect_equal(
    result_2000$impact_table$n_high_risk[1] / result_1000$impact_table$n_high_risk[1],
    2,
    tolerance = 0.1
  )

  # But rates (PPV, NPV) should stay the same
  expect_equal(
    result_2000$impact_table$ppv[1],
    result_1000$impact_table$ppv[1],
    tolerance = 0.01
  )
})

# ===== Test 15: Clinical Impact - Threshold Effects =====

test_that("Clinical impact varies with threshold", {

  result <- OmicSelector_clinical_impact(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    thresholds = seq(0.1, 0.9, by = 0.2),
    population_size = 1000,
    plot = FALSE
  )

  impact_df <- result$impact_table

  # As threshold increases, fewer classified as high risk
  expect_true(all(diff(impact_df$n_high_risk) <= 0))

  # PPV should generally increase with threshold
  # (more stringent threshold = fewer false positives)
  ppv_trend <- diff(impact_df$ppv)
  expect_true(sum(ppv_trend > 0) >= sum(ppv_trend < 0))
})

# ===== Test 16: S3 Print Methods =====

test_that("S3 print methods work correctly", {

  # Calibration
  calib_result <- OmicSelector_calibrate(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    method = "platt",
    plot = FALSE
  )

  expect_output(print(calib_result), "Calibration")
  expect_output(print(calib_result), "Brier Score")

  # DCA
  dca_result <- OmicSelector_decision_curve(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    plot = FALSE
  )

  expect_output(print(dca_result), "Decision Curve")
  expect_output(print(dca_result), "Net benefit")

  # Clinical Impact
  impact_result <- OmicSelector_clinical_impact(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    plot = FALSE
  )

  expect_output(print(impact_result), "Clinical Impact")
  expect_output(print(impact_result), "Population")
})

# ===== Test 17: Plotting Functions =====

test_that("Plotting functions work", {
  skip_if_not_installed("ggplot2")

  # Calibration plot
  calib_result <- OmicSelector_calibrate(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    method = "platt",
    plot = TRUE
  )

  expect_true(!is.null(calib_result$plot))
  expect_s3_class(calib_result$plot, "ggplot")

  # DCA plot
  dca_result <- OmicSelector_decision_curve(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    plot = TRUE
  )

  expect_true(!is.null(dca_result$plot))
  expect_s3_class(dca_result$plot, "ggplot")

  # Clinical impact plot
  impact_result <- OmicSelector_clinical_impact(
    predictions = calibrated_data$predictions,
    outcomes = calibrated_data$outcomes,
    plot = TRUE
  )

  expect_true(!is.null(impact_result$plot))
  expect_s3_class(impact_result$plot, "ggplot")
})

# ===== Test 18: Integration Tests =====

test_that("Clinical utility functions work together", {

  # Step 1: Calibrate
  calib_result <- OmicSelector_calibrate(
    predictions = uncalibrated_data$predictions,
    outcomes = uncalibrated_data$outcomes,
    method = "platt",
    plot = FALSE
  )

  # Step 2: DCA with calibrated predictions
  dca_result <- OmicSelector_decision_curve(
    predictions = calib_result$calibrated_predictions,
    outcomes = uncalibrated_data$outcomes,
    plot = FALSE
  )

  # Step 3: Clinical impact with calibrated predictions
  impact_result <- OmicSelector_clinical_impact(
    predictions = calib_result$calibrated_predictions,
    outcomes = uncalibrated_data$outcomes,
    plot = FALSE
  )

  # All should succeed
  expect_s3_class(calib_result, "OmicSelector_calibration")
  expect_s3_class(dca_result, "OmicSelector_dca")
  expect_s3_class(impact_result, "OmicSelector_impact")

  # Optimal thresholds should be in reasonable range
  expect_true(dca_result$optimal_threshold > 0 && dca_result$optimal_threshold < 1)
  expect_true(impact_result$optimal_threshold > 0 && impact_result$optimal_threshold < 1)
})

# ===== Test 19: Edge Cases =====

test_that("Functions handle edge cases", {

  # Perfect predictions
  perfect_pred <- calibrated_data$outcomes
  perfect_outcomes <- calibrated_data$outcomes

  calib_perfect <- OmicSelector_calibrate(
    predictions = perfect_pred,
    outcomes = perfect_outcomes,
    method = "assessment_only",
    plot = FALSE
  )

  # Brier score should be 0
  expect_equal(calib_perfect$metrics$brier_score, 0, tolerance = 0.01)

  # All predictions same
  constant_pred <- rep(0.5, length(calibrated_data$outcomes))

  result_constant <- OmicSelector_decision_curve(
    predictions = constant_pred,
    outcomes = calibrated_data$outcomes,
    thresholds = seq(0.1, 0.9, by = 0.2),
    plot = FALSE
  )

  expect_s3_class(result_constant, "OmicSelector_dca")

  # Very small dataset
  small_pred <- calibrated_data$predictions[1:20]
  small_outcomes <- calibrated_data$outcomes[1:20]

  result_small <- OmicSelector_calibrate(
    predictions = small_pred,
    outcomes = small_outcomes,
    method = "platt",
    n_bins = 4,  # Fewer bins for small sample
    plot = FALSE
  )

  expect_s3_class(result_small, "OmicSelector_calibration")
})

# ===== Test 20: Factor Outcomes =====

test_that("Functions work with factor outcomes", {

  # Convert to factor
  outcomes_factor <- factor(calibrated_data$outcomes, levels = c(0, 1), labels = c("No", "Yes"))

  # Calibration
  calib_result <- OmicSelector_calibrate(
    predictions = calibrated_data$predictions,
    outcomes = outcomes_factor,
    method = "platt",
    plot = FALSE
  )

  expect_s3_class(calib_result, "OmicSelector_calibration")

  # DCA
  dca_result <- OmicSelector_decision_curve(
    predictions = calibrated_data$predictions,
    outcomes = outcomes_factor,
    plot = FALSE
  )

  expect_s3_class(dca_result, "OmicSelector_dca")

  # Clinical Impact
  impact_result <- OmicSelector_clinical_impact(
    predictions = calibrated_data$predictions,
    outcomes = outcomes_factor,
    plot = FALSE
  )

  expect_s3_class(impact_result, "OmicSelector_impact")
})
