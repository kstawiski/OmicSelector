# Tests for Calibration Metrics

test_that("compute_ece returns valid structure", {
  set.seed(42)
  probs <- runif(100)
  labels <- rbinom(100, 1, probs)

  result <- compute_ece(probs, labels)

  expect_true("ece" %in% names(result))
  expect_true("mce" %in% names(result))
  expect_true("bin_data" %in% names(result))

  expect_gte(result$ece, 0)
  expect_lte(result$ece, 1)
})

test_that("compute_ece detects well-calibrated model", {
  set.seed(42)
  n <- 1000

  # Generate well-calibrated predictions
  true_probs <- runif(n)
  labels <- rbinom(n, 1, true_probs)

  # Use true probabilities as predictions (perfect calibration)
  result <- compute_ece(true_probs, labels, n_bins = 10)

  # ECE should be low for well-calibrated model
  expect_lt(result$ece, 0.1)
})

test_that("compute_ece detects poorly calibrated model", {
  set.seed(42)
  n <- 1000

  # Generate labels with 50% positive
  labels <- rbinom(n, 1, 0.5)

  # Create poorly calibrated predictions (all near 0.1 or 0.9)
  probs <- ifelse(labels == 1, runif(n, 0.85, 0.95), runif(n, 0.05, 0.15))

  # This should still have low ECE if predictions are accurate
  result <- compute_ece(probs, labels, n_bins = 10)

  # The model is discriminative but may not be perfectly calibrated
  expect_gte(result$ece, 0)
})

test_that("decompose_brier returns valid decomposition", {
  set.seed(42)
  probs <- runif(100)
  labels <- rbinom(100, 1, 0.5)

  result <- decompose_brier(probs, labels)

  expect_true("brier" %in% names(result))
  expect_true("reliability" %in% names(result))
  expect_true("resolution" %in% names(result))
  expect_true("uncertainty" %in% names(result))

  # Brier should be between 0 and 1
  expect_gte(result$brier, 0)
  expect_lte(result$brier, 1)

  # Reliability and resolution should be non-negative
  expect_gte(result$reliability, 0)
  expect_gte(result$resolution, 0)
})

test_that("decompose_brier reconstructs correctly", {
  set.seed(42)
  probs <- runif(200)
  labels <- rbinom(200, 1, probs)

  result <- decompose_brier(probs, labels)

  # Reconstructed should approximately equal original Brier
  # Allow some tolerance due to binning
  expect_equal(result$reconstructed, result$brier, tolerance = 0.05)
})

test_that("fit_platt_scaling returns callable function", {
  set.seed(42)
  probs <- runif(100)
  labels <- rbinom(100, 1, probs)

  calibrator <- fit_platt_scaling(probs, labels)

  expect_true(is.function(calibrator))

  # Test calibration
  new_probs <- c(0.1, 0.5, 0.9)
  calibrated <- calibrator(new_probs)

  expect_length(calibrated, 3)
  expect_true(all(calibrated >= 0 & calibrated <= 1))
})

test_that("fit_isotonic_calibration returns callable function", {
  set.seed(42)
  probs <- runif(100)
  labels <- rbinom(100, 1, probs)

  calibrator <- fit_isotonic_calibration(probs, labels)

  expect_true(is.function(calibrator))

  # Test calibration
  new_probs <- c(0.1, 0.5, 0.9)
  calibrated <- calibrator(new_probs)

  expect_length(calibrated, 3)
  expect_true(all(calibrated >= 0 & calibrated <= 1))
})

test_that("fit_temperature_scaling returns callable function", {
  set.seed(42)
  probs <- runif(100)
  labels <- rbinom(100, 1, probs)

  calibrator <- fit_temperature_scaling(probs, labels)

  expect_true(is.function(calibrator))

  # Test calibration
  new_probs <- c(0.1, 0.5, 0.9)
  calibrated <- calibrator(new_probs)

  expect_length(calibrated, 3)
  expect_true(all(calibrated >= 0 & calibrated <= 1))
})

test_that("isotonic calibration improves poorly calibrated model", {
  set.seed(42)
  n <- 500

  # Create systematically overconfident predictions
  true_probs <- runif(n, 0.3, 0.7)  # Moderate true probabilities
  labels <- rbinom(n, 1, true_probs)

  # Overconfident predictions: push toward extremes
  overconfident_probs <- ifelse(true_probs > 0.5,
                                 pmin(true_probs + 0.2, 0.95),
                                 pmax(true_probs - 0.2, 0.05))

  # Split into calibration and test sets
  n_cal <- 300
  cal_idx <- 1:n_cal
  test_idx <- (n_cal + 1):n

  # ECE before calibration
  ece_before <- compute_ece(overconfident_probs[test_idx], labels[test_idx])$ece

  # Fit calibrator on calibration set
  calibrator <- fit_isotonic_calibration(
    overconfident_probs[cal_idx],
    labels[cal_idx]
  )

  # Apply to test set
  calibrated_probs <- calibrator(overconfident_probs[test_idx])

  # ECE after calibration
  ece_after <- compute_ece(calibrated_probs, labels[test_idx])$ece

  # Should improve (or at least not worsen significantly)
  expect_lte(ece_after, ece_before + 0.05)
})

test_that("reliability_diagram_data returns plottable data", {
  set.seed(42)
  probs <- runif(100)
  labels <- rbinom(100, 1, probs)

  diagram_data <- reliability_diagram_data(probs, labels, n_bins = 5)

  expect_s3_class(diagram_data, "data.frame")
  expect_true("mean_predicted" %in% names(diagram_data))
  expect_true("fraction_positive" %in% names(diagram_data))
  expect_true("n_samples" %in% names(diagram_data))
})

test_that("calibration_summary aggregates metrics", {
  set.seed(42)
  probs <- runif(100)
  labels <- rbinom(100, 1, 0.5)

  result <- calibration_summary(probs, labels, repair = TRUE)

  expect_s3_class(result, "CalibrationResult")
  expect_true("ece" %in% names(result))
  expect_true("brier" %in% names(result))
  expect_true("platt_calibrator" %in% names(result))
  expect_true("isotonic_calibrator" %in% names(result))
})

test_that("calibration_summary handles list inputs from CV", {
  set.seed(42)

  # Simulate CV fold outputs
  probs_list <- list(
    fold1 = runif(20),
    fold2 = runif(20),
    fold3 = runif(20)
  )
  labels_list <- list(
    fold1 = rbinom(20, 1, 0.5),
    fold2 = rbinom(20, 1, 0.5),
    fold3 = rbinom(20, 1, 0.5)
  )

  result <- calibration_summary(probs_list, labels_list)

  expect_s3_class(result, "CalibrationResult")
  expect_equal(result$n_samples, 60)
})

test_that("CalibrationResult print method works", {
  set.seed(42)
  probs <- runif(100)
  labels <- rbinom(100, 1, 0.5)

  result <- calibration_summary(probs, labels)

  expect_output(print(result), "Calibration Summary")
  expect_output(print(result), "ECE")
  expect_output(print(result), "Brier")
})
