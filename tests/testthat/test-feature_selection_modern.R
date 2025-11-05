# Test suite for feature_selection_modern.R
# Tests OmicSelector_stable_features() and related functions

library(testthat)

# ===== Test Data Setup =====

# Create synthetic test data
create_test_data <- function(n_samples = 100, n_features = 50, n_informative = 10, seed = 123) {
  set.seed(seed)

  # Create informative features
  X_informative <- matrix(rnorm(n_samples * n_informative), ncol = n_informative)
  colnames(X_informative) <- paste0("important_", 1:n_informative)

  # Create noise features
  X_noise <- matrix(rnorm(n_samples * (n_features - n_informative)),
                    ncol = n_features - n_informative)
  colnames(X_noise) <- paste0("noise_", 1:(n_features - n_informative))

  # Combine
  X <- cbind(X_informative, X_noise)

  # Create outcome based on first 5 informative features
  y <- rowSums(X_informative[, 1:5]) + rnorm(n_samples, sd = 0.5)
  outcome <- ifelse(y > median(y), "High", "Low")
  outcome <- factor(outcome)

  # Create data frame
  data <- data.frame(X, outcome = outcome)

  return(data)
}

test_data <- create_test_data()

# ===== Test 1: Input Validation =====

test_that("OmicSelector_stable_features validates inputs correctly", {

  # Test missing data
  expect_error(
    OmicSelector_stable_features(
      data = NULL,
      outcome = "outcome"
    ),
    "data"
  )

  # Test missing outcome
  expect_error(
    OmicSelector_stable_features(
      data = test_data,
      outcome = "nonexistent"
    ),
    "outcome"
  )

  # Test invalid method
  expect_error(
    OmicSelector_stable_features(
      data = test_data,
      outcome = "outcome",
      method = "invalid_method"
    ),
    "should be one of"
  )

  # Test invalid selection_threshold
  expect_error(
    OmicSelector_stable_features(
      data = test_data,
      outcome = "outcome",
      selection_threshold = 1.5
    ),
    "selection_threshold"
  )

  # Test invalid subsample_rate
  expect_error(
    OmicSelector_stable_features(
      data = test_data,
      outcome = "outcome",
      subsample_rate = 0
    ),
    "subsample_rate"
  )
})

# ===== Test 2: Stability Selection Method =====

test_that("Stability selection identifies stable features", {
  skip_if_not_installed("glmnet")

  result <- OmicSelector_stable_features(
    data = test_data,
    outcome = "outcome",
    method = "stability_selection",
    n_iterations = 20,  # Reduced for testing speed
    selection_threshold = 0.5,
    parallel = FALSE,
    verbose = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_stable_features")

  # Check required elements
  expect_true("selected_features" %in% names(result))
  expect_true("stability_scores" %in% names(result))
  expect_true("selection_frequency" %in% names(result))
  expect_true("nogueira_metrics" %in% names(result))
  expect_true("method" %in% names(result))

  # Check selected features
  expect_type(result$selected_features, "character")
  expect_true(length(result$selected_features) > 0)

  # Check stability scores
  expect_type(result$stability_scores, "double")
  expect_true(all(result$stability_scores >= 0 & result$stability_scores <= 1))

  # Check that informative features are preferentially selected
  important_selected <- sum(grepl("important_", result$selected_features))
  noise_selected <- sum(grepl("noise_", result$selected_features))
  expect_true(important_selected >= noise_selected)

  # Check Nogueira metrics exist
  expect_true("stability" %in% names(result$nogueira_metrics))
  expect_true("kuncheva_index" %in% names(result$nogueira_metrics))
  expect_true("mean_jaccard" %in% names(result$nogueira_metrics))
})

# ===== Test 3: Boruta Stable Method =====

test_that("Boruta stable method works correctly", {
  skip_if_not_installed("Boruta")
  skip_if_not_installed("randomForest")

  result <- OmicSelector_stable_features(
    data = test_data,
    outcome = "outcome",
    method = "boruta_stable",
    n_iterations = 10,  # Reduced for testing
    selection_threshold = 0.5,
    parallel = FALSE,
    verbose = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_stable_features")

  # Check structure
  expect_true(length(result$selected_features) > 0)
  expect_true(all(result$selected_features %in% colnames(test_data)))

  # Method should be recorded
  expect_equal(result$method, "boruta_stable")

  # Check parameters are stored
  expect_true("parameters" %in% names(result))
  expect_equal(result$parameters$method, "boruta_stable")
})

# ===== Test 4: LASSO Stable Method =====

test_that("LASSO stable method works correctly", {
  skip_if_not_installed("glmnet")

  result <- OmicSelector_stable_features(
    data = test_data,
    outcome = "outcome",
    method = "lasso_stable",
    n_iterations = 20,
    selection_threshold = 0.6,
    parallel = FALSE,
    verbose = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_stable_features")

  # Check structure
  expect_true(length(result$selected_features) > 0)

  # Check stability scores are calculated
  expect_true(all(result$stability_scores >= 0))
  expect_true(all(result$stability_scores <= 1))

  # Selected features should have high stability
  selected_scores <- result$stability_scores[result$selected_features]
  expect_true(all(selected_scores >= result$parameters$selection_threshold))
})

# ===== Test 5: RFE Stable Method =====

test_that("RFE stable method works correctly", {
  skip_if_not_installed("caret")
  skip_if_not_installed("randomForest")

  result <- OmicSelector_stable_features(
    data = test_data,
    outcome = "outcome",
    method = "rfe_stable",
    n_iterations = 10,
    selection_threshold = 0.5,
    parallel = FALSE,
    verbose = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_stable_features")

  # Check structure
  expect_true(length(result$selected_features) > 0)
  expect_type(result$stability_scores, "double")
})

# ===== Test 6: Nogueira Stability Metrics =====

test_that("Nogueira stability metrics are calculated correctly", {
  skip_if_not_installed("glmnet")

  result <- OmicSelector_stable_features(
    data = test_data,
    outcome = "outcome",
    method = "stability_selection",
    n_iterations = 30,
    selection_threshold = 0.6,
    parallel = FALSE,
    verbose = FALSE
  )

  metrics <- result$nogueira_metrics

  # Check all required metrics exist
  expect_true("stability" %in% names(metrics))
  expect_true("kuncheva_index" %in% names(metrics))
  expect_true("mean_jaccard" %in% names(metrics))
  expect_true("median_jaccard" %in% names(metrics))
  expect_true("avg_n_features" %in% names(metrics))
  expect_true("variance_n_features" %in% names(metrics))

  # Check value ranges
  # Kuncheva index: -1 to 1
  expect_true(metrics$kuncheva_index >= -1 && metrics$kuncheva_index <= 1)

  # Jaccard: 0 to 1
  expect_true(metrics$mean_jaccard >= 0 && metrics$mean_jaccard <= 1)
  expect_true(metrics$median_jaccard >= 0 && metrics$median_jaccard <= 1)

  # Average number of features should be positive
  expect_true(metrics$avg_n_features > 0)

  # Variance should be non-negative
  expect_true(metrics$variance_n_features >= 0)
})

# ===== Test 7: Selection Threshold Effect =====

test_that("Selection threshold affects number of selected features", {
  skip_if_not_installed("glmnet")

  # Lower threshold should select more features
  result_low <- OmicSelector_stable_features(
    data = test_data,
    outcome = "outcome",
    method = "stability_selection",
    n_iterations = 20,
    selection_threshold = 0.3,
    parallel = FALSE,
    verbose = FALSE
  )

  # Higher threshold should select fewer features
  result_high <- OmicSelector_stable_features(
    data = test_data,
    outcome = "outcome",
    method = "stability_selection",
    n_iterations = 20,
    selection_threshold = 0.8,
    seed = 123,  # Same seed for comparison
    parallel = FALSE,
    verbose = FALSE
  )

  # Lower threshold should give more features
  expect_true(length(result_low$selected_features) >=
              length(result_high$selected_features))
})

# ===== Test 8: Max Features Parameter =====

test_that("Max features parameter limits selection", {
  skip_if_not_installed("glmnet")

  max_f <- 10

  result <- OmicSelector_stable_features(
    data = test_data,
    outcome = "outcome",
    method = "stability_selection",
    n_iterations = 20,
    selection_threshold = 0.3,  # Low threshold to try to select many
    max_features = max_f,
    parallel = FALSE,
    verbose = FALSE
  )

  # Should not exceed max_features
  expect_true(length(result$selected_features) <= max_f)
})

# ===== Test 9: Reproducibility with Seed =====

test_that("Results are reproducible with seed", {
  skip_if_not_installed("glmnet")

  result1 <- OmicSelector_stable_features(
    data = test_data,
    outcome = "outcome",
    method = "stability_selection",
    n_iterations = 20,
    seed = 42,
    parallel = FALSE,
    verbose = FALSE
  )

  result2 <- OmicSelector_stable_features(
    data = test_data,
    outcome = "outcome",
    method = "stability_selection",
    n_iterations = 20,
    seed = 42,
    parallel = FALSE,
    verbose = FALSE
  )

  # Should get identical features
  expect_identical(result1$selected_features, result2$selected_features)

  # Should get identical stability scores
  expect_equal(result1$stability_scores, result2$stability_scores)
})

# ===== Test 10: S3 Methods =====

test_that("S3 methods work correctly", {
  skip_if_not_installed("glmnet")

  result <- OmicSelector_stable_features(
    data = test_data,
    outcome = "outcome",
    method = "stability_selection",
    n_iterations = 20,
    parallel = FALSE,
    verbose = FALSE
  )

  # Test print method
  expect_output(print(result), "OmicSelector")
  expect_output(print(result), "Stable Feature Selection")

  # Test summary method (if exists)
  if ("summary.OmicSelector_stable_features" %in% methods("summary")) {
    expect_output(summary(result), "Features")
  }

  # Test plot method (if exists)
  if ("plot.OmicSelector_stable_features" %in% methods("plot")) {
    expect_silent(p <- plot(result))
    expect_true(inherits(p, "ggplot") || inherits(p, "recordedplot"))
  }
})

# ===== Test 11: Edge Cases =====

test_that("Handles edge cases gracefully", {
  skip_if_not_installed("glmnet")

  # Very small dataset
  small_data <- test_data[1:20, 1:10]

  expect_warning(
    result <- OmicSelector_stable_features(
      data = small_data,
      outcome = "outcome",
      method = "stability_selection",
      n_iterations = 5,
      parallel = FALSE,
      verbose = FALSE
    ),
    NA  # No warning expected, but handle if one occurs
  )

  # Dataset with high correlation (multicollinearity)
  corr_data <- test_data
  corr_data$corr_1 <- test_data$important_1 + rnorm(nrow(test_data), sd = 0.1)
  corr_data$corr_2 <- test_data$important_1 + rnorm(nrow(test_data), sd = 0.1)

  result_corr <- OmicSelector_stable_features(
    data = corr_data,
    outcome = "outcome",
    method = "stability_selection",
    n_iterations = 10,
    parallel = FALSE,
    verbose = FALSE
  )

  expect_s3_class(result_corr, "OmicSelector_stable_features")
})

# ===== Test 12: Integration with Other Functions =====

test_that("Results can be used in downstream analysis", {
  skip_if_not_installed("glmnet")

  result <- OmicSelector_stable_features(
    data = test_data,
    outcome = "outcome",
    method = "stability_selection",
    n_iterations = 20,
    parallel = FALSE,
    verbose = FALSE
  )

  # Selected features should be valid column names
  expect_true(all(result$selected_features %in% colnames(test_data)))

  # Should be able to subset data
  subset_data <- test_data[, c(result$selected_features, "outcome")]
  expect_equal(ncol(subset_data), length(result$selected_features) + 1)

  # Should be able to build a model (if caret is available)
  if (requireNamespace("caret", quietly = TRUE) &&
      requireNamespace("randomForest", quietly = TRUE)) {

    # Simple model to verify features work
    set.seed(123)
    model <- caret::train(
      outcome ~ .,
      data = subset_data,
      method = "rf",
      trControl = caret::trainControl(method = "cv", number = 3),
      tuneLength = 1
    )

    expect_s3_class(model, "train")
  }
})

# ===== Test 13: Parallel Processing =====

test_that("Parallel processing works without errors", {
  skip_if_not_installed("glmnet")
  skip_on_cran()  # Parallel tests can be flaky on CRAN

  # Test with parallel = TRUE
  result_parallel <- OmicSelector_stable_features(
    data = test_data,
    outcome = "outcome",
    method = "stability_selection",
    n_iterations = 20,
    parallel = TRUE,
    cores = 2,
    seed = 123,
    verbose = FALSE
  )

  # Test with parallel = FALSE
  result_sequential <- OmicSelector_stable_features(
    data = test_data,
    outcome = "outcome",
    method = "stability_selection",
    n_iterations = 20,
    parallel = FALSE,
    seed = 123,
    verbose = FALSE
  )

  # Results should be similar (within numerical precision)
  expect_equal(
    sort(result_parallel$selected_features),
    sort(result_sequential$selected_features)
  )
})

# ===== Test 14: Memory and Performance =====

test_that("Function handles moderately large datasets", {
  skip_if_not_installed("glmnet")
  skip_on_cran()  # Performance tests can be slow

  # Create larger dataset
  large_data <- create_test_data(n_samples = 200, n_features = 200, n_informative = 20)

  # Should complete without running out of memory
  expect_error(
    result <- OmicSelector_stable_features(
      data = large_data,
      outcome = "outcome",
      method = "stability_selection",
      n_iterations = 10,  # Keep low for testing
      parallel = FALSE,
      verbose = FALSE
    ),
    NA  # No error expected
  )
})

# ===== Test 15: Documentation and Help =====

test_that("Function has proper documentation", {
  # Check that help exists
  expect_true("OmicSelector_stable_features" %in% ls("package:OmicSelector"))

  # Check that function is exported
  exports <- getNamespaceExports("OmicSelector")
  expect_true("OmicSelector_stable_features" %in% exports)
})
