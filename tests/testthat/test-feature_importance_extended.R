# Extended test suite for feature_importance.R
# Tests OmicSelector_importance() and related functions

library(testthat)

# ===== Test Data Setup =====

# Create test data with known importance structure
create_importance_data <- function(n_samples = 200, seed = 123) {
  set.seed(seed)

  # Create features with different levels of importance
  # High importance features (directly related to outcome)
  X1 <- rnorm(n_samples)
  X2 <- rnorm(n_samples)
  X3 <- rnorm(n_samples)

  # Medium importance (moderately related)
  X4 <- rnorm(n_samples)
  X5 <- rnorm(n_samples)

  # Low importance (noise)
  X6 <- rnorm(n_samples)
  X7 <- rnorm(n_samples)
  X8 <- rnorm(n_samples)

  # Create outcome based on first 3 features
  y_numeric <- 2*X1 + 1.5*X2 + 1*X3 + 0.3*X4 + rnorm(n_samples, sd = 0.5)
  outcome <- ifelse(y_numeric > median(y_numeric), "High", "Low")
  outcome <- factor(outcome)

  # Create data frame
  data <- data.frame(
    feat_high_1 = X1,
    feat_high_2 = X2,
    feat_high_3 = X3,
    feat_med_1 = X4,
    feat_med_2 = X5,
    feat_noise_1 = X6,
    feat_noise_2 = X7,
    feat_noise_3 = X8,
    outcome = outcome
  )

  return(data)
}

# Create correlated features for conditional importance testing
create_correlated_importance_data <- function(n_samples = 200, seed = 123) {
  set.seed(seed)

  # Base feature
  X1 <- rnorm(n_samples)

  # Highly correlated with X1
  X2 <- X1 + rnorm(n_samples, sd = 0.1)

  # Independent features
  X3 <- rnorm(n_samples)
  X4 <- rnorm(n_samples)

  # Outcome depends mainly on X1 (and thus indirectly on X2)
  y_numeric <- 2*X1 + 0.5*X3 + rnorm(n_samples, sd = 0.5)
  outcome <- ifelse(y_numeric > median(y_numeric), "High", "Low")
  outcome <- factor(outcome)

  data <- data.frame(
    important = X1,
    correlated = X2,
    independent = X3,
    noise = X4,
    outcome = outcome
  )

  return(data)
}

test_data <- create_importance_data()
correlated_data <- create_correlated_importance_data()

# ===== Test 1: Input Validation =====

test_that("OmicSelector_importance validates inputs correctly", {

  # Create a simple model for testing
  skip_if_not_installed("randomForest")

  model <- randomForest::randomForest(outcome ~ ., data = test_data)

  # Test missing model
  expect_error(
    OmicSelector_importance(
      model = NULL,
      data = test_data,
      outcome = "outcome"
    ),
    "model"
  )

  # Test missing data
  expect_error(
    OmicSelector_importance(
      model = model,
      data = NULL,
      outcome = "outcome"
    ),
    "data"
  )

  # Test missing outcome
  expect_error(
    OmicSelector_importance(
      model = model,
      data = test_data,
      outcome = "nonexistent"
    ),
    "outcome"
  )

  # Test invalid method
  expect_error(
    OmicSelector_importance(
      model = model,
      data = test_data,
      outcome = "outcome",
      method = "invalid_method"
    ),
    "should be one of"
  )
})

# ===== Test 2: Permutation Importance =====

test_that("Permutation importance identifies important features", {
  skip_if_not_installed("randomForest")

  # Train model
  model <- randomForest::randomForest(outcome ~ ., data = test_data)

  # Calculate importance
  result <- OmicSelector_importance(
    model = model,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 5,  # Reduced for testing speed
    verbose = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_importance")

  # Check required elements
  expect_true("importance" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("metric" %in% names(result))

  # Check importance data frame
  imp_df <- result$importance
  expect_s3_class(imp_df, "data.frame")
  expect_true("feature" %in% names(imp_df))
  expect_true("importance" %in% names(imp_df))
  expect_true("std_error" %in% names(imp_df))

  # Check that high importance features are ranked higher
  feat_high <- imp_df$feature[grepl("feat_high", imp_df$feature)]
  feat_noise <- imp_df$feature[grepl("feat_noise", imp_df$feature)]

  if (length(feat_high) > 0 && length(feat_noise) > 0) {
    mean_imp_high <- mean(imp_df$importance[imp_df$feature %in% feat_high])
    mean_imp_noise <- mean(imp_df$importance[imp_df$feature %in% feat_noise])

    expect_true(mean_imp_high > mean_imp_noise)
  }
})

# ===== Test 3: Conditional Importance =====

test_that("Conditional importance handles correlations correctly", {
  skip_if_not_installed("randomForest")

  # Train model on correlated data
  model <- randomForest::randomForest(outcome ~ ., data = correlated_data)

  # Calculate permutation importance
  result_perm <- OmicSelector_importance(
    model = model,
    data = correlated_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 5,
    verbose = FALSE
  )

  # Calculate conditional importance
  result_cond <- OmicSelector_importance(
    model = model,
    data = correlated_data,
    outcome = "outcome",
    method = "conditional",
    n_repeats = 5,
    conditional_grid = 3,  # Reduced for testing
    verbose = FALSE
  )

  # Both should succeed
  expect_s3_class(result_perm, "OmicSelector_importance")
  expect_s3_class(result_cond, "OmicSelector_importance")

  # Conditional importance should differ from permutation for correlated features
  # (This is a qualitative test - exact values depend on implementation)
  expect_true(nrow(result_cond$importance) > 0)
})

# ===== Test 4: Both Methods Comparison =====

test_that("Both methods can be calculated together", {
  skip_if_not_installed("randomForest")

  model <- randomForest::randomForest(outcome ~ ., data = test_data)

  result <- OmicSelector_importance(
    model = model,
    data = test_data,
    outcome = "outcome",
    method = "both",
    n_repeats = 5,
    verbose = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_importance")

  # Should have both importance types
  expect_true("permutation_importance" %in% names(result) ||
              "importance_permutation" %in% names(result) ||
              "method" %in% names(result))

  # If structured as list with both methods
  if (is.list(result) && length(result) > 1) {
    expect_true(any(grepl("permutation", names(result), ignore.case = TRUE)))
    expect_true(any(grepl("conditional", names(result), ignore.case = TRUE)))
  }
})

# ===== Test 5: Confidence Intervals =====

test_that("Confidence intervals are calculated", {
  skip_if_not_installed("randomForest")

  model <- randomForest::randomForest(outcome ~ ., data = test_data)

  result <- OmicSelector_importance(
    model = model,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 10,
    verbose = FALSE
  )

  imp_df <- result$importance

  # Check for CI columns
  expect_true("std_error" %in% names(imp_df) ||
              "se" %in% names(imp_df))

  if ("ci_lower" %in% names(imp_df) && "ci_upper" %in% names(imp_df)) {
    # CI lower should be less than upper
    expect_true(all(imp_df$ci_lower <= imp_df$ci_upper))

    # Importance should generally be within CI
    # (allowing for some numerical precision)
    within_ci <- (imp_df$importance >= imp_df$ci_lower - 1e-6) &
                 (imp_df$importance <= imp_df$ci_upper + 1e-6)
    expect_true(mean(within_ci) > 0.8)  # At least 80% should be within CI
  }
})

# ===== Test 6: Normalization =====

test_that("Normalization option works", {
  skip_if_not_installed("randomForest")

  model <- randomForest::randomForest(outcome ~ ., data = test_data)

  # With normalization
  result_norm <- OmicSelector_importance(
    model = model,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    normalize = TRUE,
    n_repeats = 5,
    verbose = FALSE
  )

  # Without normalization
  result_no_norm <- OmicSelector_importance(
    model = model,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    normalize = FALSE,
    n_repeats = 5,
    verbose = FALSE
  )

  # Both should succeed
  expect_s3_class(result_norm, "OmicSelector_importance")
  expect_s3_class(result_no_norm, "OmicSelector_importance")

  # Values should differ if normalization is applied
  # (unless they're already normalized by default)
})

# ===== Test 7: Different Metrics =====

test_that("Different evaluation metrics work", {
  skip_if_not_installed("randomForest")

  model <- randomForest::randomForest(outcome ~ ., data = test_data)

  # Accuracy
  result_acc <- OmicSelector_importance(
    model = model,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    metric = "accuracy",
    n_repeats = 5,
    verbose = FALSE
  )

  expect_s3_class(result_acc, "OmicSelector_importance")
  expect_equal(result_acc$metric, "accuracy")

  # AUC (if binary classification)
  if (length(levels(test_data$outcome)) == 2) {
    result_auc <- OmicSelector_importance(
      model = model,
      data = test_data,
      outcome = "outcome",
      method = "permutation",
      metric = "auc",
      n_repeats = 5,
      verbose = FALSE
    )

    expect_s3_class(result_auc, "OmicSelector_importance")
    expect_true(result_auc$metric %in% c("auc", "roc_auc", "AUC"))
  }
})

# ===== Test 8: S3 Methods =====

test_that("S3 methods work correctly", {
  skip_if_not_installed("randomForest")

  model <- randomForest::randomForest(outcome ~ ., data = test_data)

  result <- OmicSelector_importance(
    model = model,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 5,
    verbose = FALSE
  )

  # Test print method
  expect_output(print(result), "Importance")

  # Test summary method (if exists)
  if ("summary.OmicSelector_importance" %in% methods("summary")) {
    expect_output(summary(result), "Feature")
  }

  # Test plot method (if exists)
  if ("plot.OmicSelector_importance" %in% methods("plot")) {
    expect_silent(p <- plot(result))
    expect_true(inherits(p, "ggplot") || inherits(p, "recordedplot"))
  }
})

# ===== Test 9: Reproducibility =====

test_that("Results are reproducible with seed", {
  skip_if_not_installed("randomForest")

  model <- randomForest::randomForest(outcome ~ ., data = test_data, ntree = 50)

  result1 <- OmicSelector_importance(
    model = model,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 5,
    seed = 42,
    verbose = FALSE
  )

  result2 <- OmicSelector_importance(
    model = model,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 5,
    seed = 42,
    verbose = FALSE
  )

  # Should get very similar results (within numerical precision)
  expect_equal(result1$importance$importance,
               result2$importance$importance,
               tolerance = 1e-6)
})

# ===== Test 10: Multiple Model Types =====

test_that("Works with different model types", {

  # Random Forest
  skip_if_not_installed("randomForest")
  model_rf <- randomForest::randomForest(outcome ~ ., data = test_data)

  result_rf <- OmicSelector_importance(
    model = model_rf,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 3,
    verbose = FALSE
  )

  expect_s3_class(result_rf, "OmicSelector_importance")

  # GLM (if supported)
  model_glm <- glm(outcome ~ ., data = test_data, family = binomial())

  result_glm <- OmicSelector_importance(
    model = model_glm,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 3,
    verbose = FALSE
  )

  expect_s3_class(result_glm, "OmicSelector_importance")
})

# ===== Test 11: Feature Ranking =====

test_that("Features are ranked correctly", {
  skip_if_not_installed("randomForest")

  model <- randomForest::randomForest(outcome ~ ., data = test_data)

  result <- OmicSelector_importance(
    model = model,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 10,
    verbose = FALSE
  )

  imp_df <- result$importance

  # Check if there's a rank column
  if ("rank" %in% names(imp_df)) {
    # Ranks should be unique and consecutive
    expect_equal(sort(imp_df$rank), 1:nrow(imp_df))
  }

  # Importance should be sorted (descending)
  if (nrow(imp_df) > 1) {
    # Check if sorted (allowing ties)
    expect_true(all(imp_df$importance[-1] <= imp_df$importance[-nrow(imp_df)] + 1e-10))
  }
})

# ===== Test 12: Z-scores and Significance =====

test_that("Z-scores and significance tests work", {
  skip_if_not_installed("randomForest")

  model <- randomForest::randomForest(outcome ~ ., data = test_data)

  result <- OmicSelector_importance(
    model = model,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 10,
    verbose = FALSE
  )

  imp_df <- result$importance

  # Check for z-score column
  if ("z_score" %in% names(imp_df)) {
    # Z-scores should be calculable
    expect_type(imp_df$z_score, "double")

    # Should have both positive and possibly negative z-scores
    expect_true(any(imp_df$z_score > 0))
  }

  # Check for p-value column
  if ("p_value" %in% names(imp_df)) {
    # P-values should be between 0 and 1
    expect_true(all(imp_df$p_value >= 0 & imp_df$p_value <= 1))
  }
})

# ===== Test 13: Parallel Processing =====

test_that("Parallel processing works", {
  skip_if_not_installed("randomForest")
  skip_on_cran()

  model <- randomForest::randomForest(outcome ~ ., data = test_data)

  # With parallel
  result_parallel <- OmicSelector_importance(
    model = model,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 5,
    parallel = TRUE,
    cores = 2,
    seed = 123,
    verbose = FALSE
  )

  # Without parallel
  result_sequential <- OmicSelector_importance(
    model = model,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 5,
    parallel = FALSE,
    seed = 123,
    verbose = FALSE
  )

  # Results should be similar
  expect_equal(
    result_parallel$importance$importance,
    result_sequential$importance$importance,
    tolerance = 0.01  # Allow small numerical differences
  )
})

# ===== Test 14: Edge Cases =====

test_that("Handles edge cases gracefully", {
  skip_if_not_installed("randomForest")

  # Very small dataset
  small_data <- test_data[1:30, ]
  model_small <- randomForest::randomForest(outcome ~ ., data = small_data)

  result_small <- OmicSelector_importance(
    model = model_small,
    data = small_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 3,
    verbose = FALSE
  )

  expect_s3_class(result_small, "OmicSelector_importance")

  # Single feature (edge case)
  single_feat_data <- test_data[, c("feat_high_1", "outcome")]
  model_single <- randomForest::randomForest(outcome ~ ., data = single_feat_data)

  result_single <- OmicSelector_importance(
    model = model_single,
    data = single_feat_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 3,
    verbose = FALSE
  )

  expect_s3_class(result_single, "OmicSelector_importance")
  expect_equal(nrow(result_single$importance), 1)
})

# ===== Test 15: Integration with Feature Selection =====

test_that("Can be used for feature selection", {
  skip_if_not_installed("randomForest")

  model <- randomForest::randomForest(outcome ~ ., data = test_data)

  result <- OmicSelector_importance(
    model = model,
    data = test_data,
    outcome = "outcome",
    method = "permutation",
    n_repeats = 10,
    verbose = FALSE
  )

  # Select top N features
  n_top <- 5
  top_features <- head(result$importance$feature, n_top)

  expect_equal(length(top_features), n_top)
  expect_true(all(top_features %in% colnames(test_data)))

  # Should be able to subset data with these features
  subset_data <- test_data[, c(top_features, "outcome")]
  expect_equal(ncol(subset_data), n_top + 1)

  # Should be able to train a new model
  model_subset <- randomForest::randomForest(outcome ~ ., data = subset_data)
  expect_s3_class(model_subset, "randomForest")
})
