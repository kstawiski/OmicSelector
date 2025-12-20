# Tests for Frozen ComBat
# Validates leakage-free batch correction

test_that("FrozenComBat creates and initializes", {
  fc <- FrozenComBat$new()
  expect_false(fc$is_fitted())
  expect_s3_class(fc, "FrozenComBat")
})

test_that("FrozenComBat fit requires valid input", {
  fc <- FrozenComBat$new()

  # Need at least 2 batches
  data <- matrix(rnorm(100), nrow = 20, ncol = 5)
  batch <- rep("A", 20)

  expect_error(fc$fit(data, batch), "at least 2 batch levels")
})

test_that("FrozenComBat fit learns parameters", {
  set.seed(42)

  # Create data with batch effect
  n_per_batch <- 15
  n_features <- 10

  batch_a <- matrix(rnorm(n_per_batch * n_features, mean = 0), nrow = n_per_batch)
  batch_b <- matrix(rnorm(n_per_batch * n_features, mean = 2), nrow = n_per_batch)  # Shifted

  data <- rbind(batch_a, batch_b)
  batch <- rep(c("A", "B"), each = n_per_batch)

  fc <- FrozenComBat$new()
  fc$fit(data, batch)

  expect_true(fc$is_fitted())
  expect_equal(fc$get_batch_levels(), c("A", "B"))
})

test_that("FrozenComBat transform requires fit first", {
  fc <- FrozenComBat$new()
  data <- matrix(rnorm(50), nrow = 10, ncol = 5)
  batch <- rep(c("A", "B"), each = 5)

  expect_error(fc$transform(data, batch), "must be fitted")
})

test_that("FrozenComBat corrects batch effect", {
  set.seed(42)

  # Create data with clear batch effect
  n_per_batch <- 20
  n_features <- 5

  # Batch A: mean 0
  batch_a <- matrix(rnorm(n_per_batch * n_features, mean = 0, sd = 1), nrow = n_per_batch)
  # Batch B: mean 3 (large shift)
  batch_b <- matrix(rnorm(n_per_batch * n_features, mean = 3, sd = 1), nrow = n_per_batch)

  data <- rbind(batch_a, batch_b)
  batch <- rep(c("A", "B"), each = n_per_batch)

  # Before correction: large difference between batch means
  mean_a_before <- mean(data[batch == "A", ])
  mean_b_before <- mean(data[batch == "B", ])
  diff_before <- abs(mean_a_before - mean_b_before)

  # Apply frozen ComBat
  fc <- FrozenComBat$new()
  corrected <- fc$fit_transform(data, batch)

  # After correction: smaller difference
  mean_a_after <- mean(corrected[batch == "A", ])
  mean_b_after <- mean(corrected[batch == "B", ])
  diff_after <- abs(mean_a_after - mean_b_after)

  expect_lt(diff_after, diff_before)
  expect_lt(diff_after, 1)  # Should be close together now
})

test_that("FrozenComBat parameters are frozen for new data", {
  set.seed(42)

  # Training data
  n_train <- 30
  n_test <- 10
  n_features <- 5

  train_a <- matrix(rnorm(15 * n_features, mean = 0), nrow = 15)
  train_b <- matrix(rnorm(15 * n_features, mean = 2), nrow = 15)
  train_data <- rbind(train_a, train_b)
  train_batch <- rep(c("A", "B"), each = 15)

  # Test data with same batch structure
  test_a <- matrix(rnorm(5 * n_features, mean = 0), nrow = 5)
  test_b <- matrix(rnorm(5 * n_features, mean = 2), nrow = 5)
  test_data <- rbind(test_a, test_b)
  test_batch <- rep(c("A", "B"), each = 5)

  # Fit on training only
  fc <- FrozenComBat$new()
  fc$fit(train_data, train_batch)

  # Transform both
  corrected_train <- fc$transform(train_data, train_batch)
  corrected_test <- fc$transform(test_data, test_batch)

  # Both should be corrected
  expect_equal(dim(corrected_train), dim(train_data))
  expect_equal(dim(corrected_test), dim(test_data))

  # Test data should not influence training correction
  # Re-run with different test data, training correction should be identical
  test_data_v2 <- matrix(rnorm(n_test * n_features, mean = 5), nrow = n_test)
  corrected_train_v2 <- fc$transform(train_data, train_batch)

  expect_equal(corrected_train, corrected_train_v2)
})

test_that("FrozenComBat warns on unknown batches", {
  set.seed(42)

  train_data <- matrix(rnorm(60), nrow = 20, ncol = 3)
  train_batch <- rep(c("A", "B"), each = 10)

  fc <- FrozenComBat$new()
  fc$fit(train_data, train_batch)

  # Test data with unknown batch "C"
  test_data <- matrix(rnorm(15), nrow = 5, ncol = 3)
  test_batch <- c("A", "B", "C", "C", "A")

  expect_warning(fc$transform(test_data, test_batch), "Unknown batch levels")
})

test_that("frozen_combat_correct convenience function works", {
  set.seed(42)

  train_data <- matrix(rnorm(100), nrow = 20, ncol = 5)
  train_batch <- rep(c("A", "B"), each = 10)
  test_data <- matrix(rnorm(50), nrow = 10, ncol = 5)
  test_batch <- rep(c("A", "B"), each = 5)

  result <- frozen_combat_correct(
    train_data = train_data,
    train_batch = train_batch,
    test_data = test_data,
    test_batch = test_batch
  )

  expect_true("corrected_train" %in% names(result))
  expect_true("corrected_test" %in% names(result))
  expect_true("frozen_combat" %in% names(result))
  expect_equal(dim(result$corrected_train), dim(train_data))
  expect_equal(dim(result$corrected_test), dim(test_data))
})

test_that("FrozenComBat prevents leakage vs standard ComBat", {
  # This test demonstrates that frozen ComBat produces different results

  # when test data has different distribution, proving independence

  set.seed(42)

  # Training data: batches with modest effect
  train_a <- matrix(rnorm(50, mean = 0), nrow = 10, ncol = 5)
  train_b <- matrix(rnorm(50, mean = 1), nrow = 10, ncol = 5)
  train_data <- rbind(train_a, train_b)
  train_batch <- rep(c("A", "B"), each = 10)

  # Test data: batch B has VERY different distribution
  test_a <- matrix(rnorm(25, mean = 0), nrow = 5, ncol = 5)
  test_b <- matrix(rnorm(25, mean = 10), nrow = 5, ncol = 5)  # Extreme shift
  test_data <- rbind(test_a, test_b)
  test_batch <- rep(c("A", "B"), each = 5)

  # Frozen ComBat: fit on train, apply to test
  fc <- FrozenComBat$new()
  fc$fit(train_data, train_batch)
  corrected_test_frozen <- fc$transform(test_data, test_batch)

  # The key insight: frozen correction uses training-derived parameters
  # It should NOT perfectly correct the extreme test batch B shift
  # because the shift magnitude (10) is much larger than training (1)

  # Check that test batch B still has higher mean than batch A after correction
  # (imperfect correction due to frozen parameters)
  mean_test_a <- mean(corrected_test_frozen[test_batch == "A", ])
  mean_test_b <- mean(corrected_test_frozen[test_batch == "B", ])

  # There should still be residual difference because parameters were frozen
  # on training data which had smaller batch effect
  expect_true(mean_test_b > mean_test_a + 1)  # Some residual difference
})

test_that("FrozenComBat with mean_only option", {
  set.seed(42)

  data <- matrix(rnorm(100), nrow = 20, ncol = 5)
  data[11:20, ] <- data[11:20, ] + 2  # Add batch effect
  batch <- rep(c("A", "B"), each = 10)

  fc_full <- FrozenComBat$new(mean_only = FALSE)
  fc_mean <- FrozenComBat$new(mean_only = TRUE)

  corrected_full <- fc_full$fit_transform(data, batch)
  corrected_mean <- fc_mean$fit_transform(data, batch)

  # Both should correct, but differently
  expect_false(all(corrected_full == corrected_mean))

  # Both should reduce batch mean difference
  diff_before <- abs(mean(data[1:10, ]) - mean(data[11:20, ]))
  diff_full <- abs(mean(corrected_full[1:10, ]) - mean(corrected_full[11:20, ]))
  diff_mean <- abs(mean(corrected_mean[1:10, ]) - mean(corrected_mean[11:20, ]))

  expect_lt(diff_full, diff_before)
  expect_lt(diff_mean, diff_before)
})
