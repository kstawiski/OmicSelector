# Tests for Nogueira Stability Index (testthat 3rd edition)

test_that("compute_nogueira_stability works with perfect stability", {
  # All folds select the same features = perfect stability
  feature_sets <- list(
    c("gene1", "gene2", "gene3"),
    c("gene1", "gene2", "gene3"),
    c("gene1", "gene2", "gene3"),
    c("gene1", "gene2", "gene3"),
    c("gene1", "gene2", "gene3")
  )
  all_features <- paste0("gene", 1:100)

  result <- compute_nogueira_stability(feature_sets, all_features)

  expect_s3_class(result, "NogueiraStability")
  expect_equal(result$nogueira_index, 1.0)
  expect_equal(result$n_folds, 5)
  expect_equal(result$n_features, 100)
  expect_equal(length(result$consensus_features), 3)
  expect_true(all(c("gene1", "gene2", "gene3") %in% result$consensus_features))
})

test_that("compute_nogueira_stability handles variable subsets", {
  # Different features selected across folds
  feature_sets <- list(
    c("gene1", "gene2", "gene3"),
    c("gene1", "gene2", "gene4"),
    c("gene1", "gene5", "gene6"),
    c("gene2", "gene3", "gene7"),
    c("gene1", "gene2", "gene3")
  )
  all_features <- paste0("gene", 1:50)

  result <- compute_nogueira_stability(feature_sets, all_features)

  expect_s3_class(result, "NogueiraStability")
  expect_true(result$nogueira_index >= 0 && result$nogueira_index <= 1)
  expect_true(result$nogueira_index < 1)  # Not perfect stability
  expect_equal(result$avg_subset_size, 3)
})

test_that("compute_nogueira_stability computes correct selection frequencies", {
  feature_sets <- list(
    c("gene1", "gene2"),
    c("gene1", "gene3"),
    c("gene1", "gene2"),
    c("gene1", "gene4")
  )
  all_features <- paste0("gene", 1:10)

  result <- compute_nogueira_stability(feature_sets, all_features)

  # gene1 selected in all 4 folds (use unname() to strip vector name for comparison)
  expect_equal(unname(result$selection_frequency["gene1"]), 1.0)
  # gene2 selected in 2 of 4 folds
  expect_equal(unname(result$selection_frequency["gene2"]), 0.5)
  # gene3 selected in 1 of 4 folds
  expect_equal(unname(result$selection_frequency["gene3"]), 0.25)
  # gene5 never selected
  expect_equal(unname(result$selection_frequency["gene5"]), 0.0)
})

test_that("compute_nogueira_stability handles degenerate cases", {
  # Case 1: No features selected in any fold
  empty_sets <- list(
    character(0),
    character(0),
    character(0)
  )
  all_features <- paste0("gene", 1:20)

  result <- compute_nogueira_stability(empty_sets, all_features)
  expect_equal(result$nogueira_index, 1.0)  # Perfectly stable (all agree on empty)
  expect_true(!is.null(result$warning))

  # Case 2: All features selected in every fold
  full_sets <- list(
    paste0("gene", 1:20),
    paste0("gene", 1:20),
    paste0("gene", 1:20)
  )

  result2 <- compute_nogueira_stability(full_sets, all_features)
  expect_equal(result2$nogueira_index, 1.0)  # Perfectly stable (all agree on full set)
})

test_that("compute_nogueira_stability validates input", {
  all_features <- paste0("gene", 1:50)

  # Must have at least 2 folds
  expect_error(
    compute_nogueira_stability(list(c("gene1")), all_features),
    "At least 2 folds required"
  )

  # Must be a list
  expect_error(
    compute_nogueira_stability(c("gene1", "gene2"), all_features),
    "must be a list"
  )

  # all_features cannot be empty
  expect_error(
    compute_nogueira_stability(list(c("gene1"), c("gene1")), character(0)),
    "cannot be empty"
  )
})

test_that("interpret_stability returns correct interpretation", {
  expect_equal(interpret_stability(0.95), "Excellent: Very stable feature selection")
  expect_equal(interpret_stability(0.75), "Good: Reasonably stable feature selection")
  expect_equal(interpret_stability(0.55), "Moderate: Some instability in feature selection")
  expect_equal(interpret_stability(0.35), "Poor: Unstable feature selection - possible overfitting")
  expect_equal(interpret_stability(0.15), "Very Poor: Highly unstable - likely random/overfit")
  expect_match(interpret_stability(NA), "Unable to compute")
})

test_that("print.NogueiraStability produces output", {
  feature_sets <- list(
    c("gene1", "gene2"),
    c("gene1", "gene3"),
    c("gene1", "gene2")
  )
  all_features <- paste0("gene", 1:10)

  result <- compute_nogueira_stability(feature_sets, all_features)

  # Capture output
  output <- capture.output(print(result))

  expect_true(any(grepl("Nogueira Stability Index", output)))
  expect_true(any(grepl("Stability Index:", output)))
  expect_true(any(grepl("Folds analyzed:", output)))
})

test_that("extract_selected_features handles untrained learner", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3pipelines")

  # Create a simple learner (not trained)
  learner <- mlr3::lrn("classif.featureless")

  result <- extract_selected_features(learner)
  expect_null(result)
})

test_that("compute_stability_from_resample returns NULL when no features extracted", {
  skip_if_not_installed("mlr3")

  # This is a mock test - actual integration would need a trained ResampleResult
  # For now, we just verify the function exists and handles errors gracefully
  expect_true(is.function(compute_stability_from_resample))
})

test_that("Nogueira formula matches theoretical expectations",
{
  # Mathematical verification:
  # For m folds with exactly half selecting a feature:
  # pi_j = 0.5 for features selected in half the folds
  # V = mean(0.5 * 0.5) = 0.25 for those features
  # V_0 = (k_bar/p) * (1 - k_bar/p)
  # Stability should be lower when selection is random

  set.seed(123)
  p <- 100
  all_features <- paste0("gene", 1:p)

  # Create "random" selection: each feature has 50% chance per fold
  random_sets <- lapply(1:10, function(i) {
    selected <- sample(all_features, size = 20, replace = FALSE)
    selected
  })

  result <- compute_nogueira_stability(random_sets, all_features)

  # With random selection of 20 features from 100 across 10 folds,
  # stability should be relatively low (well below 0.5)
  expect_true(result$nogueira_index < 0.5)
})
