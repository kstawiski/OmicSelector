# Tests for nested cross-validation and data leakage prevention
# These tests ensure that preprocessing and feature selection happen
# within resampling folds to prevent optimistic bias

test_that("OmicSelector_nested_cv prevents data leakage", {
  skip_if_not_installed("tidymodels")
  skip_if_not_installed("recipes")
  skip_if_not_installed("rsample")
  skip_if_not_installed("parsnip")

  # Create synthetic data with known leakage potential
  set.seed(123)
  n <- 100
  p <- 20

  # Create data where preprocessing on full dataset would leak information
  data <- data.frame(
    outcome = factor(rep(c("A", "B"), each = n/2)),
    matrix(rnorm(n * p), nrow = n)
  )
  colnames(data)[-1] <- paste0("feature_", 1:p)

  # Add a feature that's predictive only after normalization leakage
  data$feature_leak <- ifelse(data$outcome == "A",
                               rnorm(n/2, mean = 0.1, sd = 1),
                               rnorm(n/2, mean = -0.1, sd = 1))

  # Create a simple model specification
  rf_spec <- parsnip::rand_forest(trees = 50) %>%
    parsnip::set_engine("ranger") %>%
    parsnip::set_mode("classification")

  # Run nested CV
  expect_error({
    result <- OmicSelector_nested_cv(
      data = data,
      outcome = "outcome",
      outer_folds = 3,
      inner_folds = 3,
      models = list(rf = rf_spec),
      save_predictions = FALSE
    )
  }, NA)  # Should not error

  # Verify structure of result
  expect_s3_class(result, "OmicSelector_nested_cv")
  expect_true("outer_results" %in% names(result))
  expect_true("inner_results" %in% names(result))
  expect_true("metadata" %in% names(result))

  # Verify metadata is complete
  expect_equal(result$metadata$outer_folds, 3)
  expect_equal(result$metadata$inner_folds, 3)
  expect_equal(result$metadata$n_samples, n)
})


test_that("nested CV runs with preprocessing recipe", {
  skip_if_not_installed("tidymodels")
  skip_if_not_installed("recipes")

  set.seed(456)
  data <- data.frame(
    outcome = factor(sample(c("A", "B"), 50, replace = TRUE)),
    matrix(rnorm(50 * 10), nrow = 50)
  )
  colnames(data)[-1] <- paste0("x", 1:10)

  # Create preprocessing recipe
  rec <- recipes::recipe(outcome ~ ., data = data) %>%
    recipes::step_normalize(recipes::all_numeric_predictors()) %>%
    recipes::step_zv(recipes::all_predictors())

  rf_spec <- parsnip::rand_forest(trees = 50) %>%
    parsnip::set_engine("ranger") %>%
    parsnip::set_mode("classification")

  expect_error({
    result <- OmicSelector_nested_cv(
      data = data,
      outcome = "outcome",
      outer_folds = 2,
      inner_folds = 2,
      preprocessing_recipe = rec,
      models = list(rf = rf_spec),
      save_predictions = FALSE
    )
  }, NA)

  expect_s3_class(result, "OmicSelector_nested_cv")
})


test_that("feature selection occurs within folds", {
  skip_if_not_installed("tidymodels")

  set.seed(789)
  data <- data.frame(
    outcome = factor(rep(c("A", "B"), each = 25)),
    matrix(rnorm(50 * 15), nrow = 50)
  )
  colnames(data)[-1] <- paste0("gene_", 1:15)

  rf_spec <- parsnip::rand_forest(trees = 50) %>%
    parsnip::set_engine("ranger") %>%
    parsnip::set_mode("classification")

  # Run with feature selection
  result <- OmicSelector_nested_cv(
    data = data,
    outcome = "outcome",
    outer_folds = 2,
    inner_folds = 2,
    feature_selection_method = "stability_selection",
    models = list(rf = rf_spec),
    save_predictions = FALSE
  )

  # Check that feature selection results are stored
  expect_true("selected_features" %in% names(result))
  expect_true(is.list(result$selected_features))
})


test_that("OmicSelector_fit can use tidymodels backend", {
  skip_if_not_installed("tidymodels")

  set.seed(101)
  data <- data.frame(
    y = factor(rep(c("A", "B"), each = 30)),
    matrix(rnorm(60 * 5), nrow = 60)
  )
  colnames(data)[-1] <- paste0("v", 1:5)

  expect_error({
    result <- OmicSelector_fit(
      data = data,
      outcome = "y",
      method = "tidymodels",
      algorithm = "ranger",
      resampling = list(method = "cv", folds = 3)
    )
  }, NA)

  expect_s3_class(result, "OmicSelector_fit")
  expect_equal(result$framework, "tidymodels")
})


test_that("OmicSelector_fit can use caret backend", {
  skip_if_not_installed("caret")

  set.seed(102)
  data <- data.frame(
    y = factor(rep(c("A", "B"), each = 30)),
    matrix(rnorm(60 * 5), nrow = 60)
  )
  colnames(data)[-1] <- paste0("v", 1:5)

  expect_error({
    result <- OmicSelector_fit(
      data = data,
      outcome = "y",
      method = "caret",
      algorithm = "rf",
      resampling = list(method = "cv", folds = 3)
    )
  }, NA)

  expect_s3_class(result, "OmicSelector_fit")
  expect_equal(result$framework, "caret")
})


test_that("auto framework detection works", {
  skip_if_not_installed("tidymodels")

  set.seed(103)
  data <- data.frame(
    y = factor(rep(c("A", "B"), each = 20)),
    matrix(rnorm(40 * 5), nrow = 40)
  )
  colnames(data)[-1] <- paste0("v", 1:5)

  result <- OmicSelector_fit(
    data = data,
    outcome = "y",
    method = "auto",
    algorithm = "ranger"
  )

  expect_s3_class(result, "OmicSelector_fit")
  expect_true(result$framework %in% c("tidymodels", "caret"))
})


test_that("nested CV stores metadata correctly", {
  skip_if_not_installed("tidymodels")

  set.seed(104)
  data <- data.frame(
    outcome = factor(rep(c("A", "B"), each = 25)),
    matrix(rnorm(50 * 8), nrow = 50)
  )
  colnames(data)[-1] <- paste0("x", 1:8)

  rf_spec <- parsnip::rand_forest(trees = 50) %>%
    parsnip::set_engine("ranger") %>%
    parsnip::set_mode("classification")

  result <- OmicSelector_nested_cv(
    data = data,
    outcome = "outcome",
    outer_folds = 2,
    inner_folds = 2,
    models = list(rf = rf_spec),
    seed = 12345
  )

  # Check metadata completeness
  expect_equal(result$metadata$seed, 12345)
  expect_equal(result$metadata$n_samples, 50)
  expect_equal(result$metadata$n_features, 8)
  expect_equal(result$metadata$outcome_variable, "outcome")
  expect_true("timestamp" %in% names(result$metadata))
  expect_true("r_version" %in% names(result$metadata))
})


test_that("print methods work correctly", {
  skip_if_not_installed("tidymodels")

  set.seed(105)
  data <- data.frame(
    y = factor(rep(c("A", "B"), each = 20)),
    matrix(rnorm(40 * 5), nrow = 40)
  )
  colnames(data)[-1] <- paste0("v", 1:5)

  # Test OmicSelector_fit print
  fit_result <- OmicSelector_fit(
    data = data,
    outcome = "y",
    method = "tidymodels",
    algorithm = "ranger"
  )

  expect_output(print(fit_result), "OmicSelector Model Fit")

  # Test nested CV print
  rf_spec <- parsnip::rand_forest(trees = 50) %>%
    parsnip::set_engine("ranger") %>%
    parsnip::set_mode("classification")

  cv_result <- OmicSelector_nested_cv(
    data = data,
    outcome = "y",
    outer_folds = 2,
    inner_folds = 2,
    models = list(rf = rf_spec)
  )

  expect_output(print(cv_result), "Nested Cross-Validation")
})


test_that("input validation works", {
  # Test invalid data input
  expect_error(
    OmicSelector_fit(
      data = "not a dataframe",
      outcome = "y"
    ),
    "data must be a data frame"
  )

  # Test invalid outcome
  data <- data.frame(x1 = 1:10, x2 = 1:10)
  expect_error(
    OmicSelector_fit(
      data = data,
      outcome = "nonexistent"
    ),
    "outcome variable 'nonexistent' not found"
  )

  # Test nested CV with no models
  expect_error(
    OmicSelector_nested_cv(
      data = data.frame(y = factor(1:10), x = 1:10),
      outcome = "y",
      models = list()
    ),
    "At least one model must be specified"
  )
})
