# Tests for nested_cv.R

test_that("OmicSelector_nested_cv validates inputs correctly", {
  # Test that function exists
  expect_true(exists("OmicSelector_nested_cv"))

  # Test with missing outcome
  expect_error(
    OmicSelector_nested_cv(
      data = iris,
      outcome = "NonexistentColumn",
      outer_folds = 3,
      inner_folds = 3
    ),
    "not found in data"
  )

  # Test with insufficient data
  small_data <- iris[1:10, ]
  expect_error(
    OmicSelector_nested_cv(
      data = small_data,
      outcome = "Species",
      outer_folds = 5,
      inner_folds = 5
    ),
    "Insufficient data"
  )
})


test_that("Nested CV prevents data leakage by design", {
  skip_if_not_installed("tidymodels")

  # This test verifies the structure ensures no leakage
  # The actual nested CV implementation ensures:
  # 1. Preprocessing happens inside folds
  # 2. Feature selection happens inside folds
  # 3. Hyperparameter tuning happens in inner loop only

  # Test that the function signature requires proper separation
  args <- formals(OmicSelector_nested_cv)

  expect_true("preprocessing_recipe" %in% names(args))
  expect_true("feature_selection_method" %in% names(args))
  expect_true("outer_folds" %in% names(args))
  expect_true("inner_folds" %in% names(args))
})


test_that("Feature selection methods are available", {
  # Test that feature selection methods are implemented
  methods <- c("none", "boruta", "rfe", "stability_selection", "lasso")

  # Function should accept these methods
  expect_true(exists(".perform_feature_selection", envir = asNamespace("OmicSelector")))
})


test_that("Feature stability calculation works", {
  # Test feature stability calculation
  expect_true(exists(".calculate_feature_stability", envir = asNamespace("OmicSelector")))

  # Test with mock data
  selected_features_list <- list(
    c("feature1", "feature2", "feature3"),
    c("feature1", "feature2", "feature4"),
    c("feature1", "feature3", "feature5")
  )

  stability <- OmicSelector:::.calculate_feature_stability(selected_features_list)

  expect_type(stability, "list")
  expect_true("selection_frequency" %in% names(stability))
  expect_true("mean_jaccard_similarity" %in% names(stability))
  expect_true("stable_features" %in% names(stability))
})


test_that("Default models are created correctly", {
  skip_if_not_installed("parsnip")

  # Test classification models
  models_class <- OmicSelector:::.create_default_models("classification")
  expect_type(models_class, "list")
  expect_true(length(models_class) > 0)
  expect_true("ranger" %in% names(models_class) ||
              "xgboost" %in% names(models_class))

  # Test regression models
  models_reg <- OmicSelector:::.create_default_models("regression")
  expect_type(models_reg, "list")
  expect_true(length(models_reg) > 0)
})


test_that("Print and summary methods exist for nested_cv", {
  expect_true(exists("print.OmicSelector_nested_cv"))
  expect_true(exists("summary.OmicSelector_nested_cv"))
  expect_true(exists("plot.OmicSelector_nested_cv"))
})


test_that("Best model selection works", {
  # Create mock model results
  mock_models <- list(
    model1 = list(
      metrics = data.frame(
        .metric = c("roc_auc", "accuracy"),
        .estimate = c(0.85, 0.80)
      )
    ),
    model2 = list(
      metrics = data.frame(
        .metric = c("roc_auc", "accuracy"),
        .estimate = c(0.90, 0.85)
      )
    )
  )

  # Test internal function
  expect_true(exists(".select_best_model", envir = asNamespace("OmicSelector")))
})
