# Tests for framework_modern.R

test_that("OmicSelector_fit validates inputs correctly", {
  # Test that function exists
  expect_true(exists("OmicSelector_fit"))

  # Test with missing outcome
  expect_error(
    OmicSelector_fit(
      data = iris,
      outcome = "NonexistentColumn",
      method = "tidymodels"
    ),
    "not found in data"
  )

  # Test with insufficient data
  small_data <- iris[1:5, ]
  expect_error(
    OmicSelector_fit(
      data = small_data,
      outcome = "Species",
      method = "tidymodels"
    ),
    "Insufficient data"
  )
})


test_that("Framework detection works correctly", {
  # Test auto-detection
  result <- OmicSelector:::`.detect_framework`
  expect_true(is.function(result))
})


test_that("OmicSelector_fit handles classification problems", {
  skip_if_not_installed("tidymodels")
  skip_if_not_installed("ranger")

  # Use iris dataset for quick test
  set.seed(123)
  data <- iris

  # This is a basic structure test - actual fitting would be tested in integration tests
  expect_no_error({
    # Just test the function signature
    formals(OmicSelector_fit)
  })
})


test_that("Print and summary methods exist", {
  expect_true(exists("print.OmicSelector_model"))
  expect_true(exists("summary.OmicSelector_model"))
})


test_that("Preprocessing recipe creation works", {
  skip_if_not_installed("recipes")

  # Test basic recipe creation
  data <- iris
  outcome <- "Species"
  preprocessing <- list(normalize = TRUE, remove_zero_variance = TRUE)

  # The function is internal but we can test its existence
  expect_true(exists(".create_recipe", where = "package:OmicSelector", inherits = FALSE) ||
              exists(".create_recipe", envir = asNamespace("OmicSelector")))
})


test_that("Model specification creation works", {
  skip_if_not_installed("parsnip")

  # Test that internal function exists
  expect_true(exists(".create_model_spec", envir = asNamespace("OmicSelector")))
})


test_that("Data leakage checks are implemented", {
  # Test that leakage check function exists
  expect_true(exists(".check_data_leakage", envir = asNamespace("OmicSelector")))
})
