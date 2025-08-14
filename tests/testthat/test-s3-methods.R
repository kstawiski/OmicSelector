test_that("OmicSelector S3 object creation works", {
  
  # Create a mock results object
  results <- list(
    stamp = "test_123",
    config = get_config(),
    methods = c(1, 2, 3),
    data_info = list(
      train_samples = 100,
      train_features = 50,
      test_samples = 30,
      validation_samples = 20
    ),
    formulas = list(
      "class ~ feature1 + feature2",
      "class ~ feature3 + feature4 + feature5"
    ),
    method_results = list(
      "1" = list(features = c("feature1", "feature2")),
      "2" = list(features = c("feature3", "feature4", "feature5"))
    ),
    errors = list(),
    warnings = list(),
    started_at = Sys.time(),
    completed_at = Sys.time(),
    runtime = 5.2
  )
  
  class(results) <- c("OmicSelector", "list")
  
  # Test S3 class
  expect_s3_class(results, "OmicSelector")
  
  # Test print method
  expect_output(print(results), "OmicSelector Results")
  expect_output(print(results), "Run ID: test_123")
  expect_output(print(results), "Training samples: 100")
  
  # Test summary method
  expect_output(summary(results), "OmicSelector Results Summary")
  
  # Test formula extraction
  formulas_extracted <- formulas(results)
  expect_length(formulas_extracted, 2)
  expect_true(all(grepl("class ~", formulas_extracted)))
  
  # Test feature extraction
  features_all <- features(results)
  expect_type(features_all, "list")
  expect_equal(names(features_all), c("1", "2"))
  
  features_method1 <- features(results, method_id = 1)
  expect_equal(features_method1, c("feature1", "feature2"))
  
  # Test non-existent method
  expect_warning(features(results, method_id = 99))
})

test_that("OmicSelector plotting works", {
  skip_if_not_installed("ggplot2")
  
  # Create mock results with method data
  results <- list(
    stamp = "test_plot",
    method_results = list(
      "1" = list(features = c("f1", "f2")),
      "2" = list(features = c("f3", "f4", "f5")),
      "3" = list(features = c("f6"))
    ),
    formulas = list("formula1", "formula2", "formula3"),
    errors = list()
  )
  class(results) <- c("OmicSelector", "list")
  
  # Test plot creation (should not error)
  expect_silent(plot(results, type = "overview"))
  expect_silent(plot(results, type = "features"))
  expect_silent(plot(results, type = "methods"))
  
  # Test invalid plot type
  expect_error(plot(results, type = "invalid"))
})

test_that("Empty results object handling", {
  
  # Create empty results object
  empty_results <- list(
    stamp = "empty_test",
    method_results = list(),
    formulas = list(),
    errors = list(),
    data_info = list(
      train_samples = 0,
      train_features = 0,
      test_samples = 0,
      validation_samples = 0
    )
  )
  class(empty_results) <- c("OmicSelector", "list")
  
  # Test that methods handle empty results gracefully
  expect_output(print(empty_results), "Completed: 0")
  
  formulas_empty <- formulas(empty_results)
  expect_length(formulas_empty, 0)
  
  features_empty <- features(empty_results)
  expect_length(features_empty, 0)
})
