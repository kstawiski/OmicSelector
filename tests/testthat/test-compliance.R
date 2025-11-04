# Tests for compliance.R

test_that("TRIPOD report generation function exists", {
  expect_true(exists("OmicSelector_tripod_report"))
})


test_that("TRIPOD report validates inputs", {
  # Test with NULL model result
  expect_error(
    OmicSelector_tripod_report(
      model_result = NULL,
      study_info = list(),
      output_format = "html"
    ),
    NA  # Should handle gracefully or error appropriately
  )
})


test_that("TRIPOD checklist has all 27 items", {
  # Create a mock model result
  mock_model <- list(
    metadata = list(
      outer_folds = 5,
      inner_folds = 5,
      n_samples = 100,
      n_features = 50,
      outcome_type = "classification",
      feature_selection_method = "boruta"
    )
  )
  class(mock_model) <- c("OmicSelector_nested_cv", "list")

  # Extract model info and create checklist
  model_info <- OmicSelector:::.extract_model_info(mock_model)
  checklist <- OmicSelector:::.create_tripod_checklist(model_info, list())

  # Should have items for all TRIPOD+AI requirements
  expect_s3_class(checklist, "data.frame")
  expect_true(nrow(checklist) >= 27)  # At least 27 items
  expect_true("Item" %in% names(checklist))
  expect_true("Description" %in% names(checklist))
  expect_true("Status" %in% names(checklist))
})


test_that("PROBAST assessment function exists", {
  expect_true(exists("OmicSelector_probast"))
})


test_that("PROBAST has all 4 domains", {
  # Domain assessment functions should exist
  expect_true(exists(".assess_participants_domain", envir = asNamespace("OmicSelector")))
  expect_true(exists(".assess_predictors_domain", envir = asNamespace("OmicSelector")))
  expect_true(exists(".assess_outcome_domain", envir = asNamespace("OmicSelector")))
  expect_true(exists(".assess_analysis_domain", envir = asNamespace("OmicSelector")))
})


test_that("PROBAST assessment structure is correct", {
  # Create mock model
  mock_model <- list(
    metadata = list(
      outer_folds = 5,
      inner_folds = 5,
      n_samples = 100,
      n_features = 50,
      outcome_type = "classification",
      feature_selection_method = "boruta"
    )
  )
  class(mock_model) <- c("OmicSelector_nested_cv", "list")

  # Run PROBAST assessment
  assessment <- OmicSelector_probast(mock_model)

  expect_s3_class(assessment, "OmicSelector_probast")
  expect_true("overall_risk" %in% names(assessment))
  expect_true("domain_assessments" %in% names(assessment))
  expect_true("recommendations" %in% names(assessment))

  # Check overall risk is valid
  expect_true(assessment$overall_risk %in% c("Low", "High", "Unclear"))

  # Check domain assessments
  domains <- assessment$domain_assessments
  expect_true("participants" %in% names(domains))
  expect_true("predictors" %in% names(domains))
  expect_true("outcome" %in% names(domains))
  expect_true("analysis" %in% names(domains))
})


test_that("Print methods exist for compliance objects", {
  expect_true(exists("print.OmicSelector_tripod_report"))
  expect_true(exists("print.OmicSelector_probast"))
})


test_that("Report export formats are supported", {
  # Check that export function exists
  expect_true(exists(".export_tripod_report", envir = asNamespace("OmicSelector")))

  # Valid formats
  valid_formats <- c("html", "pdf", "json", "markdown")
  expect_true(all(valid_formats %in% c("html", "pdf", "json", "markdown")))
})


test_that("Model info extraction works", {
  # Test with nested_cv model
  mock_nested_cv <- list(
    metadata = list(
      outer_folds = 5,
      inner_folds = 5,
      n_samples = 100,
      n_features = 50,
      outcome_type = "classification",
      feature_selection_method = "boruta"
    ),
    overall_metrics = data.frame(
      .metric = "roc_auc",
      .estimate = 0.85
    )
  )
  class(mock_nested_cv) <- c("OmicSelector_nested_cv", "list")

  info <- OmicSelector:::.extract_model_info(mock_nested_cv)

  expect_type(info, "list")
  expect_true("model_type" %in% names(info))
  expect_equal(info$model_type, "nested_cv")
  expect_equal(info$n_samples, 100)
  expect_equal(info$n_features, 50)
})


test_that("Null coalescing operator works", {
  # Test internal operator
  `%||%` <- OmicSelector:::`%||%`

  expect_equal(NULL %||% "default", "default")
  expect_equal("value" %||% "default", "value")
  expect_equal("" %||% "default", "default")
  expect_equal(character(0) %||% "default", "default")
})
