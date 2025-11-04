# Tests for TRIPOD+AI and PROBAST+AI compliance functions

test_that("TRIPOD+AI checklist is initialized correctly", {
  # Test internal function for checklist initialization
  checklist <- OmicSelector:::.initialize_tripod_checklist()

  expect_s3_class(checklist, "data.frame")
  expect_equal(nrow(checklist), 27)  # 27 TRIPOD+AI items
  expect_true(all(c("section", "item", "description", "status", "page") %in% colnames(checklist)))
  expect_true(all(checklist$status == "Not assessed"))
})


test_that("OmicSelector_tripod_report generates report structure", {
  skip_if_not_installed("tidymodels")

  # Create a mock model result
  set.seed(123)
  data <- data.frame(
    outcome = factor(rep(c("A", "B"), each = 25)),
    matrix(rnorm(50 * 5), nrow = 50)
  )
  colnames(data)[-1] <- paste0("x", 1:5)

  # Fit a simple model
  result <- OmicSelector_fit(
    data = data,
    outcome = "outcome",
    method = "tidymodels",
    algorithm = "ranger"
  )

  # Generate TRIPOD report
  study_info <- list(
    title = "Test Study",
    study_design = "cross-sectional",
    outcome_definition = "Binary outcome"
  )

  tripod_report <- OmicSelector_tripod_report(
    model_result = result,
    study_info = study_info,
    output_format = "json",
    output_file = NULL
  )

  expect_s3_class(tripod_report, "OmicSelector_tripod_report")
  expect_true("checklist" %in% names(tripod_report))
  expect_true("report" %in% names(tripod_report))
  expect_true("completeness_score" %in% names(tripod_report))
  expect_true("metadata" %in% names(tripod_report))

  # Check completeness score structure
  expect_true("n_complete" %in% names(tripod_report$completeness_score))
  expect_true("n_total" %in% names(tripod_report$completeness_score))
  expect_true("percentage" %in% names(tripod_report$completeness_score))
})


test_that("OmicSelector_tripod_report print method works", {
  skip_if_not_installed("tidymodels")

  set.seed(124)
  data <- data.frame(
    outcome = factor(rep(c("A", "B"), each = 20)),
    matrix(rnorm(40 * 5), nrow = 40)
  )
  colnames(data)[-1] <- paste0("x", 1:5)

  result <- OmicSelector_fit(
    data = data,
    outcome = "outcome",
    method = "tidymodels"
  )

  tripod_report <- OmicSelector_tripod_report(
    model_result = result,
    study_info = list(title = "Test"),
    output_format = "json"
  )

  expect_output(print(tripod_report), "TRIPOD\\+AI Compliance Report")
  expect_output(print(tripod_report), "Completeness")
})


test_that("TRIPOD report can be exported to JSON", {
  skip_if_not_installed("tidymodels")
  skip_if_not_installed("jsonlite")

  set.seed(125)
  data <- data.frame(
    outcome = factor(rep(c("A", "B"), each = 20)),
    matrix(rnorm(40 * 5), nrow = 40)
  )
  colnames(data)[-1] <- paste0("x", 1:5)

  result <- OmicSelector_fit(
    data = data,
    outcome = "outcome",
    method = "tidymodels"
  )

  temp_file <- tempfile(fileext = ".json")

  tripod_report <- OmicSelector_tripod_report(
    model_result = result,
    study_info = list(title = "Test"),
    output_format = "json",
    output_file = temp_file
  )

  expect_true(file.exists(temp_file))

  # Clean up
  unlink(temp_file)
})


test_that("TRIPOD report can be exported to HTML", {
  skip_if_not_installed("tidymodels")

  set.seed(126)
  data <- data.frame(
    outcome = factor(rep(c("A", "B"), each = 20)),
    matrix(rnorm(40 * 5), nrow = 40)
  )
  colnames(data)[-1] <- paste0("x", 1:5)

  result <- OmicSelector_fit(
    data = data,
    outcome = "outcome",
    method = "tidymodels"
  )

  temp_file <- tempfile(fileext = ".html")

  tripod_report <- OmicSelector_tripod_report(
    model_result = result,
    study_info = list(title = "Test"),
    output_format = "html",
    output_file = temp_file
  )

  expect_true(file.exists(temp_file))
  html_content <- readLines(temp_file, warn = FALSE)
  expect_true(any(grepl("TRIPOD", html_content)))

  # Clean up
  unlink(temp_file)
})


test_that("PROBAST structure is initialized correctly", {
  probast <- OmicSelector:::.initialize_probast()

  expect_type(probast, "list")
  expect_true(all(c("participants", "predictors", "outcome", "analysis") %in% names(probast)))

  # Check each domain has required elements
  for (domain in names(probast)) {
    expect_true("questions" %in% names(probast[[domain]]))
    expect_true("risk" %in% names(probast[[domain]]))
  }
})


test_that("OmicSelector_probast generates assessment", {
  skip_if_not_installed("tidymodels")

  set.seed(127)
  data <- data.frame(
    outcome = factor(rep(c("A", "B"), each = 25)),
    matrix(rnorm(50 * 5), nrow = 50)
  )
  colnames(data)[-1] <- paste0("x", 1:5)

  rf_spec <- parsnip::rand_forest(trees = 50) %>%
    parsnip::set_engine("ranger") %>%
    parsnip::set_mode("classification")

  # Run nested CV (should result in lower risk due to proper validation)
  result <- OmicSelector_nested_cv(
    data = data,
    outcome = "outcome",
    outer_folds = 2,
    inner_folds = 2,
    models = list(rf = rf_spec)
  )

  probast <- OmicSelector_probast(
    model_result = result,
    output_format = "json"
  )

  expect_s3_class(probast, "OmicSelector_probast")
  expect_true("domain_assessments" %in% names(probast))
  expect_true("overall_risk" %in% names(probast))
  expect_true("recommendations" %in% names(probast))

  # Check overall risk is one of expected values
  expect_true(probast$overall_risk %in% c("Low", "Moderate", "High"))
})


test_that("OmicSelector_probast print method works", {
  skip_if_not_installed("tidymodels")

  set.seed(128)
  data <- data.frame(
    outcome = factor(rep(c("A", "B"), each = 20)),
    matrix(rnorm(40 * 5), nrow = 40)
  )
  colnames(data)[-1] <- paste0("x", 1:5)

  result <- OmicSelector_fit(
    data = data,
    outcome = "outcome",
    method = "tidymodels"
  )

  probast <- OmicSelector_probast(model_result = result)

  expect_output(print(probast), "PROBAST\\+AI")
  expect_output(print(probast), "Overall Risk")
})


test_that("PROBAST can be exported to JSON", {
  skip_if_not_installed("tidymodels")
  skip_if_not_installed("jsonlite")

  set.seed(129)
  data <- data.frame(
    outcome = factor(rep(c("A", "B"), each = 20)),
    matrix(rnorm(40 * 5), nrow = 40)
  )
  colnames(data)[-1] <- paste0("x", 1:5)

  result <- OmicSelector_fit(
    data = data,
    outcome = "outcome",
    method = "tidymodels"
  )

  temp_file <- tempfile(fileext = ".json")

  probast <- OmicSelector_probast(
    model_result = result,
    output_format = "json",
    output_file = temp_file
  )

  expect_true(file.exists(temp_file))

  # Clean up
  unlink(temp_file)
})


test_that("automatic PROBAST assessment detects nested CV", {
  # Test that nested CV is automatically recognized as lower risk

  model_info_nested <- list(
    resampling = list(
      method = "nested_cv",
      outer_folds = 5,
      inner_folds = 5
    )
  )

  assessment <- OmicSelector:::.automatic_probast_assessment(model_info_nested)

  expect_equal(assessment$analysis$risk, "Low")
  expect_true("note" %in% names(assessment$analysis))
})


test_that("OmicSelector_model_card creates card structure", {
  skip_if_not_installed("tidymodels")

  set.seed(130)
  data <- data.frame(
    outcome = factor(rep(c("A", "B"), each = 20)),
    matrix(rnorm(40 * 5), nrow = 40)
  )
  colnames(data)[-1] <- paste0("x", 1:5)

  result <- OmicSelector_fit(
    data = data,
    outcome = "outcome",
    method = "tidymodels"
  )

  card <- OmicSelector_model_card(
    model_result = result,
    model_details = list(name = "Test Model", version = "1.0"),
    intended_use = "For research purposes only",
    limitations = "Small sample size"
  )

  expect_s3_class(card, "OmicSelector_model_card")
  expect_true("model_details" %in% names(card))
  expect_true("intended_use" %in% names(card))
  expect_true("limitations" %in% names(card))
  expect_true("ethical_considerations" %in% names(card))
})


test_that("model card can be exported", {
  skip_if_not_installed("tidymodels")
  skip_if_not_installed("jsonlite")

  set.seed(131)
  data <- data.frame(
    outcome = factor(rep(c("A", "B"), each = 20)),
    matrix(rnorm(40 * 5), nrow = 40)
  )
  colnames(data)[-1] <- paste0("x", 1:5)

  result <- OmicSelector_fit(
    data = data,
    outcome = "outcome",
    method = "tidymodels"
  )

  temp_file <- tempfile(fileext = ".json")

  card <- OmicSelector_model_card(
    model_result = result,
    model_details = list(name = "Test"),
    output_file = temp_file
  )

  expect_true(file.exists(temp_file))

  # Clean up
  unlink(temp_file)
})


test_that("completeness score calculation is correct", {
  checklist <- data.frame(
    status = c(rep("Complete", 15), rep("Not assessed", 12)),
    stringsAsFactors = FALSE
  )

  completeness <- OmicSelector:::.calculate_completeness(checklist)

  expect_equal(completeness$n_complete, 15)
  expect_equal(completeness$n_total, 27)
  expect_equal(completeness$percentage, round(100 * 15/27, 1))
})


test_that("overall PROBAST risk calculation is correct", {
  # All low risk
  assessment_low <- list(
    participants = list(risk = "Low"),
    predictors = list(risk = "Low"),
    outcome = list(risk = "Low"),
    analysis = list(risk = "Low")
  )
  expect_equal(OmicSelector:::.calculate_overall_risk(assessment_low), "Low")

  # One high risk
  assessment_high <- assessment_low
  assessment_high$analysis$risk <- "High"
  expect_equal(OmicSelector:::.calculate_overall_risk(assessment_high), "High")

  # One moderate risk
  assessment_mod <- assessment_low
  assessment_mod$predictors$risk <- "Moderate"
  expect_equal(OmicSelector:::.calculate_overall_risk(assessment_mod), "Moderate")

  # Unknown risk
  assessment_unknown <- assessment_low
  assessment_unknown$outcome$risk <- "Unknown"
  expect_equal(OmicSelector:::.calculate_overall_risk(assessment_unknown), "Moderate")
})
