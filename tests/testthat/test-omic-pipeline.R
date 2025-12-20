# Tests for OmicPipeline R6 Class

context("OmicPipeline R6 Class")

# Skip if mlr3 not available
skip_if_not_installed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    testthat::skip(paste("Package", pkg, "not available"))
  }
}

test_that("OmicPipeline can be created from data frame", {
  # Create simple test data
  set.seed(42)
  n <- 100
  data <- data.frame(
    feature_1 = rnorm(n),
    feature_2 = rnorm(n),
    feature_3 = rnorm(n),
    outcome = factor(sample(c("Case", "Control"), n, replace = TRUE))
  )

  # Skip if R6 not available
  skip_if_not_installed("R6")

  # Create pipeline
  pipeline <- OmicPipeline$new(
    data = data,
    target = "outcome",
    positive = "Case"
  )

  # Check structure
  expect_s3_class(pipeline, "OmicPipeline")
  expect_equal(length(pipeline$get_feature_names()), 3)
  expect_true("feature_1" %in% pipeline$get_feature_names())
})

test_that("OmicPipeline validates input data", {
  # Missing target column
  data <- data.frame(x = 1:10, y = 1:10)

  skip_if_not_installed("R6")

  expect_error(
    OmicPipeline$new(data = data, target = "outcome"),
    "not found"  # Matches both old and new error messages
  )
})

test_that("OmicPipeline creates mlr3 Task", {
  skip_if_not_installed("R6")
  skip_if_not_installed("mlr3")

  set.seed(42)
  n <- 50
  data <- data.frame(
    f1 = rnorm(n),
    f2 = rnorm(n),
    outcome = factor(sample(c("A", "B"), n, replace = TRUE))
  )

  pipeline <- OmicPipeline$new(data = data, target = "outcome", positive = "A")
  task <- pipeline$get_task()

  expect_s3_class(task, "TaskClassif")
  expect_equal(task$target_names, "outcome")
  expect_equal(task$nrow, n)
})

test_that("OmicPipeline factory method works", {
  skip_if_not_installed("R6")

  set.seed(42)
  data <- data.frame(
    x1 = rnorm(30),
    x2 = rnorm(30),
    y = factor(rep(c("A", "B"), 15))
  )

  pipeline <- omic_pipeline(data, target = "y", positive = "A")
  expect_s3_class(pipeline, "OmicPipeline")
})

test_that("OmicPipeline creates GraphLearner with proper structure", {
  skip_if_not_installed("R6")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3pipelines")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("mlr3filters")

  set.seed(42)
  n <- 100
  data <- data.frame(
    matrix(rnorm(n * 20), nrow = n),
    outcome = factor(sample(c("Case", "Control"), n, replace = TRUE))
  )
  names(data)[1:20] <- paste0("feature_", 1:20)

  pipeline <- OmicPipeline$new(data = data, target = "outcome", positive = "Case")

  # Create graph learner
  learner <- pipeline$create_graph_learner(
    filter = "anova",
    model = "ranger",
    n_features = 10,
    impute_method = "median",
    scale = TRUE
  )

  expect_s3_class(learner, "GraphLearner")

  # The graph should contain imputation, scaling, filter, and learner
  graph <- learner$graph
  expect_true(length(graph$ids()) >= 3)  # At minimum: impute, filter, learner
})

test_that("OmicPipeline handles patient grouping", {
  skip_if_not_installed("R6")
  skip_if_not_installed("mlr3")

  # Data with patient IDs
  set.seed(42)
  data <- data.frame(
    patient_id = rep(1:20, each = 3),  # 20 patients, 3 samples each
    x1 = rnorm(60),
    x2 = rnorm(60),
    outcome = factor(rep(c("A", "B"), 30))
  )

  pipeline <- OmicPipeline$new(
    data = data,
    target = "outcome",
    positive = "A",
    patient_id = "patient_id"
  )

  # patient_id should not be in features
  expect_false("patient_id" %in% pipeline$get_feature_names())
  expect_equal(length(pipeline$get_feature_names()), 2)
})
