# Tests for BenchmarkService R6 Class

context("BenchmarkService R6 Class")

skip_if_not_installed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    testthat::skip(paste("Package", pkg, "not available"))
  }
}

test_that("BenchmarkService can be created from OmicPipeline", {
  skip_if_not_installed("R6")
  skip_if_not_installed("mlr3")

  set.seed(42)
  n <- 100
  data <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    outcome = factor(sample(c("A", "B"), n, replace = TRUE))
  )

  pipeline <- OmicPipeline$new(data = data, target = "outcome", positive = "A")
  service <- BenchmarkService$new(task = pipeline, outer_folds = 5, inner_folds = 3)

  expect_s3_class(service, "BenchmarkService")
})

test_that("BenchmarkService convenience function works", {
  skip_if_not_installed("R6")
  skip_if_not_installed("mlr3")

  set.seed(42)
  data <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50),
    y = factor(rep(c("A", "B"), 25))
  )

  pipeline <- omic_pipeline(data, target = "y", positive = "A")
  service <- omic_benchmark(pipeline, outer_folds = 3, inner_folds = 2)

  expect_s3_class(service, "BenchmarkService")
})

test_that("BenchmarkService requires learners before running", {
  skip_if_not_installed("R6")
  skip_if_not_installed("mlr3")

  set.seed(42)
  data <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50),
    y = factor(rep(c("A", "B"), 25))
  )

  pipeline <- omic_pipeline(data, target = "y", positive = "A")
  service <- omic_benchmark(pipeline, outer_folds = 3)

  # Should error without learners
  expect_error(service$run(), "No learners added")
})

test_that("BenchmarkService warns about plain learners", {
  skip_if_not_installed("R6")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")

  set.seed(42)
  data <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50),
    y = factor(rep(c("A", "B"), 25))
  )

  pipeline <- omic_pipeline(data, target = "y", positive = "A")
  service <- omic_benchmark(pipeline, outer_folds = 3)

  # Adding a plain learner should warn
  plain_learner <- mlr3::lrn("classif.featureless")
  expect_warning(
    service$add_learner(plain_learner),
    "Plain learner detected"
  )
})

test_that("BenchmarkService with GraphLearner does not warn", {
  skip_if_not_installed("R6")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3pipelines")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("mlr3filters")

  set.seed(42)
  n <- 100
  data <- data.frame(
    matrix(rnorm(n * 10), nrow = n),
    outcome = factor(sample(c("A", "B"), n, replace = TRUE))
  )
  names(data)[1:10] <- paste0("f", 1:10)

  pipeline <- OmicPipeline$new(data = data, target = "outcome", positive = "A")
  graph_learner <- pipeline$create_graph_learner(
    filter = "anova", model = "ranger", n_features = 5
  )

  service <- omic_benchmark(pipeline, outer_folds = 3)

  # GraphLearner should not produce warning (messages are OK for info)
  expect_no_warning(service$add_learner(graph_learner, id = "test_learner"))
})
