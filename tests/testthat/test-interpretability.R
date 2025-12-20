# Tests for Interpretability Module

test_that("check_feature_correlations detects correlated pairs", {
  set.seed(42)

  # Create data with correlated features
  n <- 100
  x1 <- rnorm(n)
  x2 <- x1 + rnorm(n, sd = 0.1)  # Highly correlated with x1
  x3 <- rnorm(n)  # Independent

  data <- data.frame(x1 = x1, x2 = x2, x3 = x3)

  result <- check_feature_correlations(data, threshold = 0.7)

  expect_s3_class(result, "CorrelationCheck")
  expect_true(result$warning)
  expect_equal(result$n_pairs, 1)
  expect_true("x1" %in% result$correlated_pairs$feature1 ||
              "x1" %in% result$correlated_pairs$feature2)
})

test_that("check_feature_correlations returns no warning for uncorrelated", {
  set.seed(42)

  # Create uncorrelated data
  data <- data.frame(
    x1 = rnorm(100),
    x2 = rnorm(100),
    x3 = rnorm(100)
  )

  result <- check_feature_correlations(data, threshold = 0.7)

  expect_false(result$warning)
  expect_equal(result$n_pairs, 0)
})

test_that("check_feature_correlations handles single column", {
  data <- data.frame(x = rnorm(100))

  result <- check_feature_correlations(data)

  expect_equal(result$n_pairs, 0)
  expect_false(result$warning)
})

test_that("feature_importance computes valid scores", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")

  set.seed(42)

  # Create simple dataset
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- factor(ifelse(x1 > 0, "A", "B"))  # x1 is predictive, x2 is noise

  data <- data.frame(x1 = x1, x2 = x2, y = y)

  # Train simple model
  task <- mlr3::TaskClassif$new("test", data, target = "y")
  learner <- mlr3::lrn("classif.rpart", predict_type = "prob")
  learner$train(task)

  # Create explainer
  explainer <- create_explainer(
    learner = learner,
    data = data[, c("x1", "x2")],
    y = data$y,
    label = "test_model"
  )

  # Compute importance
  importance <- feature_importance(explainer, n_permutations = 5)

  expect_s3_class(importance, "FeatureImportance")
  expect_true("feature" %in% names(importance))
  expect_true("importance" %in% names(importance))
  expect_equal(nrow(importance), 2)

  # x1 should be more important than x2
  x1_imp <- importance$importance[importance$feature == "x1"]
  x2_imp <- importance$importance[importance$feature == "x2"]
  expect_gt(x1_imp, x2_imp)
})

test_that("partial_dependence computes valid profiles", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")

  set.seed(42)

  # Create dataset
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- factor(ifelse(x1 > 0, "A", "B"))

  data <- data.frame(x1 = x1, x2 = x2, y = y)

  # Train model
  task <- mlr3::TaskClassif$new("test", data, target = "y")
  learner <- mlr3::lrn("classif.rpart", predict_type = "prob")
  learner$train(task)

  # Create explainer
  explainer <- create_explainer(
    learner = learner,
    data = data[, c("x1", "x2")],
    y = data$y
  )

  # Compute partial dependence
  pdp <- partial_dependence(explainer, features = "x1", n_grid = 10)

  expect_s3_class(pdp, "PartialDependence")
  expect_true("value" %in% names(pdp))
  expect_true("prediction" %in% names(pdp))
  expect_equal(nrow(pdp), 10)

  # Predictions should be between 0 and 1
  expect_true(all(pdp$prediction >= 0 & pdp$prediction <= 1))
})

test_that("CorrelationCheck print method works", {
  set.seed(42)

  x1 <- rnorm(50)
  data <- data.frame(x1 = x1, x2 = x1 + rnorm(50, sd = 0.1), x3 = rnorm(50))

  result <- check_feature_correlations(data)

  expect_output(print(result), "Feature Correlation Check")
  expect_output(print(result), "WARNING")
})

test_that("FeatureImportance print method works", {
  skip_if_not_installed("mlr3")

  set.seed(42)

  data <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50),
    y = factor(sample(c("A", "B"), 50, replace = TRUE))
  )

  task <- mlr3::TaskClassif$new("test", data, target = "y")
  learner <- mlr3::lrn("classif.rpart", predict_type = "prob")
  learner$train(task)

  explainer <- create_explainer(
    learner = learner,
    data = data[, c("x1", "x2")],
    y = data$y
  )

  importance <- feature_importance(explainer, n_permutations = 3)

  expect_output(print(importance), "Feature Importance")
})
