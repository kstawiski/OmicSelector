library(OmicSelector)

test_that("nested CV tuning runs and returns stability", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("mlr3tuning")
  skip_if_not_installed("mlr3filters")
  skip_if_not_installed("rpart")

  df <- data.frame(
    x1 = rnorm(40),
    x2 = rnorm(40),
    x3 = rnorm(40),
    x4 = rnorm(40),
    outcome = factor(rep(c("A", "B"), each = 20))
  )

  pipeline <- OmicPipeline$new(df, target = "outcome", positive = "A")
  learner <- pipeline$create_graph_learner(
    filter = "anova",
    model = "rpart",
    n_features = 2
  )

  result <- pipeline$benchmark(
    learners = learner,
    outer_folds = 2,
    inner_folds = 2,
    seed = 123,
    parallel = FALSE,
    threads = 1
  )

  expect_s3_class(result, "NestedCVResult")
  expect_true(!is.null(result$stability))
})
