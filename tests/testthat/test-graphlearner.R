library(OmicSelector)

test_that("create_graph_learner returns a GraphLearner with a filter", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")

  df <- data.frame(
    x1 = rnorm(30),
    x2 = rnorm(30),
    x3 = rnorm(30),
    outcome = factor(rep(c("A", "B"), each = 15))
  )

  pipeline <- OmicPipeline$new(df, target = "outcome", positive = "A")
  learner <- pipeline$create_graph_learner(
    filter = "anova",
    model = "rpart",
    n_features = 2
  )

  expect_s3_class(learner, "GraphLearner")
  has_filter <- any(vapply(learner$graph$pipeops, inherits, logical(1), "PipeOpFilter"))
  expect_true(has_filter)
})
