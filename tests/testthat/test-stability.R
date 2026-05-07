library(OmicSelector)

test_that("extract_selected_features returns selected feature names", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("mlr3filters")
  skip_if_not_installed("rpart")

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

  learner$train(pipeline$get_task())
  features <- extract_selected_features(learner)

  expect_type(features, "character")
  expect_true(length(features) > 0)
})
