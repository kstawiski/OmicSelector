library(OmicSelector)

test_that("OmicPipeline$new builds a task and features", {
  skip_if_not_installed("mlr3")

  df <- data.frame(
    x1 = rnorm(20),
    x2 = rnorm(20),
    outcome = factor(rep(c("A", "B"), each = 10))
  )

  pipeline <- OmicPipeline$new(
    data = df,
    target = "outcome",
    positive = "A"
  )

  expect_s3_class(pipeline, "OmicPipeline")
  expect_equal(pipeline$get_task()$nrow, 20)
  expect_true(all(c("x1", "x2") %in% pipeline$get_feature_names()))
})
