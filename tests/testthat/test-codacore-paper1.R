library(OmicSelector)

test_that("CoDaCoRe fallback fit predicts probabilities", {
  set.seed(1)
  X <- matrix(rnorm(14 * 20), nrow = 20, ncol = 14)
  y <- factor(rep(c("control", "case"), each = 10), levels = c("control", "case"))
  X[y == "case", 1:4] <- X[y == "case", 1:4] + 0.75

  fit <- suppressWarnings(
    train_codacore(
      X,
      y,
      backend = "fallback",
      class_weight = FALSE,
      verbose = FALSE,
      seed = 1
    )
  )
  pred <- suppressWarnings(predict_codacore(fit, X[1:5, , drop = FALSE]))

  expect_s3_class(fit, "omicselector_codacore")
  expect_length(pred, 5)
  expect_true(all(pred >= 0 & pred <= 1))
})

test_that("CoDaCoRe learner registers in mlr3 and predicts", {
  skip_if_not_installed("mlr3")

  task_df <- data.frame(
    matrix(rnorm(14 * 24), nrow = 24, ncol = 14),
    outcome = factor(rep(c("control", "case"), each = 12), levels = c("control", "case"))
  )
  task_df[task_df$outcome == "case", 1:4] <- task_df[task_df$outcome == "case", 1:4] + 0.75

  task <- mlr3::TaskClassif$new(
    id = "codacore_smoke",
    backend = task_df,
    target = "outcome",
    positive = "case"
  )

  expect_true("classif.codacore" %in% mlr3::mlr_learners$keys())

  learner <- mlr3::lrn(
    "classif.codacore",
    backend = "fallback",
    class_weight = FALSE,
    verbose = FALSE,
    seed = 1
  )
  suppressWarnings(learner$train(task))
  pred <- suppressWarnings(learner$predict(task))

  expect_s3_class(pred, "PredictionClassif")
  expect_equal(nrow(pred$prob), task$nrow)
})
