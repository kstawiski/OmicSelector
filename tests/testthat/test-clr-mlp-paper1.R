library(OmicSelector)

test_that("CLR transform is invariant to uniform additive shifts", {
  x <- matrix(c(5, 3, 1, 4, 2, 0), nrow = 2, byrow = TRUE)
  clr_a <- clr_transform(x)
  clr_b <- clr_transform(x + 11)

  expect_equal(clr_a, clr_b, tolerance = 1e-12)
  expect_equal(rowMeans(clr_a), c(0, 0), tolerance = 1e-12)
})

test_that("CLR MLP fallback fit predicts probabilities", {
  set.seed(1)
  X <- matrix(rnorm(14 * 20), nrow = 20, ncol = 14)
  y <- factor(rep(c("control", "case"), each = 10), levels = c("control", "case"))
  X[y == "case", 1:3] <- X[y == "case", 1:3] + 0.8

  fit <- suppressWarnings(
    train_clr_mlp(X, y, backend = "glm", class_weight = FALSE, verbose = FALSE, seed = 1)
  )
  pred <- suppressWarnings(predict_clr_mlp(fit, X[1:5, , drop = FALSE]))

  expect_s3_class(fit, "omicselector_clr_mlp")
  expect_length(pred, 5)
  expect_true(all(pred >= 0 & pred <= 1))
})

test_that("CLR learner registers in mlr3 and predicts", {
  skip_if_not_installed("mlr3")

  task_df <- data.frame(
    matrix(rnorm(14 * 24), nrow = 24, ncol = 14),
    outcome = factor(rep(c("control", "case"), each = 12), levels = c("control", "case"))
  )
  task_df[task_df$outcome == "case", 1:3] <- task_df[task_df$outcome == "case", 1:3] + 0.8

  task <- mlr3::TaskClassif$new(
    id = "clr_mlp_smoke",
    backend = task_df,
    target = "outcome",
    positive = "case"
  )

  expect_true("classif.clr_mlp" %in% mlr3::mlr_learners$keys())

  learner <- mlr3::lrn(
    "classif.clr_mlp",
    backend = "glm",
    class_weight = FALSE,
    verbose = FALSE,
    seed = 1
  )
  suppressWarnings(learner$train(task))
  pred <- suppressWarnings(learner$predict(task))

  expect_s3_class(pred, "PredictionClassif")
  expect_equal(nrow(pred$prob), task$nrow)
})
