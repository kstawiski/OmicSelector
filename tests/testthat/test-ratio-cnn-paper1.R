library(OmicSelector)

test_that("ratio CNN training and multiseed helper work when torch is available", {
  skip_if_not_installed("torch")
  skip_if_not_installed("coro")
  skip_if_not(torch::torch_is_installed())

  set.seed(1)
  X <- matrix(rnorm(4 * 20), nrow = 20, ncol = 4)
  y <- factor(rep(c("control", "case"), each = 10), levels = c("control", "case"))
  X[y == "case", 1:2] <- X[y == "case", 1:2] + 1.5

  fit <- train_ratio_cnn(
    X_train = X,
    y_train = y,
    X_test = X[1:4, , drop = FALSE],
    epochs = 1,
    batch_size = 4,
    verbose = FALSE,
    seed = 1
  )
  multi <- train_ratio_cnn_multiseed(
    X_train = X,
    y_train = y,
    X_test = X[1:4, , drop = FALSE],
    seeds = c(1, 2, 3),
    epochs = 1,
    batch_size = 4,
    verbose = FALSE
  )

  expect_length(fit$predictions, 4)
  expect_equal(dim(multi$seed_predictions), c(4, 3))
  expect_length(multi$predictions, 4)
})

test_that("ratio CNN learner registers in mlr3", {
  skip_if_not_installed("mlr3")
  expect_true("classif.ratio_cnn" %in% mlr3::mlr_learners$keys())
})

