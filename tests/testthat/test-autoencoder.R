library(testthat)
library(OmicSelector)

test_that("autoencoder_fit/encode works on small data", {
  skip_if_not_installed("torch")

  x <- matrix(runif(60), nrow = 12, ncol = 5)
  model <- autoencoder_fit(
    x,
    latent_dim = 2,
    hidden_layers = c(6, 4),
    epochs = 2,
    batch_size = 4,
    validation_split = 0.2,
    early_stopping = TRUE
  )

  encoded <- autoencoder_encode(model, x)
  expect_s3_class(encoded, "data.table")
  expect_equal(nrow(encoded), nrow(x))
  expect_equal(ncol(encoded), 2)
})

test_that("autoencoder pipeop returns expected shape", {
  skip_if_not_installed("torch")
  skip_if_not_installed("mlr3pipelines")
  skip_if_not_installed("mlr3learners")

  x <- data.frame(matrix(runif(60), nrow = 12, ncol = 5))
  po <- create_autoencoder_pipeop(latent_dim = 3, epochs = 2, batch_size = 4)
  task <- mlr3::TaskClassif$new("ae", backend = data.frame(x, y = factor(rep(c("a", "b"), 6))), target = "y")
  graph <- po %>>% mlr3pipelines::po("learner", mlr3::lrn("classif.rpart", predict_type = "prob"))
  glrn <- mlr3pipelines::GraphLearner$new(graph)
  glrn$train(task)
  pred <- glrn$predict(task)
  expect_s3_class(pred, "PredictionClassif")
})
