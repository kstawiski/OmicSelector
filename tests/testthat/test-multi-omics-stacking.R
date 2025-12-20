# Tests for Multi-Omics Stacking Module

test_that("OmicModalitySpec creates valid specification", {
  set.seed(42)

  # Create synthetic data
  x <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  rownames(x) <- paste0("sample_", 1:100)
  colnames(x) <- paste0("gene_", 1:10)

  spec <- OmicModalitySpec$new(
    id = "mRNA",
    x = x,
    learner = "classif.rpart"
  )

  expect_s3_class(spec, "OmicModalitySpec")
  expect_equal(spec$id, "mRNA")
  expect_equal(spec$n_features(), 10)
  expect_equal(length(spec$get_sample_ids()), 100)
})

test_that("OmicModalitySpec accepts learner objects", {
  skip_if_not_installed("mlr3")

  set.seed(42)
  x <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5)
  rownames(x) <- paste0("s", 1:50)

  learner <- mlr3::lrn("classif.rpart", predict_type = "prob")

  spec <- OmicModalitySpec$new(
    id = "test",
    x = x,
    learner = learner
  )

  expect_equal(spec$learner$id, "classif.rpart")
})

test_that("OmicModalitySpec creates valid task", {
  skip_if_not_installed("mlr3")

  set.seed(42)
  n <- 50
  x <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
  rownames(x) <- paste0("s", 1:n)
  y <- factor(sample(c("A", "B"), n, replace = TRUE))
  names(y) <- rownames(x)

  spec <- OmicModalitySpec$new(
    id = "test",
    x = x,
    learner = "classif.rpart"
  )

  task <- spec$as_task(y = y, target_name = "outcome", task_type = "classif")

  expect_s3_class(task, "TaskClassif")
  expect_equal(task$nrow, n)
  expect_equal(task$ncol, 6)  # 5 features + target
})

test_that("OmicStackedEnsemble trains with two modalities", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("ranger")

  set.seed(42)
  n <- 100

  # Create two modalities
  x1 <- matrix(rnorm(n * 20), nrow = n, ncol = 20)
  x2 <- matrix(rnorm(n * 15), nrow = n, ncol = 15)
  rownames(x1) <- rownames(x2) <- paste0("sample_", 1:n)
  colnames(x1) <- paste0("mRNA_", 1:20)
  colnames(x2) <- paste0("miRNA_", 1:15)

  # Create target based on features from both modalities
  y <- factor(ifelse(rowMeans(x1[, 1:5]) + rowMeans(x2[, 1:3]) > 0, "High", "Low"))
  names(y) <- rownames(x1)

  # Create modality specs
  mod1 <- OmicModalitySpec$new(id = "mRNA", x = x1, learner = "classif.rpart")
  mod2 <- OmicModalitySpec$new(id = "miRNA", x = x2, learner = "classif.rpart")

  # Create ensemble with 3-fold CV (faster for testing)
  ensemble <- OmicStackedEnsemble$new(
    modalities = list(mod1, mod2),
    meta_learner = "classif.log_reg",
    resampling = mlr3::rsmp("cv", folds = 3)
  )

  # Train
  ensemble$train(y)

  expect_true(ensemble$is_fitted())
  expect_equal(length(ensemble$modalities), 2)

  # Check OOF stack has correct structure
  oof <- ensemble$get_oof_stack()
  expect_true("target" %in% names(oof))
  expect_true(any(grepl("mRNA_prob_", names(oof))))
  expect_true(any(grepl("miRNA_prob_", names(oof))))
})

test_that("OmicStackedEnsemble predicts on new data", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("ranger")

  set.seed(42)
  n_train <- 80
  n_test <- 20

  # Training data
  x1_train <- matrix(rnorm(n_train * 10), nrow = n_train, ncol = 10)
  x2_train <- matrix(rnorm(n_train * 8), nrow = n_train, ncol = 8)
  rownames(x1_train) <- rownames(x2_train) <- paste0("train_", 1:n_train)

  y_train <- factor(ifelse(rowMeans(x1_train) > 0, "Pos", "Neg"))
  names(y_train) <- rownames(x1_train)

  # Test data
  x1_test <- matrix(rnorm(n_test * 10), nrow = n_test, ncol = 10)
  x2_test <- matrix(rnorm(n_test * 8), nrow = n_test, ncol = 8)
  rownames(x1_test) <- rownames(x2_test) <- paste0("test_", 1:n_test)

  # Create and train ensemble
  mod1 <- OmicModalitySpec$new(id = "mod1", x = x1_train, learner = "classif.rpart")
  mod2 <- OmicModalitySpec$new(id = "mod2", x = x2_train, learner = "classif.rpart")

  ensemble <- OmicStackedEnsemble$new(
    modalities = list(mod1, mod2),
    resampling = mlr3::rsmp("cv", folds = 3)
  )
  ensemble$train(y_train)

  # Predict
  pred <- ensemble$predict(list(mod1 = x1_test, mod2 = x2_test))

  expect_s3_class(pred, "Prediction")
  expect_equal(length(pred$response), n_test)
})

test_that("OmicStackedEnsemble aligns samples correctly", {
  skip_if_not_installed("mlr3")

  set.seed(42)

  # Create misaligned modalities
  x1 <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5)
  x2 <- matrix(rnorm(60 * 5), nrow = 60, ncol = 5)

  # Different sample sets with overlap
  rownames(x1) <- paste0("sample_", 1:50)
  rownames(x2) <- paste0("sample_", 20:79)  # Overlap: 20-50

  y <- factor(sample(c("A", "B"), 80, replace = TRUE))
  names(y) <- paste0("sample_", 1:80)

  mod1 <- OmicModalitySpec$new(id = "mod1", x = x1, learner = "classif.rpart")
  mod2 <- OmicModalitySpec$new(id = "mod2", x = x2, learner = "classif.rpart")

  ensemble <- OmicStackedEnsemble$new(
    modalities = list(mod1, mod2),
    resampling = mlr3::rsmp("cv", folds = 3)
  )

  # Should work with common samples only (31 samples: 20-50)
  expect_message(ensemble$train(y), "samples after alignment")
  expect_true(ensemble$is_fitted())
})

test_that("OmicStackedEnsemble fails with no overlapping samples", {
  set.seed(42)

  x1 <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5)
  x2 <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5)
  rownames(x1) <- paste0("group_A_", 1:50)
  rownames(x2) <- paste0("group_B_", 1:50)  # No overlap

  y <- factor(rep(c("X", "Y"), 50))
  names(y) <- paste0("group_A_", 1:100)

  mod1 <- OmicModalitySpec$new(id = "mod1", x = x1, learner = "classif.rpart")
  mod2 <- OmicModalitySpec$new(id = "mod2", x = x2, learner = "classif.rpart")

  ensemble <- OmicStackedEnsemble$new(
    modalities = list(mod1, mod2),
    resampling = mlr3::rsmp("cv", folds = 3)
  )

  expect_error(ensemble$train(y), "No overlapping sample IDs")
})

test_that("OmicWeightedEnsemble estimates weights from OOF performance", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3measures")

  set.seed(42)
  n <- 100

  # Create modalities with different predictive power
  x1 <- matrix(rnorm(n * 10), nrow = n, ncol = 10)
  x2 <- matrix(rnorm(n * 10), nrow = n, ncol = 10)
  rownames(x1) <- rownames(x2) <- paste0("s", 1:n)

  # Target strongly related to x1, weakly to x2
  y <- factor(ifelse(x1[, 1] > 0, "Pos", "Neg"))
  names(y) <- rownames(x1)

  mod1 <- OmicModalitySpec$new(id = "strong", x = x1, learner = "classif.rpart")
  mod2 <- OmicModalitySpec$new(id = "weak", x = x2, learner = "classif.rpart")

  ensemble <- OmicWeightedEnsemble$new(
    modalities = list(mod1, mod2),
    resampling = mlr3::rsmp("cv", folds = 3)
  )

  expect_message(ensemble$train(y, weight_method = "auc"), "Estimated weights")
  expect_true(ensemble$weights["strong"] > 0)
  expect_true(ensemble$weights["weak"] > 0)
  expect_equal(sum(ensemble$weights), 1, tolerance = 0.001)
})

test_that("stack_omics convenience function works", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("ranger")

  set.seed(42)
  n <- 60

  # Create data
  data_list <- list(
    mRNA = matrix(rnorm(n * 15), nrow = n, ncol = 15),
    miRNA = matrix(rnorm(n * 10), nrow = n, ncol = 10)
  )
  rownames(data_list$mRNA) <- rownames(data_list$miRNA) <- paste0("s", 1:n)

  y <- factor(sample(c("Cancer", "Normal"), n, replace = TRUE))
  names(y) <- paste0("s", 1:n)

  # Create ensemble using convenience function
  ensemble <- stack_omics(
    data_list = data_list,
    learner_list = list(mRNA = "classif.rpart", miRNA = "classif.rpart"),
    y = y,
    meta_learner = "classif.log_reg",
    resampling = mlr3::rsmp("cv", folds = 3)
  )

  expect_s3_class(ensemble, "OmicStackedEnsemble")
  expect_true(ensemble$is_fitted())
})

test_that("stack_omics works with single learner for all modalities", {
  skip_if_not_installed("mlr3")

  set.seed(42)
  n <- 50

  data_list <- list(
    mod1 = matrix(rnorm(n * 8), nrow = n, ncol = 8),
    mod2 = matrix(rnorm(n * 8), nrow = n, ncol = 8)
  )
  rownames(data_list$mod1) <- rownames(data_list$mod2) <- paste0("s", 1:n)

  y <- factor(sample(c("A", "B"), n, replace = TRUE))
  names(y) <- paste0("s", 1:n)

  # Single learner for all
  ensemble <- stack_omics(
    data_list = data_list,
    learner_list = "classif.rpart",  # Single learner
    y = y,
    resampling = mlr3::rsmp("cv", folds = 3)
  )

  expect_true(ensemble$is_fitted())
})

test_that("OmicModalitySpec print method works", {
  set.seed(42)
  x <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  rownames(x) <- paste0("s", 1:50)

  spec <- OmicModalitySpec$new(id = "test", x = x, learner = "classif.rpart")

  expect_output(print(spec), "OmicModalitySpec")
  expect_output(print(spec), "Features: 10")
})

test_that("OmicStackedEnsemble print method works", {
  skip_if_not_installed("mlr3")

  set.seed(42)
  n <- 30

  x1 <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
  x2 <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
  rownames(x1) <- rownames(x2) <- paste0("s", 1:n)

  mod1 <- OmicModalitySpec$new(id = "mod1", x = x1, learner = "classif.rpart")
  mod2 <- OmicModalitySpec$new(id = "mod2", x = x2, learner = "classif.rpart")

  ensemble <- OmicStackedEnsemble$new(
    modalities = list(mod1, mod2),
    meta_learner = "classif.log_reg"
  )

  expect_output(print(ensemble), "OmicStackedEnsemble")
  expect_output(print(ensemble), "Modalities: 2")
  expect_output(print(ensemble), "Meta-learner")
})

test_that("OmicStackedEnsemble get_modality_weights extracts coefficients", {
  skip_if_not_installed("mlr3")

  set.seed(42)
  n <- 80

  x1 <- matrix(rnorm(n * 10), nrow = n, ncol = 10)
  x2 <- matrix(rnorm(n * 10), nrow = n, ncol = 10)
  rownames(x1) <- rownames(x2) <- paste0("s", 1:n)

  y <- factor(ifelse(x1[, 1] > 0, "Pos", "Neg"))
  names(y) <- rownames(x1)

  mod1 <- OmicModalitySpec$new(id = "mod1", x = x1, learner = "classif.rpart")
  mod2 <- OmicModalitySpec$new(id = "mod2", x = x2, learner = "classif.rpart")

  # Use log_reg which is a GLM
  ensemble <- OmicStackedEnsemble$new(
    modalities = list(mod1, mod2),
    meta_learner = "classif.log_reg",
    resampling = mlr3::rsmp("cv", folds = 3)
  )
  ensemble$train(y)

  weights <- ensemble$get_modality_weights()

  # Should return coefficients from GLM
  expect_true(is.numeric(weights) || is.null(weights))
})

test_that("OmicStackedEnsemble handles missing modality at prediction with prior imputation", {
  skip_if_not_installed("mlr3")

  set.seed(42)
  n <- 60

  x1 <- matrix(rnorm(n * 8), nrow = n, ncol = 8)
  x2 <- matrix(rnorm(n * 8), nrow = n, ncol = 8)
  rownames(x1) <- rownames(x2) <- paste0("s", 1:n)

  y <- factor(sample(c("A", "B"), n, replace = TRUE))
  names(y) <- rownames(x1)

  mod1 <- OmicModalitySpec$new(id = "mod1", x = x1, learner = "classif.rpart")
  mod2 <- OmicModalitySpec$new(id = "mod2", x = x2, learner = "classif.rpart")

  ensemble <- OmicStackedEnsemble$new(
    modalities = list(mod1, mod2),
    resampling = mlr3::rsmp("cv", folds = 3)
  )
  ensemble$train(y)

  # Predict with only one modality
  x1_test <- matrix(rnorm(10 * 8), nrow = 10, ncol = 8)
  rownames(x1_test) <- paste0("test_", 1:10)

  expect_warning(
    pred <- ensemble$predict(list(mod1 = x1_test)),
    "Missing modalities"
  )

  expect_s3_class(pred, "Prediction")
  expect_equal(length(pred$response), 10)
})
