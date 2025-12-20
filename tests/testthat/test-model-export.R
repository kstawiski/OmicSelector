# Tests for Model Export Module

test_that("export_onnx creates metadata file for any learner", {
  skip_if_not_installed("mlr3")

  set.seed(42)

  # Train a model
  task <- mlr3::tsk("iris")
  learner <- mlr3::lrn("classif.rpart")
  learner$train(task)

  # Export
  tmp_dir <- tempdir()
  onnx_path <- file.path(tmp_dir, "test_model.onnx")

  result <- export_onnx(learner, onnx_path, task, export_preprocessing = FALSE)

  # Should create metadata file
  expect_true(file.exists(result$metadata_path))

  # Read and check metadata
  metadata <- jsonlite::read_json(result$metadata_path)
  expect_equal(metadata$task_type, "classif")
  expect_equal(metadata$n_features, 4)
  expect_true("feature_names" %in% names(metadata))
})

test_that("export_onnx works with xgboost learner", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("xgboost")

  set.seed(42)

  # Train xgboost model
  task <- mlr3::tsk("iris")
  learner <- mlr3::lrn("classif.xgboost", nrounds = 10, verbose = 0)
  learner$train(task)

  # Export
  tmp_dir <- tempdir()
  onnx_path <- file.path(tmp_dir, "xgb_model.onnx")

  result <- export_onnx(learner, onnx_path, task)

  expect_true(result$success)
  expect_equal(result$learner_type, "xgboost")
  expect_true(file.exists(result$onnx_path))
})

test_that("export_onnx works with glmnet learner", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("glmnet")

  set.seed(42)

  # Train glmnet model
  task <- mlr3::tsk("iris")
  learner <- mlr3::lrn("classif.glmnet")
  learner$train(task)

  # Export
  tmp_dir <- tempdir()
  onnx_path <- file.path(tmp_dir, "glmnet_model.onnx")

  result <- export_onnx(learner, onnx_path, task)

  expect_true(result$success)
  expect_equal(result$learner_type, "glmnet")

  # Check coefficients file exists
  expect_true(file.exists(result$onnx_path))

  # Read and check coefficients
  coefs <- jsonlite::read_json(result$onnx_path)
  expect_equal(coefs$model_type, "glmnet")
  expect_true("coefficients" %in% names(coefs))
})

test_that("export_onnx warns for unsupported learners", {
  skip_if_not_installed("mlr3")

  set.seed(42)

  task <- mlr3::tsk("iris")
  learner <- mlr3::lrn("classif.rpart")
  learner$train(task)

  tmp_dir <- tempdir()
  onnx_path <- file.path(tmp_dir, "rpart_model.onnx")

  expect_warning(
    result <- export_onnx(learner, onnx_path, task),
    "not supported"
  )

  expect_false(result$success)
})

test_that(".create_export_metadata generates correct structure", {
  skip_if_not_installed("mlr3")

  task <- mlr3::tsk("iris")
  learner <- mlr3::lrn("classif.rpart")
  learner$train(task)

  metadata <- .create_export_metadata(learner, task)

  expect_true("model_id" %in% names(metadata))
  expect_true("task_type" %in% names(metadata))
  expect_true("feature_names" %in% names(metadata))
  expect_true("target_levels" %in% names(metadata))
  expect_equal(length(metadata$feature_names), 4)
  expect_equal(length(metadata$target_levels), 3)  # iris has 3 classes
})

test_that("export_bundle creates complete bundle", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("jsonlite")

  set.seed(42)

  task <- mlr3::tsk("iris")
  learner <- mlr3::lrn("classif.rpart")
  learner$train(task)

  tmp_dir <- file.path(tempdir(), "test_bundle")
  if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)

  result <- export_bundle(
    learner = learner,
    model_name = "test_model",
    task = task,
    output_dir = tmp_dir,
    create_api = TRUE
  )

  expect_equal(result$model_name, "test_model")
  expect_true(dir.exists(result$output_dir))

  # Check files exist
  expect_true(file.exists(result$files$rds))
  expect_true(file.exists(result$files$metadata))
  expect_true(file.exists(result$files$readme))
  expect_true(file.exists(result$files$api))

  # Cleanup
  unlink(tmp_dir, recursive = TRUE)
})

test_that("load_bundle retrieves saved model", {
  skip_if_not_installed("mlr3")

  set.seed(42)

  task <- mlr3::tsk("iris")
  learner <- mlr3::lrn("classif.rpart")
  learner$train(task)

  tmp_dir <- file.path(tempdir(), "load_test_bundle")
  if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)

  result <- export_bundle(
    learner = learner,
    model_name = "load_test",
    task = task,
    output_dir = tmp_dir
  )

  # Load it back
  loaded <- load_bundle(tmp_dir, "load_test")

  expect_true("learner" %in% names(loaded))
  expect_true("task_info" %in% names(loaded))
  expect_s3_class(loaded$learner, "Learner")
  expect_equal(loaded$task_info$task_type, "classif")

  # Verify predictions work
  new_data <- task$data()[1:5, task$feature_names, with = FALSE]
  pred <- loaded$learner$predict_newdata(as.data.frame(new_data))
  expect_equal(length(pred$response), 5)

  # Cleanup
  unlink(tmp_dir, recursive = TRUE)
})

test_that("export_vetiver creates valid vetiver model", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("vetiver")
  skip_if_not_installed("pins")

  set.seed(42)

  task <- mlr3::tsk("iris")
  learner <- mlr3::lrn("classif.rpart", predict_type = "prob")
  learner$train(task)

  v <- export_vetiver(
    learner = learner,
    model_name = "test_vetiver_model",
    task = task,
    description = "Test model"
  )

  expect_s3_class(v, "vetiver_model")
  expect_equal(v$model_name, "test_vetiver_model")
  expect_true("features" %in% names(v$metadata))
})

test_that("export_vetiver pins to board when provided", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("vetiver")
  skip_if_not_installed("pins")

  set.seed(42)

  task <- mlr3::tsk("iris")
  learner <- mlr3::lrn("classif.rpart", predict_type = "prob")
  learner$train(task)

  tmp_dir <- file.path(tempdir(), "test_pin_board")
  if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)

  board <- pins::board_folder(tmp_dir, versioned = TRUE)

  v <- export_vetiver(
    learner = learner,
    model_name = "pinned_model",
    task = task,
    board = board
  )

  # Check pin exists
  pins_list <- pins::pin_list(board)
  expect_true("pinned_model" %in% pins_list)

  # Cleanup
  unlink(tmp_dir, recursive = TRUE)
})

test_that(".extract_preprocessing_manifest works with GraphLearner", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3pipelines")

  set.seed(42)

  task <- mlr3::tsk("iris")

  # Create pipeline
  graph <- mlr3pipelines::po("scale") %>>%
    mlr3::lrn("classif.rpart")

  glrn <- mlr3pipelines::as_learner(graph)
  glrn$train(task)

  manifest <- .extract_preprocessing_manifest(glrn, task)

  expect_true("preprocessing_steps" %in% names(manifest))
  expect_true("feature_names" %in% names(manifest))
  expect_equal(length(manifest$feature_names), 4)

  # Should have scale step
  step_ids <- names(manifest$preprocessing_steps)
  expect_true(any(grepl("scale", step_ids, ignore.case = TRUE)))
})

test_that("export fails for untrained learner", {
  skip_if_not_installed("mlr3")

  learner <- mlr3::lrn("classif.rpart")
  task <- mlr3::tsk("iris")

  expect_error(
    export_onnx(learner, "model.onnx", task),
    "must be trained"
  )

  expect_error(
    export_vetiver(learner, "test", task),
    "must be trained"
  )
})
