#' @title Model Export for OmicSelector
#'
#' @description
#' Provides model export functionality for deployment in non-R environments.
#' Supports vetiver for R-based deployment/versioning and ONNX export for
#' cross-language inference.
#'
#' @details
#' ## Export Strategy
#'
#' OmicSelector supports two complementary export paths:
#'
#' 1. **vetiver (R deployment)**: Full model including preprocessing pipeline.
#'    Best for R-based serving (Plumber, Shiny, RSConnect).
#'
#' 2. **ONNX (cross-language)**: Estimator-only export with preprocessing manifest.
#'    Best for Python/Java/etc deployment. Requires implementing preprocessing
#'    in target language.
#'
#' ## Supported Learners for ONNX
#'
#' - `classif.xgboost` / `regr.xgboost`: Well supported via xgboost native export
#' - `classif.glmnet` / `regr.glmnet`: Supported (linear coefficients export)
#' - Other learners: Export preprocessing manifest only; use vetiver for deployment
#'
#' @name model-export
NULL


#' @title Export Model as Vetiver for Deployment
#'
#' @description
#' Wraps an mlr3 Learner or GraphLearner in a vetiver model object for
#' versioning and deployment. Vetiver provides model versioning, API generation,
#' and monitoring capabilities.
#'
#' @param learner A trained mlr3 Learner or GraphLearner
#' @param model_name Name for the vetiver model (required)
#' @param task The mlr3 Task used for training (for metadata extraction)
#' @param description Optional description string
#' @param metadata Optional list of additional metadata
#' @param board Optional pins board for immediate deployment. If NULL, returns
#'        the vetiver model without pinning.
#'
#' @return A vetiver_model object (or invisibly if board is provided)
#'
#' @examples
#' \dontrun{
#' library(vetiver)
#' library(pins)
#'
#' # Train model
#' task <- tsk("iris")
#' learner <- lrn("classif.ranger", predict_type = "prob")
#' learner$train(task)
#'
#' # Export as vetiver
#' v <- export_vetiver(
#'   learner = learner,
#'   model_name = "iris_classifier",
#'   task = task,
#'   description = "Random Forest classifier for iris dataset"
#' )
#'
#' # Pin to a board
#' board <- board_folder("models", versioned = TRUE)
#' vetiver_pin_write(board, v)
#'
#' # Or export and pin in one step
#' export_vetiver(learner, "iris_classifier", task, board = board)
#' }
#'
#' @export
export_vetiver <- function(learner,
                            model_name,
                            task,
                            description = NULL,
                            metadata = NULL,
                            board = NULL) {

  if (!requireNamespace("vetiver", quietly = TRUE)) {
    stop("Package 'vetiver' is required for export_vetiver(). Install with: install.packages('vetiver')")
  }

  if (!requireNamespace("pins", quietly = TRUE)) {
    stop("Package 'pins' is required for export_vetiver(). Install with: install.packages('pins')")
  }

  # Validate inputs
  if (!inherits(learner, "Learner")) {
    stop("learner must be an mlr3 Learner or GraphLearner object")
  }

  if (is.null(learner$model)) {
    stop("learner must be trained before export")
  }

  # Extract feature information
  feature_names <- task$feature_names
  target_name <- task$target_names

  # Create prototype (0-row data.frame with correct types)
  train_data <- task$data()
  x_ptype <- as.data.frame(train_data)[0, feature_names, drop = FALSE]

  # Create custom predict function for mlr3
  predict_mlr3 <- function(model, new_data) {
    # model is the mlr3 learner object
    new_data <- data.table::as.data.table(new_data)

    # Ensure columns are in correct order
    new_data <- new_data[, feature_names, with = FALSE]

    # Use mlr3's predict_newdata
    pred <- model$predict_newdata(new_data)

    # Return standardized output
    if (task$task_type == "classif") {
      if ("prob" %in% pred$predict_types) {
        # Return probabilities with .pred_class
        probs <- as.data.frame(pred$prob)
        names(probs) <- paste0(".pred_", names(probs))
        probs$.pred_class <- pred$response
        probs
      } else {
        data.frame(.pred_class = pred$response)
      }
    } else {
      data.frame(.pred = pred$response)
    }
  }

  # Build metadata
  model_metadata <- list(
    features = feature_names,
    target = target_name,
    task_type = task$task_type,
    n_features = length(feature_names),
    learner_id = learner$id,
    created_at = Sys.time(),
    omic_selector_version = as.character(utils::packageVersion("OmicSelector"))
  )

  if (!is.null(metadata)) {
    model_metadata <- c(model_metadata, metadata)
  }

  # Create vetiver model
  v <- vetiver::vetiver_model(
    model = learner,
    model_name = model_name,
    description = description,
    metadata = model_metadata,
    ptype = x_ptype,
    predict_function = predict_mlr3
  )

  # Pin if board provided
  if (!is.null(board)) {
    vetiver::vetiver_pin_write(board, v)
    message(sprintf("Model '%s' pinned to board", model_name))
    invisible(v)
  } else {
    v
  }
}


#' @title Export Model to ONNX Format
#'
#' @description
#' Exports the estimator component of an mlr3 learner to ONNX format for
#' cross-language deployment. Supports xgboost and glmnet learners.
#'
#' @param learner A trained mlr3 Learner or GraphLearner
#' @param path File path for the ONNX model (with .onnx extension)
#' @param task The mlr3 Task used for training (for metadata)
#' @param export_preprocessing Logical, whether to also export preprocessing manifest
#'
#' @return List with export status and paths to exported files
#'
#' @details
#' ONNX export is currently supported for:
#' - xgboost: via native xgboost save + Python conversion (most reliable)
#' - glmnet: via coefficient extraction (linear models)
#'
#' For xgboost, this function exports the native .xgb format. A separate Python
#' script can convert to ONNX. For glmnet, coefficients are exported as JSON.
#'
#' For other learners, only the preprocessing manifest is exported. Use vetiver
#' for R-based deployment of these models.
#'
#' @examples
#' \dontrun{
#' # Train xgboost model
#' task <- tsk("iris")
#' learner <- lrn("classif.xgboost", nrounds = 50)
#' learner$train(task)
#'
#' # Export to ONNX-compatible format
#' export_onnx(learner, "model.onnx", task)
#' }
#'
#' @export
export_onnx <- function(learner,
                         path,
                         task,
                         export_preprocessing = TRUE) {

  # Validate
  if (!inherits(learner, "Learner")) {
    stop("learner must be an mlr3 Learner or GraphLearner object")
  }

  if (is.null(learner$model)) {
    stop("learner must be trained before export")
  }

  results <- list(
    success = FALSE,
    onnx_path = NULL,
    manifest_path = NULL,
    learner_type = NULL,
    message = NULL
  )

  # Determine base path
  base_path <- tools::file_path_sans_ext(path)

  # Extract the actual model from GraphLearner if needed
  model <- .extract_estimator(learner)
  learner_id <- learner$id

  # Determine learner type and export accordingly
  if (inherits(model, "xgb.Booster")) {
    results$learner_type <- "xgboost"
    results <- .export_xgboost(model, base_path, task, results)

  } else if (inherits(model, "glmnet") || inherits(model, "cv.glmnet")) {
    results$learner_type <- "glmnet"
    results <- .export_glmnet(model, base_path, task, results)

  } else if (inherits(model, "ranger")) {
    results$learner_type <- "ranger"
    results$message <- "ONNX export not directly supported for ranger. Use export_vetiver() for R deployment."
    warning(results$message)

  } else {
    results$learner_type <- class(model)[1]
    results$message <- sprintf(
      "ONNX export not supported for learner type '%s'. Use export_vetiver() for R deployment.",
      results$learner_type
    )
    warning(results$message)
  }

  # Export preprocessing manifest if requested
  if (export_preprocessing && inherits(learner, "GraphLearner")) {
    manifest_path <- paste0(base_path, "_preprocessing.json")
    manifest <- .extract_preprocessing_manifest(learner, task)
    jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE)
    results$manifest_path <- manifest_path
    message(sprintf("Preprocessing manifest exported to: %s", manifest_path))
  }

  # Export feature metadata
  metadata_path <- paste0(base_path, "_metadata.json")
  metadata <- .create_export_metadata(learner, task)
  jsonlite::write_json(metadata, metadata_path, auto_unbox = TRUE, pretty = TRUE)
  results$metadata_path <- metadata_path

  results
}


#' Extract estimator from GraphLearner or Learner
#' @keywords internal
.extract_estimator <- function(learner) {
  if (inherits(learner, "GraphLearner")) {
    # Find the learner PipeOp in the graph
    graph <- learner$graph
    pipeops <- graph$pipeops

    # Find PipeOp that is a LearnerClassif or LearnerRegr
    for (op_id in pipeops$ids()) {
      op <- pipeops$get(op_id)
      if (inherits(op, "PipeOpLearner")) {
        # Return the trained model from the learner
        if (!is.null(op$state) && !is.null(op$state$model)) {
          return(op$state$model)
        }
      }
    }

    # Fallback: return the whole model
    learner$model
  } else {
    learner$model
  }
}


#' Export xgboost model
#' @keywords internal
.export_xgboost <- function(model, base_path, task, results) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    results$message <- "xgboost package required for export"
    return(results)
  }

  # Export native xgboost format
  xgb_path <- paste0(base_path, ".xgb")
  xgboost::xgb.save(model, xgb_path)

  # Also export as JSON (more portable)
  json_path <- paste0(base_path, "_xgb.json")
  xgboost::xgb.save(model, json_path)

  results$success <- TRUE
  results$onnx_path <- xgb_path
  results$message <- sprintf(
    "xgboost model exported to: %s\nTo convert to ONNX, use Python: onnxmltools.convert_xgboost()",
    xgb_path
  )
  message(results$message)

  results
}


#' Export glmnet model (coefficients)
#' @keywords internal
.export_glmnet <- function(model, base_path, task, results) {

  # Extract coefficients
  if (inherits(model, "cv.glmnet")) {
    coefs <- coef(model, s = "lambda.min")
  } else {
    coefs <- coef(model)
  }

  # Convert sparse matrix to list
  coef_list <- list()
  for (i in seq_len(ncol(coefs))) {
    coef_matrix <- coefs[, i, drop = FALSE]
    coef_list[[colnames(coefs)[i]]] <- list(
      names = rownames(coef_matrix),
      values = as.vector(coef_matrix)
    )
  }

  # Export as JSON (can be loaded in Python)
  json_path <- paste0(base_path, "_glmnet_coefs.json")
  export_data <- list(
    model_type = "glmnet",
    task_type = task$task_type,
    coefficients = coef_list,
    lambda = if (inherits(model, "cv.glmnet")) model$lambda.min else model$lambda,
    alpha = model$call$alpha
  )

  jsonlite::write_json(export_data, json_path, auto_unbox = TRUE, pretty = TRUE)

  results$success <- TRUE
  results$onnx_path <- json_path
  results$message <- sprintf(
    "glmnet coefficients exported to: %s\nImplement linear prediction in target language.",
    json_path
  )
  message(results$message)

  results
}


#' Extract preprocessing manifest from GraphLearner
#' @keywords internal
.extract_preprocessing_manifest <- function(learner, task) {
  if (!inherits(learner, "GraphLearner")) {
    return(list(preprocessing_steps = list()))
  }

  graph <- learner$graph
  pipeops <- graph$pipeops

  steps <- list()

  for (op_id in pipeops$ids()) {
    op <- pipeops$get(op_id)

    # Skip learner PipeOps
    if (inherits(op, "PipeOpLearner")) {
      next
    }

    step <- list(
      id = op_id,
      class = class(op)[1]
    )

    # Extract state information based on PipeOp type
    if (!is.null(op$state)) {
      if (inherits(op, "PipeOpScale")) {
        step$type <- "scale"
        step$center <- op$state$center
        step$scale <- op$state$scale
      } else if (grepl("Impute", class(op)[1])) {
        step$type <- "imputation"
        step$model <- op$state$model
      } else if (inherits(op, "PipeOpFilter")) {
        step$type <- "filter"
        step$selected_features <- op$state$selected
      } else {
        step$type <- "unknown"
        step$state_available <- TRUE
      }
    }

    steps[[op_id]] <- step
  }

  list(
    preprocessing_steps = steps,
    feature_names = task$feature_names,
    target_name = task$target_names
  )
}


#' Create export metadata
#' @keywords internal
.create_export_metadata <- function(learner, task) {
  list(
    model_id = learner$id,
    task_type = task$task_type,
    feature_names = task$feature_names,
    n_features = length(task$feature_names),
    target_name = task$target_names,
    target_levels = if (task$task_type == "classif") task$class_names else NULL,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    r_version = R.version.string,
    omic_selector_version = as.character(utils::packageVersion("OmicSelector"))
  )
}


#' @title Create Complete Export Bundle
#'
#' @description
#' Creates a complete export bundle containing:
#' - Vetiver model for R deployment
#' - ONNX-compatible model files (if supported)
#' - Preprocessing manifest
#' - Feature metadata
#'
#' @param learner A trained mlr3 Learner or GraphLearner
#' @param model_name Name for the model bundle
#' @param task The mlr3 Task used for training
#' @param output_dir Directory for export files
#' @param description Optional description
#' @param create_api Logical, whether to generate Plumber API script
#'
#' @return List with paths to all exported files
#'
#' @examples
#' \dontrun{
#' # Create complete export bundle
#' bundle <- export_bundle(
#'   learner = trained_learner,
#'   model_name = "biomarker_classifier",
#'   task = my_task,
#'   output_dir = "model_export",
#'   create_api = TRUE
#' )
#' }
#'
#' @export
export_bundle <- function(learner,
                           model_name,
                           task,
                           output_dir,
                           description = NULL,
                           create_api = FALSE) {

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  results <- list(
    model_name = model_name,
    output_dir = output_dir,
    files = list()
  )

  # 1. Export ONNX-compatible format
  onnx_path <- file.path(output_dir, paste0(model_name, ".onnx"))
  onnx_result <- export_onnx(learner, onnx_path, task, export_preprocessing = TRUE)

  results$files$onnx <- onnx_result$onnx_path
  results$files$preprocessing <- onnx_result$manifest_path
  results$files$metadata <- onnx_result$metadata_path
  results$onnx_support <- onnx_result$success

  # 2. Export vetiver model (if available)
  if (requireNamespace("vetiver", quietly = TRUE) && requireNamespace("pins", quietly = TRUE)) {
    vetiver_dir <- file.path(output_dir, "vetiver")
    if (!dir.exists(vetiver_dir)) {
      dir.create(vetiver_dir)
    }

    board <- pins::board_folder(vetiver_dir, versioned = TRUE)
    v <- export_vetiver(learner, model_name, task, description, board = board)
    results$files$vetiver <- vetiver_dir
    results$vetiver_created <- TRUE
  } else {
    results$vetiver_created <- FALSE
    message("vetiver/pins not available. Skipping vetiver export.")
  }

  # 3. Save R model object as RDS (fallback)
  rds_path <- file.path(output_dir, paste0(model_name, ".rds"))
  saveRDS(
    list(
      learner = learner,
      task_info = list(
        feature_names = task$feature_names,
        target_name = task$target_names,
        task_type = task$task_type
      )
    ),
    rds_path
  )
  results$files$rds <- rds_path

  # 4. Generate Plumber API script
  if (create_api) {
    api_path <- file.path(output_dir, "plumber.R")
    .generate_plumber_api(model_name, task, api_path)
    results$files$api <- api_path
    message(sprintf("Plumber API script created: %s", api_path))
  }

  # 5. Create README
  readme_path <- file.path(output_dir, "README.md")
  .create_export_readme(model_name, results, task, readme_path)
  results$files$readme <- readme_path

  message(sprintf("\nExport bundle created in: %s", output_dir))
  message("Contents:")
  for (name in names(results$files)) {
    if (!is.null(results$files[[name]])) {
      message(sprintf("  - %s: %s", name, basename(results$files[[name]])))
    }
  }

  invisible(results)
}


#' Generate Plumber API script
#' @keywords internal
.generate_plumber_api <- function(model_name, task, api_path) {
  api_script <- sprintf('
# Plumber API for %s
# Generated by OmicSelector

library(plumber)

# Load model
model_data <- readRDS("%s.rds")
learner <- model_data$learner
feature_names <- model_data$task_info$feature_names

#* @apiTitle %s API
#* @apiDescription OmicSelector model prediction API

#* Health check
#* @get /health
function() {
  list(status = "ok", model = "%s")
}

#* Get model metadata
#* @get /metadata
function() {
  list(
    model_name = "%s",
    features = feature_names,
    task_type = model_data$task_info$task_type,
    target = model_data$task_info$target_name
  )
}

#* Make predictions
#* @post /predict
#* @param .req The request object
function(.req) {
  # Parse input
  input <- jsonlite::fromJSON(.req$postBody)

  # Convert to data.table
  new_data <- data.table::as.data.table(input)

  # Ensure correct column order
  new_data <- new_data[, feature_names, with = FALSE]

  # Predict
  pred <- learner$predict_newdata(new_data)

  # Return results
  if ("prob" %%in%% pred$predict_types) {
    list(
      predictions = pred$response,
      probabilities = as.data.frame(pred$prob)
    )
  } else {
    list(predictions = pred$response)
  }
}
', model_name, model_name, model_name, model_name, model_name)

  writeLines(api_script, api_path)
}


#' Create export README
#' @keywords internal
.create_export_readme <- function(model_name, results, task, readme_path) {
  readme <- sprintf('# Model Export: %s

Generated by OmicSelector on %s

## Model Information

- **Model Name**: %s
- **Task Type**: %s
- **Features**: %d
- **Target**: %s

## Files

| File | Description |
|------|-------------|
| %s.rds | R model object (for R deployment) |
| %s_metadata.json | Feature and target metadata |
| %s_preprocessing.json | Preprocessing pipeline manifest |

## R Deployment (Recommended)

```r
# Load model
model_data <- readRDS("%s.rds")
learner <- model_data$learner

# Predict on new data
new_data <- data.frame(...)  # Your new data
pred <- learner$predict_newdata(new_data)
```

## Python Deployment

For Python deployment, use the preprocessing manifest and coefficient files
to implement the prediction pipeline.

## API Deployment

If a Plumber API was generated, start it with:

```r
library(plumber)
pr <- plumb("plumber.R")
pr$run(port = 8000)
```
',
    model_name,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    model_name,
    task$task_type,
    length(task$feature_names),
    paste(task$target_names, collapse = ", "),
    model_name,
    model_name,
    model_name,
    model_name
  )

  writeLines(readme, readme_path)
}


#' @title Load Exported Model Bundle
#'
#' @description
#' Loads a model from an export bundle created by export_bundle().
#'
#' @param bundle_dir Directory containing the export bundle
#' @param model_name Name of the model (defaults to directory name)
#'
#' @return List with learner and task_info
#'
#' @export
load_bundle <- function(bundle_dir, model_name = NULL) {
  if (is.null(model_name)) {
    model_name <- basename(bundle_dir)
  }

  rds_path <- file.path(bundle_dir, paste0(model_name, ".rds"))

  if (!file.exists(rds_path)) {
    stop(sprintf("Model file not found: %s", rds_path))
  }

  model_data <- readRDS(rds_path)
  message(sprintf("Loaded model: %s", model_name))
  message(sprintf("  Task type: %s", model_data$task_info$task_type))
  message(sprintf("  Features: %d", length(model_data$task_info$feature_names)))

  model_data
}
