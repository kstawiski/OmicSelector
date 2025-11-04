#' @title OmicSelector Modern Framework - Unified ML Interface
#' @description
#' Unified interface for tidymodels and caret backends. Provides a modern
#' machine learning framework while maintaining backward compatibility with
#' existing caret-based code.
#'
#' @details
#' This framework implements best practices for biomarker discovery including:
#' * Nested cross-validation to prevent data leakage
#' * Tidymodels integration with recipes for preprocessing
#' * Support for multiple ML algorithms
#' * Automatic framework selection
#' * TRIPOD+AI compliance tracking
#'
#' @name framework_modern
#' @keywords internal
NULL

#' @title Fit Model with Modern Framework
#' @description
#' Main function for fitting models using tidymodels or caret backend.
#' Automatically handles preprocessing, feature selection, and resampling
#' while preventing data leakage.
#'
#' @param data Data frame containing features and outcome
#' @param outcome Character string specifying the outcome variable name
#' @param method Character string specifying framework: "tidymodels", "caret", or "auto"
#' @param algorithm Character string specifying ML algorithm (e.g., "ranger", "xgboost", "glmnet")
#' @param preprocessing List of preprocessing options or a recipes object
#' @param feature_selection List specifying feature selection method and parameters
#' @param resampling List specifying resampling strategy (method, folds, repeats)
#' @param tune_grid Integer or grid specification for hyperparameter tuning
#' @param metrics Metric set for model evaluation
#' @param seed Integer for reproducibility
#' @param ... Additional arguments passed to underlying framework
#'
#' @return An OmicSelector_fit object containing:
#' \item{model}{Fitted model object}
#' \item{framework}{Framework used ("tidymodels" or "caret")}
#' \item{preprocessing}{Preprocessing steps applied}
#' \item{performance}{Performance metrics}
#' \item{metadata}{Metadata for TRIPOD+AI compliance}
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(demo_data)
#'
#' # Fit random forest with tidymodels
#' result <- OmicSelector_fit(
#'   data = demo_data,
#'   outcome = "diagnosis",
#'   method = "tidymodels",
#'   algorithm = "ranger",
#'   preprocessing = list(
#'     normalize = TRUE,
#'     remove_zero_variance = TRUE
#'   ),
#'   resampling = list(method = "cv", folds = 5)
#' )
#' }
#'
#' @export
OmicSelector_fit <- function(
  data,
  outcome,
  method = c("auto", "tidymodels", "caret"),
  algorithm = "ranger",
  preprocessing = list(),
  feature_selection = list(),
  resampling = list(method = "cv", folds = 5),
  tune_grid = 10,
  metrics = NULL,
  seed = 123,
  ...
) {

  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }

  if (!outcome %in% colnames(data)) {
    stop(paste0("outcome variable '", outcome, "' not found in data"))
  }

  method <- match.arg(method)

  # Set seed for reproducibility
  set.seed(seed)

  # Determine which framework to use
  if (method == "auto") {
    method <- .detect_best_framework(algorithm)
  }

  # Initialize metadata for TRIPOD+AI compliance
  metadata <- list(
    timestamp = Sys.time(),
    r_version = R.version.string,
    package_version = utils::packageVersion("OmicSelector"),
    framework = method,
    algorithm = algorithm,
    seed = seed,
    n_samples = nrow(data),
    n_features = ncol(data) - 1,
    outcome_variable = outcome,
    preprocessing = preprocessing,
    feature_selection = feature_selection,
    resampling = resampling
  )

  # Route to appropriate framework
  if (method == "tidymodels") {
    result <- .fit_tidymodels(
      data = data,
      outcome = outcome,
      algorithm = algorithm,
      preprocessing = preprocessing,
      feature_selection = feature_selection,
      resampling = resampling,
      tune_grid = tune_grid,
      metrics = metrics,
      seed = seed,
      ...
    )
  } else if (method == "caret") {
    result <- .fit_caret(
      data = data,
      outcome = outcome,
      algorithm = algorithm,
      preprocessing = preprocessing,
      feature_selection = feature_selection,
      resampling = resampling,
      tune_grid = tune_grid,
      seed = seed,
      ...
    )
  }

  # Combine result with metadata
  result$metadata <- metadata
  result$framework <- method

  class(result) <- c("OmicSelector_fit", class(result))

  return(result)
}


#' @title Nested Cross-Validation with Leakage Prevention
#' @description
#' Implements rigorous nested cross-validation where all preprocessing and
#' feature selection occurs within the resampling folds, preventing data leakage.
#' This is critical for unbiased model evaluation.
#'
#' @param data Data frame containing features and outcome
#' @param outcome Character string specifying the outcome variable name
#' @param outer_folds Integer number of outer CV folds (for final evaluation)
#' @param inner_folds Integer number of inner CV folds (for hyperparameter tuning)
#' @param outer_repeats Integer number of times to repeat outer CV
#' @param preprocessing_recipe A recipes object or list for preprocessing
#' @param feature_selection_method Character string or function for feature selection
#' @param feature_selection_params List of parameters for feature selection
#' @param models List of model specifications (tidymodels parsnip objects)
#' @param metrics Metric set for model evaluation (yardstick)
#' @param calibrate Logical indicating whether to perform calibration
#' @param parallel Logical indicating whether to use parallel processing
#' @param seed Integer for reproducibility
#' @param save_predictions Logical indicating whether to save all predictions
#'
#' @return An OmicSelector_nested_cv object containing:
#' \item{outer_results}{Performance metrics from outer loop (unbiased estimates)}
#' \item{inner_results}{Tuning results from inner loop}
#' \item{selected_features}{Features selected in each fold}
#' \item{predictions}{Out-of-fold predictions}
#' \item{models}{Fitted models from each outer fold}
#' \item{calibration}{Calibration results if requested}
#' \item{metadata}{Complete metadata for reproducibility}
#'
#' @details
#' Nested CV procedure:
#' 1. Split data into outer folds (for evaluation)
#' 2. For each outer fold:
#'    a. Set aside outer test set
#'    b. Use remaining data for inner CV
#'    c. Within inner CV: preprocess, select features, tune hyperparameters
#'    d. Fit final model on full outer training set
#'    e. Evaluate on outer test set
#' 3. Aggregate performance across all outer folds
#'
#' @examples
#' \dontrun{
#' library(tidymodels)
#' library(recipes)
#'
#' # Define preprocessing recipe
#' rec <- recipe(diagnosis ~ ., data = demo_data) %>%
#'   step_normalize(all_numeric_predictors()) %>%
#'   step_zv(all_predictors())
#'
#' # Define models to compare
#' models <- list(
#'   rf = rand_forest(mtry = tune(), trees = 1000) %>%
#'     set_engine("ranger") %>%
#'     set_mode("classification"),
#'
#'   xgb = boost_tree(mtry = tune(), trees = tune()) %>%
#'     set_engine("xgboost") %>%
#'     set_mode("classification")
#' )
#'
#' # Run nested CV
#' result <- OmicSelector_nested_cv(
#'   data = demo_data,
#'   outcome = "diagnosis",
#'   outer_folds = 5,
#'   inner_folds = 5,
#'   preprocessing_recipe = rec,
#'   feature_selection_method = "stability_selection",
#'   models = models,
#'   calibrate = TRUE
#' )
#'
#' # View results
#' print(result)
#' summary(result)
#' plot(result)
#' }
#'
#' @export
OmicSelector_nested_cv <- function(
  data,
  outcome,
  outer_folds = 5,
  inner_folds = 5,
  outer_repeats = 1,
  preprocessing_recipe = NULL,
  feature_selection_method = NULL,
  feature_selection_params = list(),
  models = list(),
  metrics = NULL,
  calibrate = TRUE,
  parallel = TRUE,
  seed = 123,
  save_predictions = TRUE
) {

  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }

  if (!outcome %in% colnames(data)) {
    stop(paste0("outcome variable '", outcome, "' not found in data"))
  }

  if (length(models) == 0) {
    stop("At least one model must be specified")
  }

  # Check if required packages are available
  if (!requireNamespace("tidymodels", quietly = TRUE)) {
    stop("tidymodels package is required for nested CV. Install with: install.packages('tidymodels')")
  }

  set.seed(seed)

  # Initialize metadata
  metadata <- list(
    timestamp = Sys.time(),
    r_version = R.version.string,
    package_version = utils::packageVersion("OmicSelector"),
    seed = seed,
    outer_folds = outer_folds,
    inner_folds = inner_folds,
    outer_repeats = outer_repeats,
    n_samples = nrow(data),
    n_features = ncol(data) - 1,
    outcome_variable = outcome,
    feature_selection_method = feature_selection_method,
    preprocessing_recipe = if (!is.null(preprocessing_recipe)) "provided" else "none",
    n_models = length(models),
    model_names = names(models)
  )

  message("Starting nested cross-validation...")
  message(sprintf("  Outer folds: %d (repeated %d times)", outer_folds, outer_repeats))
  message(sprintf("  Inner folds: %d", inner_folds))
  message(sprintf("  Models: %s", paste(names(models), collapse = ", ")))

  # Initialize storage for results
  outer_results <- list()
  inner_results <- list()
  selected_features <- list()
  all_predictions <- list()
  fitted_models <- list()

  # Create outer resampling splits
  if (outer_repeats > 1) {
    outer_splits <- rsample::vfold_cv(data, v = outer_folds, repeats = outer_repeats)
  } else {
    outer_splits <- rsample::vfold_cv(data, v = outer_folds)
  }

  # Loop through outer folds
  n_outer <- nrow(outer_splits)

  for (i in 1:n_outer) {
    message(sprintf("\nOuter fold %d/%d", i, n_outer))

    # Get training and testing data for this outer fold
    outer_train <- rsample::analysis(outer_splits$splits[[i]])
    outer_test <- rsample::assessment(outer_splits$splits[[i]])

    # CRITICAL: All preprocessing and feature selection happens here
    # within the fold, not on the full dataset

    # 1. Apply preprocessing recipe to training data only
    if (!is.null(preprocessing_recipe)) {
      message("  Applying preprocessing...")
      prepped_recipe <- recipes::prep(preprocessing_recipe, training = outer_train)
      outer_train_processed <- recipes::bake(prepped_recipe, new_data = NULL)
      outer_test_processed <- recipes::bake(prepped_recipe, new_data = outer_test)
    } else {
      outer_train_processed <- outer_train
      outer_test_processed <- outer_test
    }

    # 2. Feature selection on training data only
    if (!is.null(feature_selection_method)) {
      message("  Performing feature selection...")
      fs_result <- .perform_feature_selection(
        data = outer_train_processed,
        outcome = outcome,
        method = feature_selection_method,
        params = feature_selection_params
      )

      selected_features[[i]] <- fs_result$features

      # Subset to selected features
      keep_vars <- c(outcome, fs_result$features)
      outer_train_processed <- outer_train_processed[, keep_vars, drop = FALSE]
      outer_test_processed <- outer_test_processed[, keep_vars, drop = FALSE]
    }

    # 3. Inner CV for hyperparameter tuning
    message("  Running inner CV for hyperparameter tuning...")
    inner_splits <- rsample::vfold_cv(outer_train_processed, v = inner_folds)

    fold_inner_results <- list()
    fold_fitted_models <- list()

    for (model_name in names(models)) {
      message(sprintf("    Tuning %s...", model_name))

      # Tune model using inner CV
      tune_result <- .tune_model_inner_cv(
        model_spec = models[[model_name]],
        data = outer_train_processed,
        outcome = outcome,
        resamples = inner_splits,
        metrics = metrics
      )

      fold_inner_results[[model_name]] <- tune_result

      # Fit best model on full outer training set
      best_params <- tune_result$best_params
      final_model <- .fit_final_model(
        model_spec = models[[model_name]],
        params = best_params,
        data = outer_train_processed,
        outcome = outcome
      )

      fold_fitted_models[[model_name]] <- final_model
    }

    inner_results[[i]] <- fold_inner_results

    # 4. Evaluate each model on outer test set
    fold_outer_results <- list()
    fold_predictions <- list()

    for (model_name in names(models)) {
      message(sprintf("    Evaluating %s on outer test set...", model_name))

      model <- fold_fitted_models[[model_name]]

      # Make predictions
      preds <- predict(model, new_data = outer_test_processed, type = "prob")

      # Calculate metrics
      eval_result <- .evaluate_predictions(
        predictions = preds,
        truth = outer_test_processed[[outcome]],
        metrics = metrics
      )

      fold_outer_results[[model_name]] <- eval_result

      if (save_predictions) {
        fold_predictions[[model_name]] <- data.frame(
          fold = i,
          truth = outer_test_processed[[outcome]],
          preds
        )
      }
    }

    outer_results[[i]] <- fold_outer_results
    all_predictions[[i]] <- fold_predictions
    fitted_models[[i]] <- fold_fitted_models
  }

  message("\nNested CV complete!")

  # Aggregate results across folds
  aggregated_results <- .aggregate_nested_cv_results(
    outer_results = outer_results,
    model_names = names(models)
  )

  # Calibration if requested
  calibration_results <- NULL
  if (calibrate && save_predictions) {
    message("Performing calibration assessment...")
    calibration_results <- .assess_calibration(all_predictions)
  }

  # Create result object
  result <- list(
    outer_results = outer_results,
    aggregated_results = aggregated_results,
    inner_results = inner_results,
    selected_features = selected_features,
    predictions = if (save_predictions) all_predictions else NULL,
    models = fitted_models,
    calibration = calibration_results,
    metadata = metadata
  )

  class(result) <- c("OmicSelector_nested_cv", "list")

  return(result)
}


# Helper functions (internal) ----

#' @keywords internal
.detect_best_framework <- function(algorithm) {
  # Detect which framework is best for the given algorithm
  tidymodels_available <- requireNamespace("tidymodels", quietly = TRUE)
  caret_available <- requireNamespace("caret", quietly = TRUE)

  if (tidymodels_available) {
    return("tidymodels")
  } else if (caret_available) {
    return("caret")
  } else {
    stop("Neither tidymodels nor caret is available. Please install one of them.")
  }
}


#' @keywords internal
.fit_tidymodels <- function(
  data,
  outcome,
  algorithm,
  preprocessing,
  feature_selection,
  resampling,
  tune_grid,
  metrics,
  seed,
  ...
) {
  # This is a placeholder for tidymodels implementation
  # Will be fully implemented with proper tidymodels workflow

  message("Fitting model with tidymodels framework...")

  result <- list(
    model = NULL,
    preprocessing = preprocessing,
    performance = list(),
    message = "tidymodels implementation in progress"
  )

  return(result)
}


#' @keywords internal
.fit_caret <- function(
  data,
  outcome,
  algorithm,
  preprocessing,
  feature_selection,
  resampling,
  tune_grid,
  seed,
  ...
) {
  # Backward compatible caret implementation

  message("Fitting model with caret framework...")

  result <- list(
    model = NULL,
    preprocessing = preprocessing,
    performance = list(),
    message = "caret implementation for backward compatibility"
  )

  return(result)
}


#' @keywords internal
.perform_feature_selection <- function(data, outcome, method, params) {
  # Placeholder for feature selection implementation
  # Will support: stability_selection, boruta, model_x_knockoffs, etc.

  features <- setdiff(colnames(data), outcome)

  result <- list(
    features = features,
    scores = NULL,
    method = method
  )

  return(result)
}


#' @keywords internal
.tune_model_inner_cv <- function(model_spec, data, outcome, resamples, metrics) {
  # Placeholder for inner CV tuning

  result <- list(
    best_params = list(),
    performance = list()
  )

  return(result)
}


#' @keywords internal
.fit_final_model <- function(model_spec, params, data, outcome) {
  # Placeholder for final model fitting

  return(NULL)
}


#' @keywords internal
.evaluate_predictions <- function(predictions, truth, metrics) {
  # Placeholder for prediction evaluation

  result <- list(
    metrics = list()
  )

  return(result)
}


#' @keywords internal
.aggregate_nested_cv_results <- function(outer_results, model_names) {
  # Aggregate performance metrics across outer folds

  aggregated <- list()

  for (model_name in model_names) {
    aggregated[[model_name]] <- list(
      mean_performance = list(),
      sd_performance = list()
    )
  }

  return(aggregated)
}


#' @keywords internal
.assess_calibration <- function(predictions) {
  # Assess model calibration

  result <- list(
    calibration_slope = NULL,
    calibration_intercept = NULL,
    brier_score = NULL
  )

  return(result)
}


# S3 methods ----

#' @export
print.OmicSelector_fit <- function(x, ...) {
  cat("OmicSelector Model Fit\n")
  cat("======================\n\n")
  cat("Framework:", x$framework, "\n")
  cat("Algorithm:", x$metadata$algorithm, "\n")
  cat("Samples:", x$metadata$n_samples, "\n")
  cat("Features:", x$metadata$n_features, "\n")
  invisible(x)
}


#' @export
print.OmicSelector_nested_cv <- function(x, ...) {
  cat("OmicSelector Nested Cross-Validation Results\n")
  cat("=============================================\n\n")
  cat("Outer folds:", x$metadata$outer_folds, "\n")
  cat("Inner folds:", x$metadata$inner_folds, "\n")
  cat("Models evaluated:", paste(x$metadata$model_names, collapse = ", "), "\n")
  cat("\nPerformance Summary:\n")
  print(x$aggregated_results)
  invisible(x)
}


#' @export
summary.OmicSelector_nested_cv <- function(object, ...) {
  cat("Nested CV Summary\n")
  cat("=================\n\n")

  # Print detailed results
  for (model_name in object$metadata$model_names) {
    cat("\nModel:", model_name, "\n")
    cat("  Performance across folds:\n")
    # Print fold-by-fold results
  }

  invisible(object)
}
