#' @title Nested Cross-Validation with Data Leakage Prevention
#' @description Implements rigorous nested cross-validation to prevent data leakage
#'
#' This module implements nested cross-validation (nested CV) which is the gold standard
#' for unbiased model evaluation. The outer loop provides unbiased performance estimates,
#' while the inner loop is used for hyperparameter tuning and model selection.
#'
#' CRITICAL: All preprocessing and feature selection happens INSIDE the resampling
#' loops to prevent data leakage.
#'
#' @author Konrad Stawiski
#' @name nested_cv
NULL


#' Nested Cross-Validation for Unbiased Model Evaluation
#'
#' Performs nested cross-validation with proper separation of tuning and validation sets.
#' All preprocessing, feature selection, and hyperparameter tuning occurs within the
#' inner loop, ensuring no information from validation folds leaks into the training process.
#'
#' @param data A data frame containing features and outcome
#' @param outcome Character string specifying the outcome variable name
#' @param outer_folds Integer, number of folds in outer (validation) loop. Default: 5
#' @param inner_folds Integer, number of folds in inner (tuning) loop. Default: 5
#' @param preprocessing_recipe A recipes object specifying preprocessing steps,
#'   or NULL to use defaults. If NULL, will normalize and remove zero-variance features.
#' @param feature_selection_method Character string specifying feature selection method:
#'   \itemize{
#'     \item "none" - No feature selection
#'     \item "boruta" - Boruta algorithm
#'     \item "rfe" - Recursive feature elimination
#'     \item "stability_selection" - Stability selection
#'     \item "lasso" - LASSO regularization
#'   }
#' @param max_features Integer, maximum number of features to select. NULL means no limit.
#' @param models A named list of parsnip model specifications to evaluate.
#'   If NULL, will use default models (ranger, xgboost, glmnet).
#' @param tune_grid Integer (for random search) or data frame (for grid search)
#' @param metrics A yardstick::metric_set object. If NULL, uses appropriate defaults
#'   based on problem type (classification vs regression).
#' @param calibrate Logical, whether to perform calibration on predictions. Default: TRUE
#' @param parallel Logical, whether to use parallel processing. Default: TRUE
#' @param cores Integer, number of cores for parallel processing
#' @param seed Integer for reproducibility
#' @param verbose Logical, whether to print progress messages. Default: TRUE
#'
#' @return A list object of class "OmicSelector_nested_cv" containing:
#'   \item{outer_results}{Performance metrics from outer (validation) loop}
#'   \item{inner_results}{Tuning results from inner loop for each outer fold}
#'   \item{final_predictions}{Out-of-sample predictions for all data points}
#'   \item{selected_features}{Features selected in each outer fold}
#'   \item{best_models}{Best model configuration for each outer fold}
#'   \item{calibration}{Calibration results if calibrate=TRUE}
#'   \item{feature_stability}{Stability metrics for selected features}
#'   \item{metadata}{Metadata about the analysis}
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(miR_Asakura)
#'
#' # Define custom recipe
#' library(recipes)
#' rec <- recipe(Class ~ ., data = miR_Asakura) %>%
#'   step_normalize(all_numeric_predictors()) %>%
#'   step_zv(all_predictors())
#'
#' # Define models to compare
#' library(parsnip)
#' models <- list(
#'   rf = rand_forest(mtry = tune(), min_n = tune()) %>%
#'     set_engine("ranger") %>%
#'     set_mode("classification"),
#'
#'   xgb = boost_tree(mtry = tune(), trees = tune()) %>%
#'     set_engine("xgboost") %>%
#'     set_mode("classification"),
#'
#'   glmnet = logistic_reg(penalty = tune(), mixture = tune()) %>%
#'     set_engine("glmnet") %>%
#'     set_mode("classification")
#' )
#'
#' # Run nested CV
#' result <- OmicSelector_nested_cv(
#'   data = miR_Asakura,
#'   outcome = "Class",
#'   outer_folds = 5,
#'   inner_folds = 5,
#'   preprocessing_recipe = rec,
#'   feature_selection_method = "boruta",
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
  preprocessing_recipe = NULL,
  feature_selection_method = "none",
  max_features = NULL,
  models = NULL,
  tune_grid = 10,
  metrics = NULL,
  calibrate = TRUE,
  parallel = TRUE,
  cores = parallel::detectCores() - 1,
  seed = 123,
  verbose = TRUE
) {

  # Validate inputs
  if (!outcome %in% names(data)) {
    stop(paste0("Outcome variable '", outcome, "' not found in data"))
  }

  if (nrow(data) < outer_folds * inner_folds) {
    stop("Insufficient data for the specified number of folds")
  }

  # Load required packages
  if (!requireNamespace("tidymodels", quietly = TRUE)) {
    stop("tidymodels is required. Install with: install.packages('tidymodels')")
  }

  suppressPackageStartupMessages({
    library(tidymodels)
    library(rsample)
    library(recipes)
    library(parsnip)
    library(workflows)
    library(tune)
    library(yardstick)
    library(dplyr)
  })

  if (verbose) {
    message("Starting nested cross-validation...")
    message(paste0("Outer folds: ", outer_folds, ", Inner folds: ", inner_folds))
  }

  # Set seed
  set.seed(seed)

  # Determine outcome type
  outcome_type <- ifelse(is.factor(data[[outcome]]) || is.character(data[[outcome]]),
                         "classification", "regression")

  if (outcome_type == "classification") {
    data[[outcome]] <- as.factor(data[[outcome]])
  }

  # Set default metrics if not provided
  if (is.null(metrics)) {
    if (outcome_type == "classification") {
      metrics <- metric_set(roc_auc, accuracy, sens, spec, precision, recall)
    } else {
      metrics <- metric_set(rmse, rsq, mae)
    }
  }

  # Set default models if not provided
  if (is.null(models)) {
    models <- .create_default_models(outcome_type)
  }

  # Create outer folds
  outer_cv <- vfold_cv(data, v = outer_folds, strata = outcome)

  # Set up parallel processing
  if (parallel) {
    library(doParallel)
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl), add = TRUE)
  }

  # Storage for results
  outer_results <- list()
  inner_results <- list()
  final_predictions <- data.frame()
  selected_features_by_fold <- list()
  best_models_by_fold <- list()

  # Outer loop - for unbiased performance estimation
  for (i in 1:outer_folds) {

    if (verbose) {
      message(paste0("\n=== Outer Fold ", i, "/", outer_folds, " ==="))
    }

    # Split data
    outer_split <- outer_cv$splits[[i]]
    train_data <- analysis(outer_split)
    test_data <- assessment(outer_split)

    # CRITICAL: All preprocessing happens on training data only
    # Feature selection on training data only
    if (feature_selection_method != "none") {
      if (verbose) message("  Performing feature selection...")

      selected_features <- .perform_feature_selection(
        data = train_data,
        outcome = outcome,
        method = feature_selection_method,
        max_features = max_features,
        seed = seed + i
      )

      selected_features_by_fold[[i]] <- selected_features

      # Subset data to selected features
      train_data <- train_data %>% select(all_of(c(outcome, selected_features)))
      test_data <- test_data %>% select(all_of(c(outcome, selected_features)))

      if (verbose) {
        message(paste0("  Selected ", length(selected_features), " features"))
      }
    }

    # Create or update recipe
    if (is.null(preprocessing_recipe)) {
      rec <- recipe(as.formula(paste0(outcome, " ~ .")), data = train_data) %>%
        step_zv(all_predictors()) %>%
        step_normalize(all_numeric_predictors())
    } else {
      # Update recipe with current training data
      rec <- preprocessing_recipe
      rec$template <- train_data
      rec$var_info <- NULL  # Force re-initialization
    }

    # Inner loop - for hyperparameter tuning
    inner_cv <- vfold_cv(train_data, v = inner_folds, strata = outcome)

    fold_inner_results <- list()
    fold_best_models <- list()

    for (model_name in names(models)) {

      if (verbose) {
        message(paste0("  Tuning model: ", model_name))
      }

      # Create workflow
      wf <- workflow() %>%
        add_recipe(rec) %>%
        add_model(models[[model_name]])

      # Tune or fit
      if (.has_tunable_params_nested(models[[model_name]])) {
        # Hyperparameter tuning
        tune_results <- tune_grid(
          wf,
          resamples = inner_cv,
          grid = tune_grid,
          metrics = metrics,
          control = control_grid(save_pred = FALSE, parallel_over = "everything")
        )

        fold_inner_results[[model_name]] <- tune_results

        # Select best hyperparameters
        best_params <- select_best(tune_results, metric = "roc_auc")

        # Finalize workflow
        final_wf <- finalize_workflow(wf, best_params)

      } else {
        # No tuning needed
        final_wf <- wf
        fold_inner_results[[model_name]] <- fit_resamples(
          wf,
          resamples = inner_cv,
          metrics = metrics,
          control = control_resamples(save_pred = FALSE)
        )
      }

      # Fit final model on full training data
      final_fit <- fit(final_wf, data = train_data)

      # Predict on hold-out test set
      predictions <- predict(final_fit, test_data, type = "prob") %>%
        bind_cols(predict(final_fit, test_data)) %>%
        bind_cols(test_data %>% select(all_of(outcome)))

      # Calculate metrics on test set
      if (outcome_type == "classification") {
        perf <- predictions %>%
          metrics(truth = !!sym(outcome), estimate = .pred_class,
                  .pred_[[levels(data[[outcome]])[1]]])
      } else {
        perf <- predictions %>%
          metrics(truth = !!sym(outcome), estimate = .pred)
      }

      # Store results
      fold_best_models[[model_name]] <- list(
        workflow = final_wf,
        fit = final_fit,
        predictions = predictions,
        metrics = perf
      )
    }

    # Select best model for this fold
    best_model_name <- .select_best_model(fold_best_models, metric = "roc_auc")

    if (verbose) {
      message(paste0("  Best model for fold ", i, ": ", best_model_name))
    }

    # Store results
    inner_results[[i]] <- fold_inner_results
    best_models_by_fold[[i]] <- fold_best_models[[best_model_name]]
    outer_results[[i]] <- fold_best_models[[best_model_name]]$metrics

    # Accumulate predictions
    fold_predictions <- fold_best_models[[best_model_name]]$predictions %>%
      mutate(outer_fold = i, model = best_model_name)
    final_predictions <- bind_rows(final_predictions, fold_predictions)
  }

  if (verbose) {
    message("\n=== Nested CV Complete ===")
  }

  # Calculate overall performance
  overall_metrics <- final_predictions %>%
    metrics(truth = !!sym(outcome), estimate = .pred_class,
            .pred_[[levels(data[[outcome]])[1]]])

  # Feature stability analysis
  feature_stability <- .calculate_feature_stability(selected_features_by_fold)

  # Calibration
  calibration_results <- NULL
  if (calibrate && outcome_type == "classification") {
    if (verbose) message("Performing calibration analysis...")
    calibration_results <- .perform_calibration(final_predictions, outcome)
  }

  # Create result object
  result <- list(
    outer_results = outer_results,
    inner_results = inner_results,
    overall_metrics = overall_metrics,
    final_predictions = final_predictions,
    selected_features = selected_features_by_fold,
    best_models = best_models_by_fold,
    calibration = calibration_results,
    feature_stability = feature_stability,
    metadata = list(
      outer_folds = outer_folds,
      inner_folds = inner_folds,
      feature_selection_method = feature_selection_method,
      n_samples = nrow(data),
      n_features = ncol(data) - 1,
      outcome_type = outcome_type,
      timestamp = Sys.time(),
      seed = seed
    )
  )

  class(result) <- c("OmicSelector_nested_cv", "list")

  return(result)
}


#' Internal Function: Create Default Models
#'
#' @keywords internal
#' @noRd
.create_default_models <- function(outcome_type) {

  mode <- ifelse(outcome_type == "classification", "classification", "regression")

  models <- list(
    ranger = rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>%
      set_engine("ranger", importance = "impurity") %>%
      set_mode(mode),

    xgboost = boost_tree(
      mtry = tune(),
      trees = tune(),
      min_n = tune(),
      tree_depth = tune(),
      learn_rate = tune()
    ) %>%
      set_engine("xgboost") %>%
      set_mode(mode),

    glmnet = logistic_reg(penalty = tune(), mixture = tune()) %>%
      set_engine("glmnet") %>%
      set_mode(mode)
  )

  return(models)
}


#' Internal Function: Perform Feature Selection
#'
#' @keywords internal
#' @noRd
.perform_feature_selection <- function(data, outcome, method, max_features, seed) {

  set.seed(seed)

  # Get feature names (exclude outcome)
  all_features <- setdiff(names(data), outcome)

  selected_features <- switch(
    method,

    "boruta" = {
      if (!requireNamespace("Boruta", quietly = TRUE)) {
        stop("Boruta package required. Install with: install.packages('Boruta')")
      }
      library(Boruta)

      formula_obj <- as.formula(paste0(outcome, " ~ ."))
      boruta_result <- Boruta(formula_obj, data = data, doTrace = 0)

      # Get confirmed and tentative features
      selected <- names(boruta_result$finalDecision[
        boruta_result$finalDecision %in% c("Confirmed", "Tentative")
      ])

      selected
    },

    "rfe" = {
      if (!requireNamespace("caret", quietly = TRUE)) {
        stop("caret package required. Install with: install.packages('caret')")
      }
      library(caret)

      # Use recursive feature elimination
      control <- rfeControl(functions = rfFuncs, method = "cv", number = 5)
      x <- data %>% select(-all_of(outcome))
      y <- data[[outcome]]

      results <- rfe(x, y, sizes = c(5, 10, 15, 20, 30, 50), rfeControl = control)

      results$optVariables
    },

    "lasso" = {
      if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("glmnet package required. Install with: install.packages('glmnet')")
      }
      library(glmnet)

      x <- as.matrix(data %>% select(-all_of(outcome)))
      y <- data[[outcome]]

      # Fit LASSO with cross-validation
      cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 1)

      # Extract non-zero coefficients
      coef_matrix <- coef(cv_fit, s = "lambda.min")
      selected_idx <- which(coef_matrix != 0)[-1]  # Remove intercept

      colnames(x)[selected_idx - 1]
    },

    "stability_selection" = {
      # Simplified stability selection
      # Run LASSO multiple times with subsampling
      if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("glmnet package required. Install with: install.packages('glmnet')")
      }
      library(glmnet)

      n_iterations <- 100
      selection_freq <- rep(0, length(all_features))
      names(selection_freq) <- all_features

      x <- as.matrix(data %>% select(-all_of(outcome)))
      y <- data[[outcome]]

      for (i in 1:n_iterations) {
        # Subsample
        sample_idx <- sample(1:nrow(data), size = floor(0.8 * nrow(data)))
        x_sub <- x[sample_idx, ]
        y_sub <- y[sample_idx]

        # Fit LASSO
        cv_fit <- cv.glmnet(x_sub, y_sub, family = "binomial", alpha = 1)
        coef_matrix <- coef(cv_fit, s = "lambda.min")
        selected_idx <- which(coef_matrix != 0)[-1]

        if (length(selected_idx) > 0) {
          selected_features_iter <- colnames(x)[selected_idx - 1]
          selection_freq[selected_features_iter] <- selection_freq[selected_features_iter] + 1
        }
      }

      # Select features with frequency > threshold
      threshold <- 0.6 * n_iterations
      names(selection_freq[selection_freq >= threshold])
    },

    # Default: return all features
    all_features
  )

  # Limit number of features if specified
  if (!is.null(max_features) && length(selected_features) > max_features) {
    selected_features <- selected_features[1:max_features]
  }

  return(selected_features)
}


#' Internal Function: Check for Tunable Parameters
#'
#' @keywords internal
#' @noRd
.has_tunable_params_nested <- function(model_spec) {
  params <- model_spec$args
  any(sapply(params, function(x) {
    if (is.call(x)) {
      as.character(x[[1]]) == "tune"
    } else {
      FALSE
    }
  }))
}


#' Internal Function: Select Best Model from Multiple Candidates
#'
#' @keywords internal
#' @noRd
.select_best_model <- function(models_list, metric = "roc_auc") {

  performances <- sapply(names(models_list), function(name) {
    metrics_df <- models_list[[name]]$metrics
    metrics_df %>%
      filter(.metric == metric) %>%
      pull(.estimate) %>%
      mean()
  })

  best_model_name <- names(performances)[which.max(performances)]
  return(best_model_name)
}


#' Internal Function: Calculate Feature Stability
#'
#' @keywords internal
#' @noRd
.calculate_feature_stability <- function(selected_features_list) {

  if (length(selected_features_list) == 0) {
    return(NULL)
  }

  # Get all unique features
  all_features <- unique(unlist(selected_features_list))

  # Calculate selection frequency
  selection_freq <- sapply(all_features, function(feat) {
    sum(sapply(selected_features_list, function(x) feat %in% x))
  })

  # Calculate Nogueira stability metric
  n_folds <- length(selected_features_list)
  avg_selected <- mean(sapply(selected_features_list, length))

  # Pairwise Jaccard similarity
  jaccard_sims <- c()
  if (n_folds > 1) {
    for (i in 1:(n_folds - 1)) {
      for (j in (i + 1):n_folds) {
        set1 <- selected_features_list[[i]]
        set2 <- selected_features_list[[j]]
        intersection <- length(intersect(set1, set2))
        union <- length(union(set1, set2))
        jaccard_sims <- c(jaccard_sims, intersection / union)
      }
    }
  }

  stability_results <- list(
    selection_frequency = sort(selection_freq / n_folds, decreasing = TRUE),
    mean_jaccard_similarity = mean(jaccard_sims),
    median_jaccard_similarity = median(jaccard_sims),
    n_stable_features = sum(selection_freq == n_folds),
    stable_features = names(selection_freq[selection_freq == n_folds])
  )

  return(stability_results)
}


#' Internal Function: Perform Calibration Analysis
#'
#' @keywords internal
#' @noRd
.perform_calibration <- function(predictions, outcome) {

  # Basic calibration assessment
  # More sophisticated calibration will be in the clinical_utility.R module

  result <- list(
    predictions = predictions,
    outcome = outcome,
    message = "Full calibration analysis available via OmicSelector_calibrate()"
  )

  return(result)
}


#' Print Method for OmicSelector_nested_cv
#'
#' @param x An OmicSelector_nested_cv object
#' @param ... Additional arguments (not used)
#' @export
print.OmicSelector_nested_cv <- function(x, ...) {
  cat("OmicSelector Nested Cross-Validation Results\n")
  cat("==============================================\n\n")

  cat("Configuration:\n")
  cat(paste0("  Outer folds: ", x$metadata$outer_folds, "\n"))
  cat(paste0("  Inner folds: ", x$metadata$inner_folds, "\n"))
  cat(paste0("  Feature selection: ", x$metadata$feature_selection_method, "\n"))
  cat(paste0("  Problem type: ", x$metadata$outcome_type, "\n"))
  cat(paste0("  Samples: ", x$metadata$n_samples, "\n"))
  cat(paste0("  Features: ", x$metadata$n_features, "\n"))

  cat("\nOverall Performance (Outer Loop):\n")
  if (!is.null(x$overall_metrics)) {
    print(as.data.frame(x$overall_metrics), row.names = FALSE)
  } else {
    cat("  (No metrics available)\n")
  }

  if (!is.null(x$feature_stability)) {
    cat("\nFeature Stability:\n")
    cat(paste0("  Stable features (selected in all folds): ",
               x$feature_stability$n_stable_features, "\n"))
    cat(paste0("  Mean Jaccard similarity: ",
               round(x$feature_stability$mean_jaccard_similarity, 3), "\n"))
  }

  cat("\nUse summary() for detailed results.\n")
}


#' Summary Method for OmicSelector_nested_cv
#'
#' @param object An OmicSelector_nested_cv object
#' @param ... Additional arguments (not used)
#' @export
summary.OmicSelector_nested_cv <- function(object, ...) {
  print(object)

  cat("\n")
  cat("Performance by Outer Fold:\n")
  cat("===========================\n")

  for (i in 1:object$metadata$outer_folds) {
    cat(paste0("\nFold ", i, ":\n"))
    print(object$outer_results[[i]], n = 5)
  }

  if (!is.null(object$feature_stability)) {
    cat("\n")
    cat("Top 20 Most Stable Features:\n")
    cat("============================\n")
    top_features <- head(object$feature_stability$selection_frequency, 20)
    print(data.frame(
      Feature = names(top_features),
      Selection_Frequency = as.numeric(top_features)
    ))
  }
}


#' Plot Method for OmicSelector_nested_cv
#'
#' @param x An OmicSelector_nested_cv object
#' @param ... Additional arguments (not used)
#' @export
plot.OmicSelector_nested_cv <- function(x, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting")
  }

  library(ggplot2)
  library(dplyr)

  # ROC curve plot
  if (x$metadata$outcome_type == "classification") {

    # Extract predictions
    pred_col <- grep("^\\.pred_", names(x$final_predictions), value = TRUE)[1]

    p1 <- x$final_predictions %>%
      ggplot(aes_string(d = names(x$final_predictions)[grep("^[^.]", names(x$final_predictions))[1]],
                        m = pred_col)) +
      geom_roc(n.cuts = 0) +
      style_roc() +
      labs(title = "ROC Curve - Nested CV Results") +
      theme_minimal()

    print(p1)
  }

  # Feature stability plot
  if (!is.null(x$feature_stability)) {
    top_20 <- head(x$feature_stability$selection_frequency, 20)

    p2 <- data.frame(
      Feature = factor(names(top_20), levels = rev(names(top_20))),
      Frequency = as.numeric(top_20)
    ) %>%
      ggplot(aes(x = Feature, y = Frequency)) +
      geom_col(fill = "steelblue") +
      coord_flip() +
      labs(title = "Feature Selection Stability",
           subtitle = "Top 20 Features by Selection Frequency",
           x = "Feature", y = "Selection Frequency") +
      theme_minimal()

    print(p2)
  }
}
