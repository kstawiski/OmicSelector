#' Classification Loop Feature Selection
#'
#' @description
#' Feature selection using classification performance in a loop with
#' Random Forest importance as the core selection mechanism.
#' This is method #4 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Classification Loop Feature Selection
#'
#' @description
#' Iterative feature selection using classification performance feedback.
#' Uses Random Forest importance with cross-validation to select features.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (ntree, cv_folds, iterations)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_classloop <- function(data, 
                               config = list(), 
                               max_features = 20, 
                               use_smote = FALSE, 
                               timeout_sec = 300) {
  
  # Check if required packages are available
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    warning("randomForest package not available for classloop")
    return(NULL)
  }
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract classes and features
  classes <- as.factor(data[, 1])
  features <- data[, -1, drop = FALSE]
  
  # Configuration parameters
  ntree <- config$ntree %||% 100
  cv_folds <- config$cv_folds %||% 5
  iterations <- config$iterations %||% 3
  
  # Initialize results
  feature_importance_scores <- numeric(ncol(features))
  names(feature_importance_scores) <- colnames(features)
  
  # Iterative feature selection loop
  for (iter in seq_len(iterations)) {
    
    tryCatch({
      # Train Random Forest model
      rf_model <- randomForest::randomForest(
        x = features,
        y = classes,
        ntree = ntree,
        importance = TRUE,
        keep.forest = FALSE
      )
      
      # Get importance scores
      importance_scores <- randomForest::importance(rf_model)[, "MeanDecreaseGini"]
      
      # Accumulate scores across iterations
      feature_importance_scores <- feature_importance_scores + importance_scores
      
    }, error = function(e) {
      warning("Iteration ", iter, " failed: ", e$message)
    })
  }
  
  # Average scores across iterations
  feature_importance_scores <- feature_importance_scores / iterations
  
  # Handle edge cases
  if (all(feature_importance_scores == 0) || length(feature_importance_scores) == 0) {
    warning("No features scored by classloop")
    return(NULL)
  }
  
  # Select top features
  n_select <- min(max_features, length(feature_importance_scores))
  top_indices <- order(feature_importance_scores, decreasing = TRUE)[seq_len(n_select)]
  selected_features <- names(feature_importance_scores)[top_indices]
  
  # Optional: Cross-validation performance assessment
  cv_performance <- NULL
  if (length(selected_features) > 0) {
    tryCatch({
      # Quick CV assessment of selected features
      selected_data <- data.frame(class = classes, features[, selected_features, drop = FALSE])
      
      # Simple cross-validation
      fold_size <- floor(nrow(selected_data) / cv_folds)
      accuracies <- numeric(cv_folds)
      
      for (fold in seq_len(cv_folds)) {
        start_idx <- (fold - 1) * fold_size + 1
        end_idx <- min(fold * fold_size, nrow(selected_data))
        
        if (end_idx > start_idx) {
          test_indices <- start_idx:end_idx
          train_indices <- setdiff(seq_len(nrow(selected_data)), test_indices)
          
          train_data <- selected_data[train_indices, ]
          test_data <- selected_data[test_indices, ]
          
          # Train simple model
          cv_model <- randomForest::randomForest(
            class ~ ., 
            data = train_data, 
            ntree = 50
          )
          
          # Predict and calculate accuracy
          predictions <- predict(cv_model, test_data)
          accuracies[fold] <- mean(predictions == test_data$class)
        }
      }
      
      cv_performance <- list(
        mean_accuracy = mean(accuracies, na.rm = TRUE),
        sd_accuracy = sd(accuracies, na.rm = TRUE),
        fold_accuracies = accuracies
      )
      
    }, error = function(e) {
      cv_performance <- list(error = e$message)
    })
  }
  
  return(list(
    features = selected_features,
    scores = feature_importance_scores[top_indices],
    metadata = list(
      method = "classloop",
      n_features_input = ncol(features),
      n_features_selected = length(selected_features),
      parameters_used = list(
        ntree = ntree,
        cv_folds = cv_folds,
        iterations = iterations,
        max_features = max_features,
        use_smote = use_smote
      ),
      performance_assessment = cv_performance,
      importance_statistics = list(
        mean_importance = mean(feature_importance_scores[top_indices]),
        max_importance = max(feature_importance_scores[top_indices]),
        min_importance = min(feature_importance_scores[top_indices])
      )
    )
  ))
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
tryCatch({
  register_fs_method(
    method_id = 4,
    method_name = "classloop",
    method_function = fs_method_classloop,
    category = "wrapper",
    dependencies = c("randomForest"),
    description = "Feature selection using classification performance feedback loop",
    parameters = list(
      ntree = 100,
      cv_folds = 5,
      iterations = 3
    ),
    complexity = "medium",
    supports_smote = TRUE,
    timeout_default = 300
  )
  
  cat("✓ Registered Classification Loop (classloop) method (ID: 4)\n")
  
}, error = function(e) {
  warning("Failed to register classloop method: ", e$message)
})
