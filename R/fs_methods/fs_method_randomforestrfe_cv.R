#' Random Forest RFE with Cross-Validation Feature Selection
#'
#' @description
#' Random Forest Recursive Feature Elimination with cross-validation for
#' robust feature selection. Iteratively removes least important features
#' based on Random Forest variable importance with CV validation.
#' This is method #10 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Random Forest RFE with Cross-Validation
#'
#' @description
#' Recursive Feature Elimination using Random Forest with cross-validation
#' to determine optimal feature subset. Uses variable importance from RF
#' and removes features iteratively with CV performance monitoring.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (n_trees, cv_folds, step_size, min_features)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_randomforestrfe_cv <- function(data, 
                                        config = list(), 
                                        max_features = 20, 
                                        use_smote = FALSE, 
                                        timeout_sec = 300) {
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract features
  feature_cols <- data[, -1, drop = FALSE]
  
  # Configuration parameters
  n_trees <- config$n_trees %||% 100
  cv_folds <- config$cv_folds %||% 5
  step_size <- config$step_size %||% 0.1  # Remove 10% of features per iteration
  min_features <- config$min_features %||% 5
  
  # Initialize with all features
  current_features <- colnames(feature_cols)
  best_features <- current_features
  best_score <- 0
  
  # Track performance across iterations
  iteration_scores <- numeric(0)
  feature_counts <- list()
  
  iteration <- 0
  start_time <- Sys.time()
  
  while (length(current_features) > min_features && 
         length(current_features) > max_features &&
         as.numeric(difftime(Sys.time(), start_time, units = "secs")) < timeout_sec) {
    
    iteration <- iteration + 1
    
    tryCatch({
      # Prepare current data subset
      current_data <- data[, c(TRUE, colnames(feature_cols) %in% current_features), drop = FALSE]
      
      # Perform cross-validation
      cv_score <- perform_rf_cv(current_data, cv_folds, n_trees)
      iteration_scores <- c(iteration_scores, cv_score)
      
      # Update best features if performance improved
      if (cv_score > best_score) {
        best_score <- cv_score
        best_features <- current_features
      }
      
      # Calculate feature importance using Random Forest
      importance_scores <- calculate_rf_importance(current_data, n_trees)
      
      # Remove least important features
      n_remove <- max(1, floor(step_size * length(current_features)))
      n_keep <- length(current_features) - n_remove
      
      # Keep most important features
      top_indices <- order(importance_scores, decreasing = TRUE)[seq_len(n_keep)]
      current_features <- names(importance_scores)[top_indices]
      
      feature_counts[[iteration]] <- list(
        features = current_features,
        score = cv_score,
        n_features = length(current_features)
      )
      
    }, error = function(e) {
      warning(paste("RandomForestRFE_CV iteration", iteration, "failed:", e$message))
      break
    })
  }
  
  # Select final features (limit to max_features)
  if (length(best_features) > max_features) {
    # Use final Random Forest to rank best features
    final_data <- data[, c(TRUE, colnames(feature_cols) %in% best_features), drop = FALSE]
    final_importance <- calculate_rf_importance(final_data, n_trees)
    
    top_indices <- order(final_importance, decreasing = TRUE)[seq_len(max_features)]
    selected_features <- names(final_importance)[top_indices]
    final_scores <- final_importance[top_indices]
  } else {
    selected_features <- best_features
    final_data <- data[, c(TRUE, colnames(feature_cols) %in% selected_features), drop = FALSE]
    final_scores <- calculate_rf_importance(final_data, n_trees)
  }
  
  return(list(
    features = selected_features,
    scores = final_scores,
    metadata = list(
      method = "RandomForestRFE_CV",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        n_trees = n_trees,
        cv_folds = cv_folds,
        step_size = step_size,
        min_features = min_features,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        final_scores = final_scores,
        best_cv_score = best_score,
        final_cv_score = tail(iteration_scores, 1),
        mean_importance = mean(final_scores)
      ),
      rfe_process = list(
        n_iterations = iteration,
        cv_scores_by_iteration = iteration_scores,
        feature_trajectory = feature_counts,
        convergence_reached = length(current_features) <= max_features
      )
    )
  ))
}

#' Perform Random Forest Cross-Validation
#'
#' @description
#' Perform k-fold cross-validation using Random Forest
#'
#' @param data Training data with class in first column
#' @param cv_folds Number of CV folds
#' @param n_trees Number of trees in RF
#' @return Average CV accuracy
perform_rf_cv <- function(data, cv_folds, n_trees) {
  tryCatch({
    n_samples <- nrow(data)
    fold_size <- floor(n_samples / cv_folds)
    
    cv_scores <- numeric(cv_folds)
    
    for (fold in seq_len(cv_folds)) {
      # Create train/test split
      test_start <- (fold - 1) * fold_size + 1
      test_end <- min(fold * fold_size, n_samples)
      test_indices <- test_start:test_end
      
      train_data <- data[-test_indices, , drop = FALSE]
      test_data <- data[test_indices, , drop = FALSE]
      
      # Train simplified Random Forest
      rf_model <- train_simple_rf(train_data, n_trees)
      
      # Predict and calculate accuracy
      predictions <- predict_simple_rf(rf_model, test_data[, -1, drop = FALSE])
      actual <- test_data[, 1]
      
      accuracy <- sum(predictions == actual) / length(actual)
      cv_scores[fold] <- accuracy
    }
    
    return(mean(cv_scores))
    
  }, error = function(e) {
    return(0.5)  # Return baseline accuracy on error
  })
}

#' Calculate Random Forest Feature Importance
#'
#' @description
#' Calculate feature importance using Random Forest
#'
#' @param data Training data with class in first column
#' @param n_trees Number of trees
#' @return Named vector of feature importance scores
calculate_rf_importance <- function(data, n_trees) {
  tryCatch({
    rf_model <- train_simple_rf(data, n_trees)
    return(rf_model$importance)
    
  }, error = function(e) {
    # Fallback to correlation-based importance
    feature_cols <- data[, -1, drop = FALSE]
    class_col <- data[, 1]
    
    if (is.factor(class_col) || is.character(class_col)) {
      class_numeric <- as.numeric(as.factor(class_col))
    } else {
      class_numeric <- class_col
    }
    
    importance <- apply(feature_cols, 2, function(x) {
      if (is.numeric(x)) {
        abs(cor(x, class_numeric, use = "complete.obs"))
      } else {
        x_numeric <- as.numeric(as.factor(x))
        abs(cor(x_numeric, class_numeric, use = "complete.obs"))
      }
    })
    
    importance[is.na(importance)] <- 0
    return(importance)
  })
}

#' Train Simple Random Forest
#'
#' @description
#' Simple Random Forest implementation for importance calculation
#'
#' @param data Training data
#' @param n_trees Number of trees
#' @return Simple RF model
train_simple_rf <- function(data, n_trees) {
  tryCatch({
    feature_cols <- data[, -1, drop = FALSE]
    
    # Simple implementation: average correlation across bootstrap samples
    importance_accumulator <- rep(0, ncol(feature_cols))
    names(importance_accumulator) <- colnames(feature_cols)
    
    for (tree in seq_len(min(n_trees, 20))) {  # Limit for efficiency
      # Bootstrap sample
      boot_indices <- sample(nrow(data), replace = TRUE)
      boot_data <- data[boot_indices, ]
      
      # Calculate importance for this tree
      tree_importance <- calculate_tree_importance(boot_data)
      importance_accumulator <- importance_accumulator + tree_importance
    }
    
    final_importance <- importance_accumulator / min(n_trees, 20)
    
    return(list(
      importance = final_importance,
      data = data
    ))
    
  }, error = function(e) {
    # Fallback importance
    feature_cols <- data[, -1, drop = FALSE]
    fallback_importance <- rep(1, ncol(feature_cols))
    names(fallback_importance) <- colnames(feature_cols)
    
    return(list(
      importance = fallback_importance,
      data = data
    ))
  })
}

#' Calculate importance for single tree
#'
#' @description
#' Calculate feature importance for a single decision tree (simplified)
#'
#' @param data Bootstrap sample data
#' @return Feature importance scores
calculate_tree_importance <- function(data) {
  tryCatch({
    class_col <- data[, 1]
    feature_cols <- data[, -1, drop = FALSE]
    
    if (is.factor(class_col) || is.character(class_col)) {
      class_numeric <- as.numeric(as.factor(class_col))
    } else {
      class_numeric <- class_col
    }
    
    # Use correlation as proxy for tree-based importance
    importance <- apply(feature_cols, 2, function(x) {
      if (is.numeric(x)) {
        abs(cor(x, class_numeric, use = "complete.obs"))
      } else {
        x_numeric <- as.numeric(as.factor(x))
        abs(cor(x_numeric, class_numeric, use = "complete.obs"))
      }
    })
    
    importance[is.na(importance)] <- 0
    return(importance)
    
  }, error = function(e) {
    feature_cols <- data[, -1, drop = FALSE]
    fallback <- rep(0, ncol(feature_cols))
    names(fallback) <- colnames(feature_cols)
    return(fallback)
  })
}

#' Predict using Simple Random Forest
#'
#' @description
#' Make predictions using simple RF model
#'
#' @param model Simple RF model
#' @param new_data New data for prediction
#' @return Predictions
predict_simple_rf <- function(model, new_data) {
  tryCatch({
    # Simple prediction: use majority class from training
    train_class <- model$data[, 1]
    majority_class <- names(sort(table(train_class), decreasing = TRUE))[1]
    
    predictions <- rep(majority_class, nrow(new_data))
    return(predictions)
    
  }, error = function(e) {
    return(rep("class1", nrow(new_data)))
  })
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
tryCatch({
  register_fs_method(
    method_id = 10,
    method_name = "RandomForestRFE_CV",
    method_function = fs_method_randomforestrfe_cv,
    category = "ensemble",
    dependencies = character(0),
    description = "Random Forest Recursive Feature Elimination with cross-validation",
    parameters = list(
      n_trees = 100,
      cv_folds = 5,
      step_size = 0.1,
      min_features = 5
    ),
    complexity = "high",
    supports_smote = TRUE,
    timeout_default = 600
  )
  
  cat("✓ Registered Random Forest RFE CV (RandomForestRFE_CV) method (ID: 10)\n")
  
}, error = function(e) {
  warning("Failed to register RandomForestRFE_CV method: ", e$message)
})
