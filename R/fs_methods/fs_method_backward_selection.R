#' Backward Selection Feature Selection
#'
#' @description
#' Backward Selection algorithm for feature selection. Starts with all features
#' and iteratively removes features that least impact model performance
#' until stopping criteria are met.
#' This is method #20 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Backward Selection Feature Selection
#'
#' @description
#' Implements backward selection algorithm that starts with all features
#' and iteratively removes the feature that least impacts model performance
#' until stopping criteria are met.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (degradation_threshold, cv_folds, scoring_metric)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_backward_selection <- function(data, 
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
  degradation_threshold <- config$degradation_threshold %||% 0.02
  cv_folds <- config$cv_folds %||% 3
  scoring_metric <- config$scoring_metric %||% "accuracy"
  min_features <- config$min_features %||% 2
  
  # Initialize backward selection with all features
  current_features <- colnames(feature_cols)
  
  # Track selection process
  elimination_history <- list()
  step <- 0
  
  # Get initial performance with all features
  current_score <- evaluate_feature_set_backward_cv(
    data, current_features, cv_folds, scoring_metric
  )
  
  start_time <- Sys.time()
  
  # Main backward selection loop
  while (length(current_features) > max_features && 
         length(current_features) > min_features &&
         as.numeric(difftime(Sys.time(), start_time, units = "secs")) < timeout_sec) {
    
    step <- step + 1
    worst_feature <- NULL
    best_score_after_removal <- -Inf
    smallest_degradation <- Inf
    
    # Test removal of each current feature
    for (candidate_feature in current_features) {
      tryCatch({
        # Create feature set without candidate
        reduced_set <- setdiff(current_features, candidate_feature)
        
        if (length(reduced_set) == 0) next
        
        # Evaluate reduced set
        reduced_score <- evaluate_feature_set_backward_cv(
          data, reduced_set, cv_folds, scoring_metric
        )
        
        # Calculate performance degradation
        degradation <- current_score - reduced_score
        
        # Find feature with smallest degradation (least important)
        if (degradation < smallest_degradation) {
          worst_feature <- candidate_feature
          best_score_after_removal <- reduced_score
          smallest_degradation <- degradation
        }
        
      }, error = function(e) {
        warning(paste("Backward selection evaluation failed for removing feature", candidate_feature, ":", e$message))
      })
    }
    
    # Check stopping criteria
    if (is.null(worst_feature) || smallest_degradation > degradation_threshold) {
      break
    }
    
    # Remove worst feature
    current_features <- setdiff(current_features, worst_feature)
    current_score <- best_score_after_removal
    
    # Record elimination step
    elimination_history[[step]] <- list(
      step = step,
      feature_removed = worst_feature,
      score_degradation = smallest_degradation,
      score_after_removal = current_score,
      remaining_features = current_features,
      n_features = length(current_features)
    )
  }
  
  # Final feature selection
  selected_features <- current_features
  
  # Calculate final feature scores
  if (length(selected_features) > 0) {
    final_scores <- calculate_backward_feature_scores(selected_features, data, feature_cols, elimination_history)
  } else {
    # Fallback to correlation-based selection
    class_col <- data[, 1]
    if (is.factor(class_col) || is.character(class_col)) {
      class_numeric <- as.numeric(as.factor(class_col))
    } else {
      class_numeric <- class_col
    }
    
    correlations <- abs(cor(feature_cols, class_numeric, use = "complete.obs"))
    correlations[is.na(correlations)] <- 0
    
    n_select <- min(max_features, length(correlations))
    top_indices <- order(correlations, decreasing = TRUE)[seq_len(n_select)]
    
    selected_features <- colnames(feature_cols)[top_indices]
    final_scores <- correlations[top_indices]
    names(final_scores) <- selected_features
  }
  
  return(list(
    features = selected_features,
    scores = final_scores,
    metadata = list(
      method = "BackwardSelection",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        degradation_threshold = degradation_threshold,
        cv_folds = cv_folds,
        scoring_metric = scoring_metric,
        min_features = min_features,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        final_scores = final_scores,
        final_cv_score = current_score,
        mean_score = mean(final_scores)
      ),
      elimination_process = list(
        steps_completed = step,
        elimination_history = elimination_history,
        features_eliminated = ncol(feature_cols) - length(selected_features),
        convergence_reason = if(length(selected_features) <= max_features) "target_reached" else 
                           if(length(selected_features) <= min_features) "min_features" else "degradation_threshold"
      )
    )
  ))
}

#' Evaluate Feature Set with Cross-Validation (Backward)
#'
#' @description
#' Evaluate a feature set using cross-validation for backward selection
#'
#' @param data Full dataset
#' @param feature_names Names of features to evaluate
#' @param cv_folds Number of CV folds
#' @param scoring_metric Scoring metric to use
#' @return Cross-validation score
evaluate_feature_set_backward_cv <- function(data, feature_names, cv_folds, scoring_metric) {
  if (length(feature_names) == 0) {
    return(0)
  }
  
  tryCatch({
    # Create subset data
    feature_indices <- match(feature_names, colnames(data)[-1]) + 1
    feature_indices <- feature_indices[!is.na(feature_indices)]
    
    if (length(feature_indices) == 0) {
      return(0)
    }
    
    subset_data <- data[, c(1, feature_indices), drop = FALSE]
    
    # Perform cross-validation
    cv_scores <- perform_backward_cv(subset_data, cv_folds, scoring_metric)
    
    return(mean(cv_scores, na.rm = TRUE))
    
  }, error = function(e) {
    return(0)
  })
}

#' Perform Cross-Validation for Backward Selection
#'
#' @description
#' Cross-validation implementation optimized for backward selection
#'
#' @param data Dataset with class in first column
#' @param cv_folds Number of folds
#' @param scoring_metric Scoring metric
#' @return Vector of CV scores
perform_backward_cv <- function(data, cv_folds, scoring_metric) {
  n_samples <- nrow(data)
  fold_size <- floor(n_samples / cv_folds)
  cv_scores <- numeric(cv_folds)
  
  for (fold in seq_len(cv_folds)) {
    # Create train/test split
    test_start <- (fold - 1) * fold_size + 1
    test_end <- min(fold * fold_size, n_samples)
    test_indices <- test_start:test_end
    
    if (length(test_indices) == 0) next
    
    train_data <- data[-test_indices, , drop = FALSE]
    test_data <- data[test_indices, , drop = FALSE]
    
    if (nrow(train_data) == 0 || nrow(test_data) == 0) next
    
    # Enhanced prediction based on features
    train_class <- train_data[, 1]
    test_class <- test_data[, 1]
    
    if (scoring_metric == "accuracy") {
      # Classification with feature-based similarity
      if (ncol(train_data) > 1) {
        train_features <- train_data[, -1, drop = FALSE]
        test_features <- test_data[, -1, drop = FALSE]
        
        # Simple k-NN style prediction (k=1 for speed)
        predictions <- apply(test_features, 1, function(test_row) {
          # Find most similar training example
          distances <- apply(train_features, 1, function(train_row) {
            if (is.numeric(test_row) && is.numeric(train_row)) {
              sqrt(sum((test_row - train_row)^2, na.rm = TRUE))
            } else {
              # For mixed data, use simple matching
              sum(test_row != train_row, na.rm = TRUE)
            }
          })
          
          closest_idx <- which.min(distances)
          return(train_class[closest_idx])
        })
        
        cv_scores[fold] <- sum(predictions == test_class) / length(test_class)
      } else {
        # Fallback to majority class
        predicted_class <- names(sort(table(train_class), decreasing = TRUE))[1]
        cv_scores[fold] <- sum(test_class == predicted_class) / length(test_class)
      }
      
    } else if (scoring_metric == "correlation") {
      # Multi-feature correlation score
      if (ncol(train_data) > 1) {
        train_features <- train_data[, -1, drop = FALSE]
        test_features <- test_data[, -1, drop = FALSE]
        
        # Calculate feature-class correlations
        if (is.factor(train_class) || is.character(train_class)) {
          train_class_numeric <- as.numeric(as.factor(train_class))
        } else {
          train_class_numeric <- train_class
        }
        
        feature_correlations <- apply(train_features, 2, function(x) {
          if (is.numeric(x)) {
            abs(cor(x, train_class_numeric, use = "complete.obs"))
          } else {
            x_numeric <- as.numeric(as.factor(x))
            abs(cor(x_numeric, train_class_numeric, use = "complete.obs"))
          }
        })
        
        feature_correlations[is.na(feature_correlations)] <- 0
        cv_scores[fold] <- mean(feature_correlations)
      } else {
        cv_scores[fold] <- 0
      }
      
    } else {
      # Default to accuracy
      predicted_class <- names(sort(table(train_class), decreasing = TRUE))[1]
      cv_scores[fold] <- sum(test_class == predicted_class) / length(test_class)
    }
  }
  
  return(cv_scores)
}

#' Calculate Backward Feature Scores
#'
#' @description
#' Calculate importance scores for backward selected features
#'
#' @param selected_features Selected feature names
#' @param data Full dataset
#' @param feature_cols Feature columns
#' @param elimination_history History of elimination process
#' @return Feature scores
calculate_backward_feature_scores <- function(selected_features, data, feature_cols, elimination_history) {
  scores <- numeric(length(selected_features))
  names(scores) <- selected_features
  
  class_col <- data[, 1]
  if (is.factor(class_col) || is.character(class_col)) {
    class_numeric <- as.numeric(as.factor(class_col))
  } else {
    class_numeric <- class_col
  }
  
  # Calculate base scores (correlation with class)
  for (i in seq_along(selected_features)) {
    feature_name <- selected_features[i]
    feature_data <- feature_cols[, feature_name]
    
    if (is.numeric(feature_data)) {
      base_score <- abs(cor(feature_data, class_numeric, use = "complete.obs"))
    } else {
      feature_numeric <- as.numeric(as.factor(feature_data))
      base_score <- abs(cor(feature_numeric, class_numeric, use = "complete.obs"))
    }
    
    # Boost score based on survival in elimination process
    # Features that survived longer (weren't eliminated) get higher scores
    survival_bonus <- 1.0  # All selected features survived completely
    
    # Check if feature was ever considered for elimination but survived
    if (length(elimination_history) > 0) {
      elimination_steps <- sapply(elimination_history, function(x) x$feature_removed)
      if (!(feature_name %in% elimination_steps)) {
        # Feature was never eliminated - boost its score
        survival_bonus <- 1.2
      }
    }
    
    scores[i] <- base_score * survival_bonus
  }
  
  scores[is.na(scores)] <- 0.1
  
  # Normalize to [0, 1]
  max_score <- max(scores)
  if (max_score > 0) {
    scores <- scores / max_score
  }
  
  return(scores)
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
tryCatch({
  register_fs_method(
    method_id = 20,
    method_name = "BackwardSelection",
    method_function = fs_method_backward_selection,
    category = "sequential",
    dependencies = character(0),
    description = "Backward selection algorithm with iterative feature elimination and cross-validation",
    parameters = list(
      degradation_threshold = 0.02,
      cv_folds = 3,
      scoring_metric = "accuracy",
      min_features = 2
    ),
    complexity = "medium",
    supports_smote = TRUE,
    timeout_default = 450
  )
  
  cat("✓ Registered Backward Selection (BackwardSelection) method (ID: 20)\n")
  
}, error = function(e) {
  warning("Failed to register BackwardSelection method: ", e$message)
})
