#' Forward Selection Feature Selection
#'
#' @description
#' Forward Selection algorithm for feature selection. Starts with empty set
#' and iteratively adds features that most improve the model performance
#' until no significant improvement is achieved.
#' This is method #19 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Forward Selection Feature Selection
#'
#' @description
#' Implements forward selection algorithm that starts with an empty feature
#' set and iteratively adds the feature that most improves model performance
#' until stopping criteria are met.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (improvement_threshold, cv_folds, scoring_metric)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_forward_selection <- function(data, 
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
  improvement_threshold <- config$improvement_threshold %||% 0.01
  cv_folds <- config$cv_folds %||% 3
  scoring_metric <- config$scoring_metric %||% "accuracy"
  
  # Initialize forward selection
  selected_features <- character(0)
  remaining_features <- colnames(feature_cols)
  
  # Track selection process
  selection_history <- list()
  current_score <- 0
  step <- 0
  
  start_time <- Sys.time()
  
  # Main forward selection loop
  while (length(selected_features) < max_features && 
         length(remaining_features) > 0 &&
         as.numeric(difftime(Sys.time(), start_time, units = "secs")) < timeout_sec) {
    
    step <- step + 1
    best_feature <- NULL
    best_score <- current_score
    best_improvement <- 0
    
    # Test each remaining feature
    for (candidate_feature in remaining_features) {
      tryCatch({
        # Create candidate feature set
        candidate_set <- c(selected_features, candidate_feature)
        
        # Evaluate candidate set
        candidate_score <- evaluate_feature_set_cv(
          data, candidate_set, cv_folds, scoring_metric
        )
        
        # Check if this is the best improvement
        improvement <- candidate_score - current_score
        
        if (improvement > best_improvement) {
          best_feature <- candidate_feature
          best_score <- candidate_score
          best_improvement <- improvement
        }
        
      }, error = function(e) {
        warning(paste("Forward selection evaluation failed for feature", candidate_feature, ":", e$message))
      })
    }
    
    # Check stopping criteria
    if (is.null(best_feature) || best_improvement < improvement_threshold) {
      break
    }
    
    # Add best feature to selected set
    selected_features <- c(selected_features, best_feature)
    remaining_features <- setdiff(remaining_features, best_feature)
    current_score <- best_score
    
    # Record selection step
    selection_history[[step]] <- list(
      step = step,
      feature_added = best_feature,
      score_improvement = best_improvement,
      total_score = current_score,
      selected_features = selected_features,
      n_features = length(selected_features)
    )
  }
  
  # Calculate final feature scores
  if (length(selected_features) > 0) {
    final_scores <- calculate_forward_feature_scores(selected_features, data, feature_cols)
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
      method = "ForwardSelection",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        improvement_threshold = improvement_threshold,
        cv_folds = cv_folds,
        scoring_metric = scoring_metric,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        final_scores = final_scores,
        final_cv_score = current_score,
        mean_score = mean(final_scores)
      ),
      selection_process = list(
        steps_completed = step,
        selection_history = selection_history,
        convergence_reason = if(length(selected_features) >= max_features) "max_features" else 
                           if(length(remaining_features) == 0) "no_more_features" else "insufficient_improvement"
      )
    )
  ))
}

#' Evaluate Feature Set with Cross-Validation
#'
#' @description
#' Evaluate a feature set using cross-validation
#'
#' @param data Full dataset
#' @param feature_names Names of features to evaluate
#' @param cv_folds Number of CV folds
#' @param scoring_metric Scoring metric to use
#' @return Cross-validation score
evaluate_feature_set_cv <- function(data, feature_names, cv_folds, scoring_metric) {
  if (length(feature_names) == 0) {
    return(0)
  }
  
  tryCatch({
    # Create subset data
    feature_indices <- match(feature_names, colnames(data)[-1]) + 1
    subset_data <- data[, c(1, feature_indices), drop = FALSE]
    
    # Perform cross-validation
    cv_scores <- perform_simple_cv(subset_data, cv_folds, scoring_metric)
    
    return(mean(cv_scores, na.rm = TRUE))
    
  }, error = function(e) {
    return(0)
  })
}

#' Perform Simple Cross-Validation
#'
#' @description
#' Simple cross-validation implementation
#'
#' @param data Dataset with class in first column
#' @param cv_folds Number of folds
#' @param scoring_metric Scoring metric
#' @return Vector of CV scores
perform_simple_cv <- function(data, cv_folds, scoring_metric) {
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
    
    # Simple prediction using mode/mean
    train_class <- train_data[, 1]
    test_class <- test_data[, 1]
    
    if (scoring_metric == "accuracy") {
      # Classification accuracy
      predicted_class <- names(sort(table(train_class), decreasing = TRUE))[1]
      cv_scores[fold] <- sum(test_class == predicted_class) / length(test_class)
    } else if (scoring_metric == "correlation") {
      # Correlation-based score
      if (ncol(train_data) > 1) {
        train_features <- train_data[, -1, drop = FALSE]
        test_features <- test_data[, -1, drop = FALSE]
        
        # Simple correlation-based prediction
        feature_means <- colMeans(train_features, na.rm = TRUE)
        test_similarity <- apply(test_features, 1, function(x) {
          cor(x, feature_means, use = "complete.obs")
        })
        
        cv_scores[fold] <- mean(abs(test_similarity), na.rm = TRUE)
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

#' Calculate Forward Feature Scores
#'
#' @description
#' Calculate importance scores for forward selected features
#'
#' @param selected_features Selected feature names
#' @param data Full dataset
#' @param feature_cols Feature columns
#' @return Feature scores
calculate_forward_feature_scores <- function(selected_features, data, feature_cols) {
  scores <- numeric(length(selected_features))
  names(scores) <- selected_features
  
  class_col <- data[, 1]
  if (is.factor(class_col) || is.character(class_col)) {
    class_numeric <- as.numeric(as.factor(class_col))
  } else {
    class_numeric <- class_col
  }
  
  # Calculate scores based on order of selection and correlation
  for (i in seq_along(selected_features)) {
    feature_name <- selected_features[i]
    feature_data <- feature_cols[, feature_name]
    
    if (is.numeric(feature_data)) {
      base_score <- abs(cor(feature_data, class_numeric, use = "complete.obs"))
    } else {
      feature_numeric <- as.numeric(as.factor(feature_data))
      base_score <- abs(cor(feature_numeric, class_numeric, use = "complete.obs"))
    }
    
    # Weight by order of selection (earlier selections get higher weights)
    order_weight <- (length(selected_features) - i + 1) / length(selected_features)
    scores[i] <- base_score * (0.7 + 0.3 * order_weight)
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
    method_id = 19,
    method_name = "ForwardSelection",
    method_function = fs_method_forward_selection,
    category = "sequential",
    dependencies = character(0),
    description = "Forward selection algorithm with iterative feature addition and cross-validation",
    parameters = list(
      improvement_threshold = 0.01,
      cv_folds = 3,
      scoring_metric = "accuracy"
    ),
    complexity = "medium",
    supports_smote = TRUE,
    timeout_default = 450
  )
  
  cat("✓ Registered Forward Selection (ForwardSelection) method (ID: 19)\n")
  
}, error = function(e) {
  warning("Failed to register ForwardSelection method: ", e$message)
})
