#' Random Forest Feature Importance Selection Method
#'
#' Implements feature selection using Random Forest variable importance measures.
#' This method combines multiple decision trees and ranks features based on
#' their importance across the ensemble.
#'
#' @param ramwas.data Expression data matrix (samples as rows, features as columns)  
#' @param target Target variable vector
#' @param smote_sampling Logical, whether to apply SMOTE sampling
#' @return List containing selected features and method metadata
#' @export
random_forest_importance_fs <- function(ramwas.data, target, smote_sampling = FALSE) {
  
  # Start timing
  start_time <- Sys.time()
  
  # Initialize results with error handling
  tryCatch({
    
    # Validate inputs
    if (!is.matrix(ramwas.data) && !is.data.frame(ramwas.data)) {
      stop("Input data must be a matrix or data.frame")
    }
    if (length(target) != nrow(ramwas.data)) {
      stop("Target length must match number of data rows")
    }
    
    ramwas.data <- as.matrix(ramwas.data)
    
    # Apply SMOTE if requested
    if (smote_sampling && length(unique(target)) > 1) {
      smote_result <- perform_smote_sampling(ramwas.data, target)
      ramwas.data <- smote_result$data
      target <- smote_result$target
    }
    
    # Remove constant and highly correlated features
    preprocessed <- preprocess_features_basic(ramwas.data, target)
    data_clean <- preprocessed$data
    
    if (ncol(data_clean) == 0) {
      return(create_empty_result("No features remained after preprocessing"))
    }
    
    # Perform Random Forest feature selection
    rf_result <- perform_random_forest_selection(data_clean, target)
    
    # Extract selected features
    if (length(rf_result$selected_features) == 0) {
      return(create_empty_result("No features selected by Random Forest"))
    }
    
    # Calculate timing
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      selected_features = rf_result$selected_features,
      feature_scores = rf_result$feature_scores,
      method_name = "Random_Forest_Importance",
      execution_time = execution_time,
      n_features_selected = length(rf_result$selected_features),
      n_features_input = ncol(ramwas.data),
      method_params = list(
        n_trees = rf_result$n_trees,
        mtry = rf_result$mtry,
        importance_type = rf_result$importance_type,
        selection_threshold = rf_result$selection_threshold,
        smote_sampling = smote_sampling
      ),
      preprocessing_info = list(
        constant_features_removed = preprocessed$constant_removed,
        correlated_features_removed = preprocessed$correlated_removed,
        final_features_count = ncol(data_clean)
      ),
      rf_info = rf_result$rf_info
    ))
    
  }, error = function(e) {
    # Return error result
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      selected_features = character(0),
      feature_scores = numeric(0),
      method_name = "Random_Forest_Importance",
      execution_time = execution_time,
      n_features_selected = 0,
      n_features_input = ncol(ramwas.data),
      error_message = paste("Random Forest failed:", e$message),
      method_params = list(smote_sampling = smote_sampling)
    ))
  })
}

#' Perform Random Forest Feature Selection
#'
#' Implements Random Forest with multiple importance measures
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with selected features and RF information
perform_random_forest_selection <- function(data_matrix, target_vector) {
  
  n_features <- ncol(data_matrix)
  
  # Determine if classification or regression
  is_classification <- determine_task_type(target_vector)
  
  # Set Random Forest parameters
  n_trees <- max(100, min(500, round(sqrt(n_features) * 10)))
  mtry <- max(1, round(sqrt(n_features)))
  
  # Build Random Forest ensemble
  rf_ensemble <- build_random_forest_ensemble(data_matrix, target_vector, n_trees, mtry, is_classification)
  
  # Calculate multiple importance measures
  importance_measures <- calculate_rf_importance_measures(rf_ensemble, data_matrix, target_vector, is_classification)
  
  # Combine importance measures
  combined_importance <- combine_importance_measures(importance_measures)
  
  # Select features based on combined importance
  selected_result <- select_features_by_importance(combined_importance, n_features)
  
  return(list(
    selected_features = selected_result$selected_features,
    feature_scores = selected_result$feature_scores,
    n_trees = n_trees,
    mtry = mtry,
    importance_type = "combined",
    selection_threshold = selected_result$threshold,
    rf_info = list(
      task_type = if(is_classification) "classification" else "regression",
      oob_error = rf_ensemble$oob_error,
      importance_measures = importance_measures
    )
  ))
}

#' Build Random Forest Ensemble
#'
#' Creates Random Forest ensemble with bootstrap sampling
#'
#' @param X Input data matrix
#' @param y Target vector
#' @param n_trees Number of trees
#' @param mtry Number of features to consider at each split
#' @param is_classification Whether this is a classification task
#' @return Random Forest ensemble object
build_random_forest_ensemble <- function(X, y, n_trees, mtry, is_classification) {
  
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  # Initialize ensemble
  trees <- list()
  oob_predictions <- matrix(0, nrow = n_samples, ncol = if(is_classification) length(unique(y)) else 1)
  oob_counts <- integer(n_samples)
  feature_importance <- numeric(n_features)
  
  # Build trees
  for (tree_idx in 1:n_trees) {
    
    # Bootstrap sampling
    bootstrap_idx <- sample(n_samples, n_samples, replace = TRUE)
    oob_idx <- setdiff(1:n_samples, unique(bootstrap_idx))
    
    X_bootstrap <- X[bootstrap_idx, , drop = FALSE]
    y_bootstrap <- y[bootstrap_idx]
    
    # Build decision tree
    tree <- build_decision_tree(X_bootstrap, y_bootstrap, mtry, is_classification)
    trees[[tree_idx]] <- tree
    
    # Out-of-bag predictions
    if (length(oob_idx) > 0) {
      oob_pred <- predict_tree(tree, X[oob_idx, , drop = FALSE], is_classification)
      
      if (is_classification) {
        for (i in seq_along(oob_idx)) {
          class_idx <- which(unique(y) == oob_pred[i])
          if (length(class_idx) > 0) {
            oob_predictions[oob_idx[i], class_idx] <- oob_predictions[oob_idx[i], class_idx] + 1
          }
        }
      } else {
        oob_predictions[oob_idx, 1] <- oob_predictions[oob_idx, 1] + oob_pred
      }
      
      oob_counts[oob_idx] <- oob_counts[oob_idx] + 1
    }
    
    # Accumulate feature importance
    feature_importance <- feature_importance + tree$feature_importance
  }
  
  # Calculate OOB error
  oob_error <- calculate_oob_error(oob_predictions, oob_counts, y, is_classification)
  
  # Average feature importance
  feature_importance <- feature_importance / n_trees
  
  return(list(
    trees = trees,
    feature_importance = feature_importance,
    oob_error = oob_error,
    n_trees = n_trees,
    mtry = mtry
  ))
}

#' Build Decision Tree
#'
#' Builds a single decision tree for the random forest
#'
#' @param X Training data
#' @param y Training labels
#' @param mtry Number of features to consider
#' @param is_classification Classification flag
#' @param max_depth Maximum tree depth
#' @return Decision tree object
build_decision_tree <- function(X, y, mtry, is_classification, max_depth = 10) {
  
  n_features <- ncol(X)
  
  # Initialize tree structure
  tree <- list(
    nodes = list(),
    feature_importance = numeric(n_features)
  )
  
  # Build tree recursively
  root_node <- build_tree_node(X, y, 1:n_features, mtry, is_classification, 1, max_depth, tree)
  tree$nodes[[1]] <- root_node
  
  return(tree)
}

#' Build Tree Node
#'
#' Recursively builds tree nodes
#'
#' @param X Node data
#' @param y Node labels
#' @param available_features Available features for splitting
#' @param mtry Features to consider
#' @param is_classification Classification flag
#' @param depth Current depth
#' @param max_depth Maximum depth
#' @param tree Tree object for importance tracking
#' @return Tree node
build_tree_node <- function(X, y, available_features, mtry, is_classification, depth, max_depth, tree) {
  
  n_samples <- nrow(X)
  
  # Stopping criteria
  if (depth >= max_depth || n_samples < 5 || length(unique(y)) == 1) {
    return(list(
      type = "leaf",
      prediction = if(is_classification) get_majority_class(y) else mean(y),
      n_samples = n_samples
    ))
  }
  
  # Select random subset of features
  selected_features <- sample(available_features, min(mtry, length(available_features)))
  
  # Find best split
  best_split <- find_best_split(X, y, selected_features, is_classification)
  
  if (is.null(best_split)) {
    return(list(
      type = "leaf",
      prediction = if(is_classification) get_majority_class(y) else mean(y),
      n_samples = n_samples
    ))
  }
  
  # Update feature importance
  tree$feature_importance[best_split$feature] <- tree$feature_importance[best_split$feature] + best_split$importance
  
  # Split data
  left_idx <- which(X[, best_split$feature] <= best_split$threshold)
  right_idx <- which(X[, best_split$feature] > best_split$threshold)
  
  # Build child nodes
  left_node <- build_tree_node(
    X[left_idx, , drop = FALSE], y[left_idx], available_features, 
    mtry, is_classification, depth + 1, max_depth, tree
  )
  
  right_node <- build_tree_node(
    X[right_idx, , drop = FALSE], y[right_idx], available_features, 
    mtry, is_classification, depth + 1, max_depth, tree
  )
  
  return(list(
    type = "internal",
    feature = best_split$feature,
    threshold = best_split$threshold,
    left = left_node,
    right = right_node,
    importance = best_split$importance,
    n_samples = n_samples
  ))
}

#' Find Best Split
#'
#' Finds the best feature and threshold for splitting
#'
#' @param X Node data
#' @param y Node labels
#' @param features Features to consider
#' @param is_classification Classification flag
#' @return Best split information
find_best_split <- function(X, y, features, is_classification) {
  
  best_score <- -Inf
  best_split <- NULL
  
  for (feature in features) {
    
    feature_values <- X[, feature]
    unique_values <- sort(unique(feature_values))
    
    if (length(unique_values) < 2) next
    
    # Try different thresholds
    for (i in 1:(length(unique_values) - 1)) {
      threshold <- (unique_values[i] + unique_values[i + 1]) / 2
      
      left_idx <- which(feature_values <= threshold)
      right_idx <- which(feature_values > threshold)
      
      if (length(left_idx) == 0 || length(right_idx) == 0) next
      
      # Calculate split quality
      split_score <- calculate_split_score(y, left_idx, right_idx, is_classification)
      
      if (split_score > best_score) {
        best_score <- split_score
        best_split <- list(
          feature = feature,
          threshold = threshold,
          importance = split_score,
          left_idx = left_idx,
          right_idx = right_idx
        )
      }
    }
  }
  
  return(best_split)
}

#' Calculate Split Score
#'
#' Calculates the quality of a split using appropriate metrics
#'
#' @param y Labels
#' @param left_idx Left split indices
#' @param right_idx Right split indices
#' @param is_classification Classification flag
#' @return Split score
calculate_split_score <- function(y, left_idx, right_idx, is_classification) {
  
  n_total <- length(y)
  n_left <- length(left_idx)
  n_right <- length(right_idx)
  
  if (is_classification) {
    # Use Gini impurity for classification
    gini_parent <- calculate_gini_impurity(y)
    gini_left <- calculate_gini_impurity(y[left_idx])
    gini_right <- calculate_gini_impurity(y[right_idx])
    
    weighted_gini <- (n_left / n_total) * gini_left + (n_right / n_total) * gini_right
    information_gain <- gini_parent - weighted_gini
    
    return(information_gain)
  } else {
    # Use variance reduction for regression
    var_parent <- var(y)
    var_left <- if(n_left > 1) var(y[left_idx]) else 0
    var_right <- if(n_right > 1) var(y[right_idx]) else 0
    
    weighted_var <- (n_left / n_total) * var_left + (n_right / n_total) * var_right
    variance_reduction <- var_parent - weighted_var
    
    return(variance_reduction)
  }
}

#' Calculate Gini Impurity
#'
#' Calculates Gini impurity for classification
#'
#' @param y Labels
#' @return Gini impurity
calculate_gini_impurity <- function(y) {
  
  if (length(y) == 0) return(0)
  
  class_counts <- table(y)
  class_probs <- class_counts / length(y)
  gini <- 1 - sum(class_probs^2)
  
  return(gini)
}

#' Get Majority Class
#'
#' Returns the majority class for classification
#'
#' @param y Labels
#' @return Majority class
get_majority_class <- function(y) {
  
  class_counts <- table(y)
  majority_class <- names(which.max(class_counts))
  
  return(majority_class)
}

#' Predict Tree
#'
#' Makes predictions using a single tree
#'
#' @param tree Decision tree
#' @param X Test data
#' @param is_classification Classification flag
#' @return Predictions
predict_tree <- function(tree, X, is_classification) {
  
  n_samples <- nrow(X)
  predictions <- character(n_samples)
  
  for (i in 1:n_samples) {
    predictions[i] <- predict_sample(tree$nodes[[1]], X[i, ])
  }
  
  if (!is_classification) {
    predictions <- as.numeric(predictions)
  }
  
  return(predictions)
}

#' Predict Sample
#'
#' Predicts a single sample using tree
#'
#' @param node Current tree node
#' @param sample Sample to predict
#' @return Prediction
predict_sample <- function(node, sample) {
  
  if (node$type == "leaf") {
    return(as.character(node$prediction))
  }
  
  if (sample[node$feature] <= node$threshold) {
    return(predict_sample(node$left, sample))
  } else {
    return(predict_sample(node$right, sample))
  }
}

#' Calculate OOB Error
#'
#' Calculates out-of-bag error rate
#'
#' @param oob_predictions OOB prediction matrix
#' @param oob_counts OOB count vector
#' @param y True labels
#' @param is_classification Classification flag
#' @return OOB error rate
calculate_oob_error <- function(oob_predictions, oob_counts, y, is_classification) {
  
  valid_idx <- which(oob_counts > 0)
  if (length(valid_idx) == 0) return(1.0)
  
  if (is_classification) {
    oob_pred_normalized <- oob_predictions[valid_idx, , drop = FALSE] / oob_counts[valid_idx]
    predicted_classes <- apply(oob_pred_normalized, 1, which.max)
    unique_classes <- unique(y)
    predicted_labels <- unique_classes[predicted_classes]
    
    error_rate <- mean(predicted_labels != y[valid_idx])
  } else {
    oob_pred_avg <- oob_predictions[valid_idx, 1] / oob_counts[valid_idx]
    error_rate <- mean((oob_pred_avg - y[valid_idx])^2)
  }
  
  return(error_rate)
}

#' Calculate RF Importance Measures
#'
#' Calculates multiple importance measures for Random Forest
#'
#' @param rf_ensemble Random Forest ensemble
#' @param X Input data
#' @param y Target variable
#' @param is_classification Classification flag
#' @return List of importance measures
calculate_rf_importance_measures <- function(rf_ensemble, X, y, is_classification) {
  
  # 1. Mean Decrease Impurity (from tree building)
  mdi_importance <- rf_ensemble$feature_importance
  
  # 2. Mean Decrease Accuracy (permutation importance)
  mda_importance <- calculate_permutation_importance(rf_ensemble, X, y, is_classification)
  
  # 3. Normalized importance measures
  mdi_normalized <- mdi_importance / sum(mdi_importance)
  mda_normalized <- mda_importance / sum(abs(mda_importance))
  
  return(list(
    mdi = mdi_importance,
    mda = mda_importance,
    mdi_normalized = mdi_normalized,
    mda_normalized = mda_normalized
  ))
}

#' Calculate Permutation Importance
#'
#' Calculates permutation-based feature importance
#'
#' @param rf_ensemble Random Forest ensemble
#' @param X Input data
#' @param y Target variable
#' @param is_classification Classification flag
#' @return Permutation importance scores
calculate_permutation_importance <- function(rf_ensemble, X, y, is_classification) {
  
  n_features <- ncol(X)
  
  # Get baseline predictions
  baseline_predictions <- predict_rf_ensemble(rf_ensemble, X, is_classification)
  baseline_error <- calculate_prediction_error(baseline_predictions, y, is_classification)
  
  permutation_importance <- numeric(n_features)
  
  for (feature in 1:n_features) {
    
    # Permute feature values
    X_permuted <- X
    X_permuted[, feature] <- sample(X_permuted[, feature])
    
    # Get predictions with permuted feature
    permuted_predictions <- predict_rf_ensemble(rf_ensemble, X_permuted, is_classification)
    permuted_error <- calculate_prediction_error(permuted_predictions, y, is_classification)
    
    # Importance is increase in error
    permutation_importance[feature] <- permuted_error - baseline_error
  }
  
  return(permutation_importance)
}

#' Predict RF Ensemble
#'
#' Makes predictions using the entire Random Forest ensemble
#'
#' @param rf_ensemble Random Forest ensemble
#' @param X Test data
#' @param is_classification Classification flag
#' @return Ensemble predictions
predict_rf_ensemble <- function(rf_ensemble, X, is_classification) {
  
  n_samples <- nrow(X)
  n_trees <- rf_ensemble$n_trees
  
  if (is_classification) {
    unique_classes <- unique(names(table(c())))  # Will be set properly in actual use
    n_classes <- length(unique_classes)
    vote_matrix <- matrix(0, nrow = n_samples, ncol = n_classes)
    
    for (tree_idx in 1:n_trees) {
      tree_predictions <- predict_tree(rf_ensemble$trees[[tree_idx]], X, is_classification)
      
      for (i in 1:n_samples) {
        class_idx <- which(unique_classes == tree_predictions[i])
        if (length(class_idx) > 0) {
          vote_matrix[i, class_idx] <- vote_matrix[i, class_idx] + 1
        }
      }
    }
    
    # Return majority vote
    predictions <- apply(vote_matrix, 1, which.max)
    return(unique_classes[predictions])
    
  } else {
    predictions <- numeric(n_samples)
    
    for (tree_idx in 1:n_trees) {
      tree_predictions <- predict_tree(rf_ensemble$trees[[tree_idx]], X, is_classification)
      predictions <- predictions + tree_predictions
    }
    
    return(predictions / n_trees)
  }
}

#' Calculate Prediction Error
#'
#' Calculates prediction error for importance calculation
#'
#' @param predictions Model predictions
#' @param y_true True labels
#' @param is_classification Classification flag
#' @return Prediction error
calculate_prediction_error <- function(predictions, y_true, is_classification) {
  
  if (is_classification) {
    return(mean(predictions != y_true))
  } else {
    return(mean((predictions - y_true)^2))
  }
}

#' Combine Importance Measures
#'
#' Combines multiple importance measures into final scores
#'
#' @param importance_measures List of importance measures
#' @return Combined importance scores
combine_importance_measures <- function(importance_measures) {
  
  # Combine MDI and MDA with equal weights
  mdi_norm <- importance_measures$mdi_normalized
  mda_norm <- abs(importance_measures$mda_normalized)
  
  # Handle edge cases
  if (sum(mdi_norm) == 0) mdi_norm <- rep(1, length(mdi_norm)) / length(mdi_norm)
  if (sum(mda_norm) == 0) mda_norm <- rep(1, length(mda_norm)) / length(mda_norm)
  
  combined_importance <- 0.6 * mdi_norm + 0.4 * mda_norm
  
  return(combined_importance)
}

#' Select Features by Importance
#'
#' Selects features based on importance scores
#'
#' @param importance_scores Feature importance scores
#' @param n_features Total number of features
#' @return Selected features and scores
select_features_by_importance <- function(importance_scores, n_features) {
  
  # Set selection threshold based on importance distribution
  threshold <- quantile(importance_scores, 0.7)  # Top 30% by default
  
  # Ensure minimum selection
  min_features <- max(1, round(0.1 * n_features))
  selected_idx <- which(importance_scores >= threshold)
  
  if (length(selected_idx) < min_features) {
    selected_idx <- order(importance_scores, decreasing = TRUE)[1:min_features]
    threshold <- importance_scores[selected_idx[min_features]]
  }
  
  # Ensure maximum selection
  max_features <- round(0.5 * n_features)
  if (length(selected_idx) > max_features) {
    selected_idx <- order(importance_scores, decreasing = TRUE)[1:max_features]
    threshold <- importance_scores[selected_idx[max_features]]
  }
  
  feature_names <- names(importance_scores)
  if (is.null(feature_names)) {
    feature_names <- paste0("Feature_", 1:n_features)
  }
  
  selected_features <- feature_names[selected_idx]
  feature_scores <- importance_scores[selected_idx]
  names(feature_scores) <- selected_features
  
  return(list(
    selected_features = selected_features,
    feature_scores = feature_scores,
    threshold = threshold
  ))
}

#' Determine Task Type
#'
#' Determines if the task is classification or regression
#'
#' @param target_vector Target variable
#' @return TRUE for classification, FALSE for regression
determine_task_type <- function(target_vector) {
  
  # Check if target is categorical or has few unique values
  unique_values <- length(unique(target_vector))
  total_values <- length(target_vector)
  
  # Classification if: factor, character, or few unique numeric values
  if (is.factor(target_vector) || is.character(target_vector)) {
    return(TRUE)
  }
  
  if (is.numeric(target_vector) && unique_values / total_values < 0.1 && unique_values <= 10) {
    return(TRUE)
  }
  
  return(FALSE)
}

#' Basic Feature Preprocessing
#'
#' Removes constant and highly correlated features
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @param cor_threshold Correlation threshold for removal
#' @return List with processed data and removal info
preprocess_features_basic <- function(data_matrix, target_vector, cor_threshold = 0.99) {
  
  # Remove constant features
  feature_vars <- apply(data_matrix, 2, var, na.rm = TRUE)
  constant_features <- which(feature_vars < 1e-10)
  
  if (length(constant_features) > 0) {
    data_matrix <- data_matrix[, -constant_features, drop = FALSE]
  }
  
  # Remove highly correlated features
  correlated_features <- c()
  if (ncol(data_matrix) > 1) {
    cor_matrix <- cor(data_matrix, use = "complete.obs")
    cor_matrix[is.na(cor_matrix)] <- 0
    
    for (i in 1:(ncol(cor_matrix) - 1)) {
      for (j in (i + 1):ncol(cor_matrix)) {
        if (abs(cor_matrix[i, j]) > cor_threshold) {
          correlated_features <- c(correlated_features, j)
        }
      }
    }
    
    if (length(correlated_features) > 0) {
      correlated_features <- unique(correlated_features)
      data_matrix <- data_matrix[, -correlated_features, drop = FALSE]
    }
  }
  
  return(list(
    data = data_matrix,
    constant_removed = length(constant_features),
    correlated_removed = length(correlated_features)
  ))
}

#' Create Empty Result for Random Forest
#'
#' Creates standardized empty result for failed feature selection
#'
#' @param message Error or warning message
#' @return Empty result list
create_empty_result <- function(message) {
  return(list(
    selected_features = character(0),
    feature_scores = numeric(0),
    method_name = "Random_Forest_Importance",
    execution_time = 0,
    n_features_selected = 0,
    warning_message = message
  ))
}

#' SMOTE Sampling Implementation
#'
#' Simple SMOTE implementation for handling class imbalance
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with augmented data and target
perform_smote_sampling <- function(data_matrix, target_vector) {
  
  # Get class frequencies
  class_counts <- table(target_vector)
  majority_class <- names(which.max(class_counts))
  majority_count <- max(class_counts)
  
  augmented_data <- data_matrix
  augmented_target <- target_vector
  
  # Generate synthetic samples for minority classes
  for (class_label in names(class_counts)) {
    if (class_label == majority_class) next
    
    class_indices <- which(target_vector == class_label)
    n_synthetic <- majority_count - class_counts[class_label]
    
    if (n_synthetic > 0 && length(class_indices) >= 2) {
      synthetic_samples <- generate_synthetic_samples(
        data_matrix[class_indices, , drop = FALSE], 
        n_synthetic
      )
      
      augmented_data <- rbind(augmented_data, synthetic_samples)
      augmented_target <- c(augmented_target, rep(class_label, n_synthetic))
    }
  }
  
  return(list(
    data = augmented_data,
    target = augmented_target
  ))
}

#' Generate Synthetic Samples for SMOTE
#'
#' Creates synthetic samples using k-nearest neighbors approach
#'
#' @param class_data Data for specific class
#' @param n_samples Number of synthetic samples to generate
#' @return Matrix of synthetic samples
generate_synthetic_samples <- function(class_data, n_samples) {
  
  n_class_samples <- nrow(class_data)
  n_features <- ncol(class_data)
  synthetic_samples <- matrix(0, nrow = n_samples, ncol = n_features)
  
  for (i in 1:n_samples) {
    # Randomly select a sample
    base_idx <- sample(n_class_samples, 1)
    base_sample <- class_data[base_idx, ]
    
    # Find nearest neighbor (simple approach)
    if (n_class_samples > 1) {
      neighbor_idx <- sample(setdiff(1:n_class_samples, base_idx), 1)
      neighbor_sample <- class_data[neighbor_idx, ]
      
      # Generate synthetic sample
      alpha <- runif(1)
      synthetic_samples[i, ] <- base_sample + alpha * (neighbor_sample - base_sample)
    } else {
      synthetic_samples[i, ] <- base_sample + rnorm(n_features, 0, 0.01)
    }
  }
  
  colnames(synthetic_samples) <- colnames(class_data)
  return(synthetic_samples)
}
