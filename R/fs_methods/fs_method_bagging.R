#' Bagging Feature Selection Method
#'
#' Implements feature selection using Bootstrap Aggregating (Bagging).
#' This method trains multiple models on bootstrap samples and aggregates
#' feature importance across the ensemble to achieve robust selection.
#'
#' @param ramwas.data Expression data matrix (samples as rows, features as columns)  
#' @param target Target variable vector
#' @param smote_sampling Logical, whether to apply SMOTE sampling
#' @return List containing selected features and method metadata
#' @export
bagging_fs <- function(ramwas.data, target, smote_sampling = FALSE) {
  
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
    
    # Perform Bagging feature selection
    bagging_result <- perform_bagging_selection(data_clean, target)
    
    # Extract selected features
    if (length(bagging_result$selected_features) == 0) {
      return(create_empty_result("No features selected by Bagging"))
    }
    
    # Calculate timing
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      selected_features = bagging_result$selected_features,
      feature_scores = bagging_result$feature_scores,
      method_name = "Bagging",
      execution_time = execution_time,
      n_features_selected = length(bagging_result$selected_features),
      n_features_input = ncol(ramwas.data),
      method_params = list(
        n_bags = bagging_result$n_bags,
        bag_size_fraction = bagging_result$bag_size_fraction,
        base_learner = bagging_result$base_learner,
        feature_subsample = bagging_result$feature_subsample,
        aggregation_method = bagging_result$aggregation_method,
        smote_sampling = smote_sampling
      ),
      preprocessing_info = list(
        constant_features_removed = preprocessed$constant_removed,
        correlated_features_removed = preprocessed$correlated_removed,
        final_features_count = ncol(data_clean)
      ),
      bagging_info = bagging_result$bagging_info
    ))
    
  }, error = function(e) {
    # Return error result
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      selected_features = character(0),
      feature_scores = numeric(0),
      method_name = "Bagging",
      execution_time = execution_time,
      n_features_selected = 0,
      n_features_input = ncol(ramwas.data),
      error_message = paste("Bagging failed:", e$message),
      method_params = list(smote_sampling = smote_sampling)
    ))
  })
}

#' Perform Bagging Feature Selection
#'
#' Implements Bootstrap Aggregating with feature importance aggregation
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with selected features and bagging information
perform_bagging_selection <- function(data_matrix, target_vector) {
  
  n_samples <- nrow(data_matrix)
  n_features <- ncol(data_matrix)
  
  # Determine task type
  is_classification <- determine_task_type(target_vector)
  
  # Set bagging parameters
  n_bags <- max(20, min(100, round(sqrt(n_samples))))
  bag_size_fraction <- 0.8
  feature_subsample <- min(1.0, max(0.5, sqrt(n_features) / n_features))
  
  # Standardize features
  data_scaled <- scale(data_matrix)
  
  # Handle any NAs from scaling
  if (any(is.na(data_scaled))) {
    non_na_cols <- apply(data_scaled, 2, function(x) !any(is.na(x)))
    data_scaled <- data_scaled[, non_na_cols, drop = FALSE]
    n_features <- ncol(data_scaled)
  }
  
  # Train bagging ensemble
  bagging_ensemble <- train_bagging_ensemble(data_scaled, target_vector, n_bags, 
                                           bag_size_fraction, feature_subsample, is_classification)
  
  # Aggregate feature importance
  aggregated_importance <- aggregate_bag_importance(bagging_ensemble, n_features)
  
  # Select features based on aggregated importance
  selected_result <- select_features_by_bagging_importance(aggregated_importance, n_features)
  
  return(list(
    selected_features = selected_result$selected_features,
    feature_scores = selected_result$feature_scores,
    n_bags = n_bags,
    bag_size_fraction = bag_size_fraction,
    base_learner = "decision_tree",
    feature_subsample = feature_subsample,
    aggregation_method = "weighted_average",
    bagging_info = list(
      task_type = if(is_classification) "classification" else "regression",
      oob_error = bagging_ensemble$oob_error,
      bag_performances = bagging_ensemble$bag_performances,
      feature_selection_frequency = bagging_ensemble$feature_selection_frequency
    )
  ))
}

#' Train Bagging Ensemble
#'
#' Trains ensemble of models using bootstrap sampling
#'
#' @param X Input data matrix
#' @param y Target vector
#' @param n_bags Number of bootstrap bags
#' @param bag_size_fraction Fraction of samples in each bag
#' @param feature_subsample Fraction of features to subsample
#' @param is_classification Classification flag
#' @return Trained bagging ensemble
train_bagging_ensemble <- function(X, y, n_bags, bag_size_fraction, feature_subsample, is_classification) {
  
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  # Initialize storage
  bag_models <- list()
  bag_importance <- matrix(0, nrow = n_bags, ncol = n_features)
  bag_performances <- numeric(n_bags)
  feature_selection_frequency <- integer(n_features)
  
  # Out-of-bag prediction tracking
  oob_predictions <- numeric(n_samples)
  oob_counts <- integer(n_samples)
  
  # Train each bag
  for (bag in 1:n_bags) {
    
    # Bootstrap sampling
    bag_size <- round(bag_size_fraction * n_samples)
    bag_idx <- sample(n_samples, bag_size, replace = TRUE)
    oob_idx <- setdiff(1:n_samples, unique(bag_idx))
    
    # Feature subsampling
    n_features_use <- max(1, round(feature_subsample * n_features))
    feature_idx <- sort(sample(n_features, n_features_use))
    
    # Extract bag data
    X_bag <- X[bag_idx, feature_idx, drop = FALSE]
    y_bag <- y[bag_idx]
    
    # Train base learner
    base_model <- train_bagging_base_learner(X_bag, y_bag, is_classification)
    
    # Calculate feature importance for this bag
    bag_importance_raw <- calculate_base_learner_importance(base_model, X_bag, y_bag, is_classification)
    
    # Map back to full feature space
    bag_importance[bag, feature_idx] <- bag_importance_raw
    feature_selection_frequency[feature_idx] <- feature_selection_frequency[feature_idx] + 1
    
    # Store model
    bag_models[[bag]] <- list(
      model = base_model,
      feature_indices = feature_idx,
      oob_indices = oob_idx
    )
    
    # Out-of-bag evaluation
    if (length(oob_idx) > 0) {
      oob_pred <- predict_bagging_base_learner(base_model, X[oob_idx, feature_idx, drop = FALSE], is_classification)
      
      if (is_classification) {
        # For classification, use majority voting
        for (i in seq_along(oob_idx)) {
          if (oob_counts[oob_idx[i]] == 0) {
            oob_predictions[oob_idx[i]] <- as.numeric(oob_pred[i] == unique(y)[1])
          } else {
            # Simple averaging for binary classification
            current_pred <- oob_predictions[oob_idx[i]] * oob_counts[oob_idx[i]]
            new_pred <- as.numeric(oob_pred[i] == unique(y)[1])
            oob_predictions[oob_idx[i]] <- (current_pred + new_pred) / (oob_counts[oob_idx[i]] + 1)
          }
          oob_counts[oob_idx[i]] <- oob_counts[oob_idx[i]] + 1
        }
      } else {
        # For regression, average predictions
        for (i in seq_along(oob_idx)) {
          if (oob_counts[oob_idx[i]] == 0) {
            oob_predictions[oob_idx[i]] <- oob_pred[i]
          } else {
            current_pred <- oob_predictions[oob_idx[i]] * oob_counts[oob_idx[i]]
            oob_predictions[oob_idx[i]] <- (current_pred + oob_pred[i]) / (oob_counts[oob_idx[i]] + 1)
          }
          oob_counts[oob_idx[i]] <- oob_counts[oob_idx[i]] + 1
        }
      }
      
      # Calculate bag performance on its OOB data
      if (is_classification) {
        oob_pred_binary <- ifelse(oob_pred == unique(y)[1], 0, 1)
        y_oob_binary <- ifelse(y[oob_idx] == unique(y)[1], 0, 1)
        bag_performances[bag] <- mean(oob_pred_binary == y_oob_binary)
      } else {
        bag_performances[bag] <- 1 / (1 + mean((oob_pred - y[oob_idx])^2))  # Convert MSE to performance score
      }
    } else {
      bag_performances[bag] <- 0.5  # Default performance
    }
  }
  
  # Calculate overall OOB error
  valid_oob_idx <- which(oob_counts > 0)
  if (length(valid_oob_idx) > 0) {
    if (is_classification) {
      final_oob_pred <- ifelse(oob_predictions[valid_oob_idx] > 0.5, unique(y)[2], unique(y)[1])
      oob_error <- mean(final_oob_pred != y[valid_oob_idx])
    } else {
      oob_error <- mean((oob_predictions[valid_oob_idx] - y[valid_oob_idx])^2)
    }
  } else {
    oob_error <- 1.0
  }
  
  return(list(
    bag_models = bag_models,
    bag_importance = bag_importance,
    bag_performances = bag_performances,
    oob_error = oob_error,
    feature_selection_frequency = feature_selection_frequency
  ))
}

#' Train Bagging Base Learner
#'
#' Trains a single base learner for bagging (decision tree)
#'
#' @param X Training data
#' @param y Training labels
#' @param is_classification Classification flag
#' @return Trained base learner
train_bagging_base_learner <- function(X, y, is_classification) {
  
  if (is_classification) {
    # Use classification tree
    tree <- build_bagging_classification_tree(X, y, max_depth = 6)
  } else {
    # Use regression tree
    tree <- build_bagging_regression_tree(X, y, max_depth = 6)
  }
  
  return(list(
    tree = tree,
    is_classification = is_classification
  ))
}

#' Build Bagging Classification Tree
#'
#' Builds classification tree for bagging
#'
#' @param X Training data
#' @param y Training labels
#' @param max_depth Maximum tree depth
#' @return Classification tree
build_bagging_classification_tree <- function(X, y, max_depth = 6) {
  
  tree <- build_bagging_classification_node(X, y, 1, max_depth)
  return(tree)
}

#' Build Bagging Classification Node
#'
#' Recursively builds classification tree nodes for bagging
#'
#' @param X Node data
#' @param y Node labels
#' @param depth Current depth
#' @param max_depth Maximum depth
#' @return Tree node
build_bagging_classification_node <- function(X, y, depth, max_depth) {
  
  n_samples <- nrow(X)
  
  # Stopping criteria
  if (depth >= max_depth || n_samples < 3 || length(unique(y)) == 1) {
    return(list(
      type = "leaf",
      prediction = get_majority_class(y),
      n_samples = n_samples,
      class_distribution = table(y)
    ))
  }
  
  # Find best split
  best_split <- find_best_bagging_classification_split(X, y)
  
  if (is.null(best_split)) {
    return(list(
      type = "leaf",
      prediction = get_majority_class(y),
      n_samples = n_samples,
      class_distribution = table(y)
    ))
  }
  
  # Split data
  left_idx <- which(X[, best_split$feature] <= best_split$threshold)
  right_idx <- which(X[, best_split$feature] > best_split$threshold)
  
  # Build child nodes
  left_node <- build_bagging_classification_node(X[left_idx, , drop = FALSE], y[left_idx], depth + 1, max_depth)
  right_node <- build_bagging_classification_node(X[right_idx, , drop = FALSE], y[right_idx], depth + 1, max_depth)
  
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

#' Find Best Bagging Classification Split
#'
#' Finds the best split for bagging classification tree
#'
#' @param X Node data
#' @param y Node labels
#' @return Best split information
find_best_bagging_classification_split <- function(X, y) {
  
  n_features <- ncol(X)
  best_score <- -Inf
  best_split <- NULL
  
  # Consider only a subset of features for additional randomness
  features_to_try <- sample(n_features, max(1, round(sqrt(n_features))))
  
  for (feature in features_to_try) {
    
    feature_values <- X[, feature]
    unique_values <- sort(unique(feature_values))
    
    if (length(unique_values) < 2) next
    
    # Try random thresholds for additional randomness
    n_thresholds <- min(10, length(unique_values) - 1)
    threshold_indices <- sample(1:(length(unique_values) - 1), n_thresholds)
    
    for (i in threshold_indices) {
      threshold <- (unique_values[i] + unique_values[i + 1]) / 2
      
      left_idx <- which(feature_values <= threshold)
      right_idx <- which(feature_values > threshold)
      
      if (length(left_idx) == 0 || length(right_idx) == 0) next
      
      # Calculate Gini impurity reduction
      parent_gini <- calculate_gini_impurity(y)
      left_gini <- calculate_gini_impurity(y[left_idx])
      right_gini <- calculate_gini_impurity(y[right_idx])
      
      weighted_gini <- (length(left_idx) * left_gini + length(right_idx) * right_gini) / length(y)
      gini_reduction <- parent_gini - weighted_gini
      
      if (gini_reduction > best_score) {
        best_score <- gini_reduction
        best_split <- list(
          feature = feature,
          threshold = threshold,
          importance = gini_reduction,
          left_idx = left_idx,
          right_idx = right_idx
        )
      }
    }
  }
  
  return(best_split)
}

#' Build Bagging Regression Tree
#'
#' Builds regression tree for bagging
#'
#' @param X Training data
#' @param y Training targets
#' @param max_depth Maximum tree depth
#' @return Regression tree
build_bagging_regression_tree <- function(X, y, max_depth = 6) {
  
  tree <- build_bagging_regression_node(X, y, 1, max_depth)
  return(tree)
}

#' Build Bagging Regression Node
#'
#' Recursively builds regression tree nodes for bagging
#'
#' @param X Node data
#' @param y Node targets
#' @param depth Current depth
#' @param max_depth Maximum depth
#' @return Tree node
build_bagging_regression_node <- function(X, y, depth, max_depth) {
  
  n_samples <- nrow(X)
  
  # Stopping criteria
  if (depth >= max_depth || n_samples < 3 || var(y) < 1e-6) {
    return(list(
      type = "leaf",
      prediction = mean(y),
      n_samples = n_samples,
      variance = var(y)
    ))
  }
  
  # Find best split
  best_split <- find_best_bagging_regression_split(X, y)
  
  if (is.null(best_split)) {
    return(list(
      type = "leaf",
      prediction = mean(y),
      n_samples = n_samples,
      variance = var(y)
    ))
  }
  
  # Split data
  left_idx <- which(X[, best_split$feature] <= best_split$threshold)
  right_idx <- which(X[, best_split$feature] > best_split$threshold)
  
  # Build child nodes
  left_node <- build_bagging_regression_node(X[left_idx, , drop = FALSE], y[left_idx], depth + 1, max_depth)
  right_node <- build_bagging_regression_node(X[right_idx, , drop = FALSE], y[right_idx], depth + 1, max_depth)
  
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

#' Find Best Bagging Regression Split
#'
#' Finds the best split for bagging regression tree
#'
#' @param X Node data
#' @param y Node targets
#' @return Best split information
find_best_bagging_regression_split <- function(X, y) {
  
  n_features <- ncol(X)
  best_score <- -Inf
  best_split <- NULL
  
  # Consider only a subset of features for additional randomness
  features_to_try <- sample(n_features, max(1, round(sqrt(n_features))))
  
  for (feature in features_to_try) {
    
    feature_values <- X[, feature]
    unique_values <- sort(unique(feature_values))
    
    if (length(unique_values) < 2) next
    
    # Try random thresholds for additional randomness
    n_thresholds <- min(10, length(unique_values) - 1)
    threshold_indices <- sample(1:(length(unique_values) - 1), n_thresholds)
    
    for (i in threshold_indices) {
      threshold <- (unique_values[i] + unique_values[i + 1]) / 2
      
      left_idx <- which(feature_values <= threshold)
      right_idx <- which(feature_values > threshold)
      
      if (length(left_idx) == 0 || length(right_idx) == 0) next
      
      # Calculate variance reduction
      total_var <- var(y)
      left_var <- if(length(left_idx) > 1) var(y[left_idx]) else 0
      right_var <- if(length(right_idx) > 1) var(y[right_idx]) else 0
      
      weighted_var <- (length(left_idx) * left_var + length(right_idx) * right_var) / length(y)
      variance_reduction <- total_var - weighted_var
      
      if (variance_reduction > best_score) {
        best_score <- variance_reduction
        best_split <- list(
          feature = feature,
          threshold = threshold,
          importance = variance_reduction,
          left_idx = left_idx,
          right_idx = right_idx
        )
      }
    }
  }
  
  return(best_split)
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
#' Returns the majority class
#'
#' @param y Labels
#' @return Majority class
get_majority_class <- function(y) {
  
  class_counts <- table(y)
  majority_class <- names(which.max(class_counts))
  
  return(majority_class)
}

#' Predict Bagging Base Learner
#'
#' Makes predictions using a bagging base learner
#'
#' @param model Trained base learner
#' @param X Test data
#' @param is_classification Classification flag
#' @return Predictions
predict_bagging_base_learner <- function(model, X, is_classification) {
  
  n_samples <- nrow(X)
  predictions <- if(is_classification) character(n_samples) else numeric(n_samples)
  
  for (i in 1:n_samples) {
    predictions[i] <- predict_bagging_tree_sample(model$tree, X[i, ])
  }
  
  if (!is_classification) {
    predictions <- as.numeric(predictions)
  }
  
  return(predictions)
}

#' Predict Bagging Tree Sample
#'
#' Predicts a single sample using bagging tree
#'
#' @param node Current tree node
#' @param sample Sample to predict
#' @return Prediction
predict_bagging_tree_sample <- function(node, sample) {
  
  if (node$type == "leaf") {
    return(as.character(node$prediction))
  }
  
  if (sample[node$feature] <= node$threshold) {
    return(predict_bagging_tree_sample(node$left, sample))
  } else {
    return(predict_bagging_tree_sample(node$right, sample))
  }
}

#' Calculate Base Learner Importance
#'
#' Calculates feature importance for a single base learner
#'
#' @param model Trained base learner
#' @param X Training data
#' @param y Training targets
#' @param is_classification Classification flag
#' @return Feature importance scores
calculate_base_learner_importance <- function(model, X, y, is_classification) {
  
  n_features <- ncol(X)
  importance <- calculate_bagging_tree_importance(model$tree, n_features)
  
  # Normalize importance
  if (sum(importance) > 0) {
    importance <- importance / sum(importance)
  } else {
    importance <- rep(1/n_features, n_features)
  }
  
  return(importance)
}

#' Calculate Bagging Tree Importance
#'
#' Calculates feature importance for bagging tree
#'
#' @param tree Bagging tree
#' @param n_features Number of features
#' @return Feature importance vector
calculate_bagging_tree_importance <- function(tree, n_features) {
  
  importance <- numeric(n_features)
  
  if (tree$type == "internal") {
    importance[tree$feature] <- tree$importance
    
    # Recursively add importance from child nodes
    left_importance <- calculate_bagging_tree_importance(tree$left, n_features)
    right_importance <- calculate_bagging_tree_importance(tree$right, n_features)
    
    importance <- importance + left_importance + right_importance
  }
  
  return(importance)
}

#' Aggregate Bag Importance
#'
#' Aggregates feature importance across all bags
#'
#' @param bagging_ensemble Bagging ensemble results
#' @param n_features Number of features
#' @return Aggregated feature importance
aggregate_bag_importance <- function(bagging_ensemble, n_features) {
  
  bag_importance <- bagging_ensemble$bag_importance
  bag_performances <- bagging_ensemble$bag_performances
  feature_selection_frequency <- bagging_ensemble$feature_selection_frequency
  
  # Weight importance by bag performance
  performance_weights <- bag_performances / sum(bag_performances)
  
  # Weighted average of bag importances
  weighted_importance <- colSums(bag_importance * performance_weights)
  
  # Adjust by selection frequency (features selected more often get bonus)
  frequency_bonus <- feature_selection_frequency / max(feature_selection_frequency)
  
  # Combine weighted importance with frequency bonus
  final_importance <- 0.8 * weighted_importance + 0.2 * frequency_bonus
  
  # Normalize
  if (sum(final_importance) > 0) {
    final_importance <- final_importance / sum(final_importance)
  } else {
    final_importance <- rep(1/n_features, n_features)
  }
  
  return(final_importance)
}

#' Select Features by Bagging Importance
#'
#' Selects features based on bagging importance scores
#'
#' @param importance_scores Feature importance scores
#' @param n_features Total number of features
#' @return Selected features and scores
select_features_by_bagging_importance <- function(importance_scores, n_features) {
  
  # Set selection threshold based on importance distribution
  threshold <- quantile(importance_scores, 0.75)  # Top 25% by default
  
  # Ensure minimum selection
  min_features <- max(1, round(0.05 * n_features))
  selected_idx <- which(importance_scores >= threshold)
  
  if (length(selected_idx) < min_features) {
    selected_idx <- order(importance_scores, decreasing = TRUE)[1:min_features]
    threshold <- importance_scores[selected_idx[min_features]]
  }
  
  # Ensure maximum selection
  max_features <- round(0.4 * n_features)
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

#' Create Empty Result for Bagging
#'
#' Creates standardized empty result for failed feature selection
#'
#' @param message Error or warning message
#' @return Empty result list
create_empty_result <- function(message) {
  return(list(
    selected_features = character(0),
    feature_scores = numeric(0),
    method_name = "Bagging",
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
