#' Extra Trees Feature Selection Method
#'
#' Implements feature selection using Extremely Randomized Trees (Extra Trees).
#' This method builds an ensemble of unpruned decision trees with random
#' threshold selection and provides feature importance through averaging.
#'
#' @param ramwas.data Expression data matrix (samples as rows, features as columns)  
#' @param target Target variable vector
#' @param smote_sampling Logical, whether to apply SMOTE sampling
#' @return List containing selected features and method metadata
#' @export
extra_trees_fs <- function(ramwas.data, target, smote_sampling = FALSE) {
  
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
    
    # Perform Extra Trees feature selection
    extra_trees_result <- perform_extra_trees_selection(data_clean, target)
    
    # Extract selected features
    if (length(extra_trees_result$selected_features) == 0) {
      return(create_empty_result("No features selected by Extra Trees"))
    }
    
    # Calculate timing
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      selected_features = extra_trees_result$selected_features,
      feature_scores = extra_trees_result$feature_scores,
      method_name = "Extra_Trees",
      execution_time = execution_time,
      n_features_selected = length(extra_trees_result$selected_features),
      n_features_input = ncol(ramwas.data),
      method_params = list(
        n_trees = extra_trees_result$n_trees,
        max_features = extra_trees_result$max_features,
        min_samples_split = extra_trees_result$min_samples_split,
        min_samples_leaf = extra_trees_result$min_samples_leaf,
        bootstrap = extra_trees_result$bootstrap,
        smote_sampling = smote_sampling
      ),
      preprocessing_info = list(
        constant_features_removed = preprocessed$constant_removed,
        correlated_features_removed = preprocessed$correlated_removed,
        final_features_count = ncol(data_clean)
      ),
      extra_trees_info = extra_trees_result$extra_trees_info
    ))
    
  }, error = function(e) {
    # Return error result
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      selected_features = character(0),
      feature_scores = numeric(0),
      method_name = "Extra_Trees",
      execution_time = execution_time,
      n_features_selected = 0,
      n_features_input = ncol(ramwas.data),
      error_message = paste("Extra Trees failed:", e$message),
      method_params = list(smote_sampling = smote_sampling)
    ))
  })
}

#' Perform Extra Trees Feature Selection
#'
#' Implements Extremely Randomized Trees with feature importance calculation
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with selected features and Extra Trees information
perform_extra_trees_selection <- function(data_matrix, target_vector) {
  
  n_samples <- nrow(data_matrix)
  n_features <- ncol(data_matrix)
  
  # Determine task type
  is_classification <- determine_task_type(target_vector)
  
  # Set Extra Trees parameters
  n_trees <- max(50, min(200, round(sqrt(n_features) * 8)))
  max_features <- max(1, round(sqrt(n_features)))
  min_samples_split <- max(2, round(0.01 * n_samples))
  min_samples_leaf <- max(1, round(0.005 * n_samples))
  bootstrap <- FALSE  # Extra Trees typically use full dataset
  
  # Standardize features
  data_scaled <- scale(data_matrix)
  
  # Handle any NAs from scaling
  if (any(is.na(data_scaled))) {
    non_na_cols <- apply(data_scaled, 2, function(x) !any(is.na(x)))
    data_scaled <- data_scaled[, non_na_cols, drop = FALSE]
    n_features <- ncol(data_scaled)
    max_features <- max(1, round(sqrt(n_features)))
  }
  
  # Train Extra Trees ensemble
  extra_trees_ensemble <- train_extra_trees_ensemble(data_scaled, target_vector, n_trees, 
                                                   max_features, min_samples_split, 
                                                   min_samples_leaf, bootstrap, is_classification)
  
  # Calculate feature importance
  feature_importance <- calculate_extra_trees_importance(extra_trees_ensemble, data_scaled, target_vector, is_classification)
  
  # Select features based on importance
  selected_result <- select_features_by_extra_trees_importance(feature_importance, n_features)
  
  return(list(
    selected_features = selected_result$selected_features,
    feature_scores = selected_result$feature_scores,
    n_trees = n_trees,
    max_features = max_features,
    min_samples_split = min_samples_split,
    min_samples_leaf = min_samples_leaf,
    bootstrap = bootstrap,
    extra_trees_info = list(
      task_type = if(is_classification) "classification" else "regression",
      oob_score = extra_trees_ensemble$oob_score,
      feature_usage_frequency = extra_trees_ensemble$feature_usage_frequency,
      tree_depths = extra_trees_ensemble$tree_depths
    )
  ))
}

#' Train Extra Trees Ensemble
#'
#' Trains ensemble of extremely randomized trees
#'
#' @param X Input data matrix
#' @param y Target vector
#' @param n_trees Number of trees
#' @param max_features Maximum features per split
#' @param min_samples_split Minimum samples for split
#' @param min_samples_leaf Minimum samples per leaf
#' @param bootstrap Whether to use bootstrap sampling
#' @param is_classification Classification flag
#' @return Trained Extra Trees ensemble
train_extra_trees_ensemble <- function(X, y, n_trees, max_features, min_samples_split, 
                                     min_samples_leaf, bootstrap, is_classification) {
  
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  # Initialize storage
  trees <- list()
  feature_importance <- matrix(0, nrow = n_trees, ncol = n_features)
  feature_usage_frequency <- integer(n_features)
  tree_depths <- integer(n_trees)
  
  # Out-of-bag tracking (if bootstrap is used)
  if (bootstrap) {
    oob_predictions <- matrix(0, nrow = n_samples, ncol = if(is_classification) length(unique(y)) else 1)
    oob_counts <- integer(n_samples)
  } else {
    oob_predictions <- NULL
    oob_counts <- NULL
  }
  
  # Build trees
  for (tree_idx in 1:n_trees) {
    
    # Sample selection
    if (bootstrap) {
      sample_idx <- sample(n_samples, n_samples, replace = TRUE)
      oob_idx <- setdiff(1:n_samples, unique(sample_idx))
    } else {
      sample_idx <- 1:n_samples
      oob_idx <- c()
    }
    
    X_tree <- X[sample_idx, , drop = FALSE]
    y_tree <- y[sample_idx]
    
    # Build extremely randomized tree
    tree <- build_extra_tree(X_tree, y_tree, max_features, min_samples_split, 
                           min_samples_leaf, is_classification)
    
    trees[[tree_idx]] <- tree
    
    # Calculate feature importance for this tree
    tree_importance <- calculate_tree_feature_importance_extra(tree, n_features)
    feature_importance[tree_idx, ] <- tree_importance
    
    # Track feature usage
    used_features <- get_tree_used_features(tree)
    feature_usage_frequency[used_features] <- feature_usage_frequency[used_features] + 1
    
    # Track tree depth
    tree_depths[tree_idx] <- calculate_tree_depth(tree)
    
    # Out-of-bag predictions
    if (bootstrap && length(oob_idx) > 0) {
      oob_pred <- predict_extra_tree(tree, X[oob_idx, , drop = FALSE], is_classification)
      
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
  }
  
  # Calculate OOB score
  if (bootstrap && !is.null(oob_counts)) {
    oob_score <- calculate_extra_trees_oob_score(oob_predictions, oob_counts, y, is_classification)
  } else {
    oob_score <- NULL
  }
  
  return(list(
    trees = trees,
    feature_importance = feature_importance,
    feature_usage_frequency = feature_usage_frequency,
    tree_depths = tree_depths,
    oob_score = oob_score,
    n_trees = n_trees
  ))
}

#' Build Extra Tree
#'
#' Builds a single extremely randomized tree
#'
#' @param X Training data
#' @param y Training labels
#' @param max_features Maximum features per split
#' @param min_samples_split Minimum samples for split
#' @param min_samples_leaf Minimum samples per leaf
#' @param is_classification Classification flag
#' @return Extra tree
build_extra_tree <- function(X, y, max_features, min_samples_split, min_samples_leaf, is_classification) {
  
  tree <- build_extra_tree_node(X, y, max_features, min_samples_split, min_samples_leaf, 
                               is_classification, depth = 1, max_depth = 20)
  return(tree)
}

#' Build Extra Tree Node
#'
#' Recursively builds extremely randomized tree nodes
#'
#' @param X Node data
#' @param y Node labels
#' @param max_features Maximum features per split
#' @param min_samples_split Minimum samples for split
#' @param min_samples_leaf Minimum samples per leaf
#' @param is_classification Classification flag
#' @param depth Current depth
#' @param max_depth Maximum depth
#' @return Tree node
build_extra_tree_node <- function(X, y, max_features, min_samples_split, min_samples_leaf, 
                                is_classification, depth, max_depth) {
  
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  # Stopping criteria
  if (depth >= max_depth || 
      n_samples < min_samples_split || 
      n_samples < 2 * min_samples_leaf ||
      length(unique(y)) == 1) {
    
    return(list(
      type = "leaf",
      prediction = if(is_classification) get_majority_class(y) else mean(y),
      n_samples = n_samples,
      impurity = if(is_classification) calculate_gini_impurity(y) else var(y)
    ))
  }
  
  # Find extremely randomized split
  best_split <- find_extra_random_split(X, y, max_features, min_samples_leaf, is_classification)
  
  if (is.null(best_split)) {
    return(list(
      type = "leaf",
      prediction = if(is_classification) get_majority_class(y) else mean(y),
      n_samples = n_samples,
      impurity = if(is_classification) calculate_gini_impurity(y) else var(y)
    ))
  }
  
  # Split data
  left_idx <- which(X[, best_split$feature] <= best_split$threshold)
  right_idx <- which(X[, best_split$feature] > best_split$threshold)
  
  # Build child nodes
  left_node <- build_extra_tree_node(X[left_idx, , drop = FALSE], y[left_idx], 
                                   max_features, min_samples_split, min_samples_leaf,
                                   is_classification, depth + 1, max_depth)
  
  right_node <- build_extra_tree_node(X[right_idx, , drop = FALSE], y[right_idx],
                                    max_features, min_samples_split, min_samples_leaf,
                                    is_classification, depth + 1, max_depth)
  
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

#' Find Extra Random Split
#'
#' Finds extremely randomized split (random threshold, random features)
#'
#' @param X Node data
#' @param y Node labels
#' @param max_features Maximum features to consider
#' @param min_samples_leaf Minimum samples per leaf
#' @param is_classification Classification flag
#' @return Best random split
find_extra_random_split <- function(X, y, max_features, min_samples_leaf, is_classification) {
  
  # Randomly select features to consider
  features_to_try <- sample(ncol(X), min(max_features, ncol(X)))
  
  best_score <- -Inf
  best_split <- NULL
  
  for (feature in features_to_try) {
    
    feature_values <- X[, feature]
    min_val <- min(feature_values)
    max_val <- max(feature_values)
    
    # Skip if all values are the same
    if (min_val == max_val) next
    
    # Generate random threshold
    threshold <- runif(1, min_val, max_val)
    
    left_idx <- which(feature_values <= threshold)
    right_idx <- which(feature_values > threshold)
    
    # Check minimum samples constraint
    if (length(left_idx) < min_samples_leaf || length(right_idx) < min_samples_leaf) next
    
    # Calculate split score
    if (is_classification) {
      # Gini impurity reduction
      parent_gini <- calculate_gini_impurity(y)
      left_gini <- calculate_gini_impurity(y[left_idx])
      right_gini <- calculate_gini_impurity(y[right_idx])
      
      weighted_gini <- (length(left_idx) * left_gini + length(right_idx) * right_gini) / length(y)
      score <- parent_gini - weighted_gini
    } else {
      # Variance reduction
      parent_var <- var(y)
      left_var <- if(length(left_idx) > 1) var(y[left_idx]) else 0
      right_var <- if(length(right_idx) > 1) var(y[right_idx]) else 0
      
      weighted_var <- (length(left_idx) * left_var + length(right_idx) * right_var) / length(y)
      score <- parent_var - weighted_var
    }
    
    if (score > best_score) {
      best_score <- score
      best_split <- list(
        feature = feature,
        threshold = threshold,
        importance = score,
        left_idx = left_idx,
        right_idx = right_idx
      )
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

#' Predict Extra Tree
#'
#' Makes predictions using an extra tree
#'
#' @param tree Extra tree
#' @param X Test data
#' @param is_classification Classification flag
#' @return Predictions
predict_extra_tree <- function(tree, X, is_classification) {
  
  n_samples <- nrow(X)
  predictions <- if(is_classification) character(n_samples) else numeric(n_samples)
  
  for (i in 1:n_samples) {
    predictions[i] <- predict_extra_tree_sample(tree, X[i, ])
  }
  
  if (!is_classification) {
    predictions <- as.numeric(predictions)
  }
  
  return(predictions)
}

#' Predict Extra Tree Sample
#'
#' Predicts a single sample using extra tree
#'
#' @param node Current tree node
#' @param sample Sample to predict
#' @return Prediction
predict_extra_tree_sample <- function(node, sample) {
  
  if (node$type == "leaf") {
    return(as.character(node$prediction))
  }
  
  if (sample[node$feature] <= node$threshold) {
    return(predict_extra_tree_sample(node$left, sample))
  } else {
    return(predict_extra_tree_sample(node$right, sample))
  }
}

#' Calculate Tree Feature Importance Extra
#'
#' Calculates feature importance for an extra tree
#'
#' @param tree Extra tree
#' @param n_features Number of features
#' @return Feature importance vector
calculate_tree_feature_importance_extra <- function(tree, n_features) {
  
  importance <- numeric(n_features)
  
  if (tree$type == "internal") {
    # Weight by number of samples affected
    importance[tree$feature] <- tree$importance * tree$n_samples
    
    # Recursively add importance from child nodes
    left_importance <- calculate_tree_feature_importance_extra(tree$left, n_features)
    right_importance <- calculate_tree_feature_importance_extra(tree$right, n_features)
    
    importance <- importance + left_importance + right_importance
  }
  
  return(importance)
}

#' Get Tree Used Features
#'
#' Gets all features used in a tree
#'
#' @param tree Tree object
#' @return Vector of used feature indices
get_tree_used_features <- function(tree) {
  
  used_features <- c()
  
  if (tree$type == "internal") {
    used_features <- c(used_features, tree$feature)
    
    # Recursively get features from child nodes
    left_features <- get_tree_used_features(tree$left)
    right_features <- get_tree_used_features(tree$right)
    
    used_features <- c(used_features, left_features, right_features)
  }
  
  return(unique(used_features))
}

#' Calculate Tree Depth
#'
#' Calculates the depth of a tree
#'
#' @param tree Tree object
#' @return Tree depth
calculate_tree_depth <- function(tree) {
  
  if (tree$type == "leaf") {
    return(1)
  }
  
  left_depth <- calculate_tree_depth(tree$left)
  right_depth <- calculate_tree_depth(tree$right)
  
  return(1 + max(left_depth, right_depth))
}

#' Calculate Extra Trees OOB Score
#'
#' Calculates out-of-bag score for Extra Trees
#'
#' @param oob_predictions OOB prediction matrix
#' @param oob_counts OOB count vector
#' @param y True labels
#' @param is_classification Classification flag
#' @return OOB score
calculate_extra_trees_oob_score <- function(oob_predictions, oob_counts, y, is_classification) {
  
  valid_idx <- which(oob_counts > 0)
  if (length(valid_idx) == 0) return(0.0)
  
  if (is_classification) {
    oob_pred_normalized <- oob_predictions[valid_idx, , drop = FALSE] / oob_counts[valid_idx]
    predicted_classes <- apply(oob_pred_normalized, 1, which.max)
    unique_classes <- unique(y)
    predicted_labels <- unique_classes[predicted_classes]
    
    accuracy <- mean(predicted_labels == y[valid_idx])
    return(accuracy)
  } else {
    oob_pred_avg <- oob_predictions[valid_idx, 1] / oob_counts[valid_idx]
    r2 <- 1 - sum((y[valid_idx] - oob_pred_avg)^2) / sum((y[valid_idx] - mean(y[valid_idx]))^2)
    return(max(0, r2))  # Ensure non-negative R²
  }
}

#' Calculate Extra Trees Importance
#'
#' Calculates comprehensive feature importance for Extra Trees
#'
#' @param extra_trees_ensemble Trained Extra Trees ensemble
#' @param X Input data
#' @param y Target variable
#' @param is_classification Classification flag
#' @return Feature importance scores
calculate_extra_trees_importance <- function(extra_trees_ensemble, X, y, is_classification) {
  
  n_features <- ncol(X)
  
  # 1. Mean impurity decrease importance
  tree_importance <- extra_trees_ensemble$feature_importance
  mean_importance <- colMeans(tree_importance)
  
  # 2. Feature usage frequency
  usage_frequency <- extra_trees_ensemble$feature_usage_frequency
  
  # 3. Permutation importance (on subset for efficiency)
  perm_importance <- calculate_extra_trees_permutation_importance(extra_trees_ensemble, X, y, is_classification)
  
  # Normalize each importance measure
  mean_norm <- if(sum(mean_importance) > 0) mean_importance / sum(mean_importance) else rep(0, n_features)
  usage_norm <- if(sum(usage_frequency) > 0) usage_frequency / sum(usage_frequency) else rep(0, n_features)
  perm_norm <- if(sum(abs(perm_importance)) > 0) abs(perm_importance) / sum(abs(perm_importance)) else rep(0, n_features)
  
  # Combine importance measures
  combined_importance <- 0.5 * mean_norm + 0.3 * usage_norm + 0.2 * perm_norm
  
  return(combined_importance)
}

#' Calculate Extra Trees Permutation Importance
#'
#' Calculates permutation importance for Extra Trees
#'
#' @param extra_trees_ensemble Trained Extra Trees ensemble
#' @param X Input data
#' @param y Target variable
#' @param is_classification Classification flag
#' @return Permutation importance scores
calculate_extra_trees_permutation_importance <- function(extra_trees_ensemble, X, y, is_classification) {
  
  n_features <- ncol(X)
  
  # Use subset for efficiency
  subset_size <- min(200, nrow(X))
  subset_idx <- sample(nrow(X), subset_size)
  X_subset <- X[subset_idx, , drop = FALSE]
  y_subset <- y[subset_idx]
  
  # Get baseline predictions
  baseline_pred <- predict_extra_trees_ensemble(extra_trees_ensemble, X_subset, is_classification)
  baseline_score <- calculate_prediction_score(baseline_pred, y_subset, is_classification)
  
  perm_importance <- numeric(n_features)
  
  for (feature in 1:n_features) {
    
    # Permute feature
    X_perm <- X_subset
    X_perm[, feature] <- sample(X_perm[, feature])
    
    # Get predictions with permuted feature
    perm_pred <- predict_extra_trees_ensemble(extra_trees_ensemble, X_perm, is_classification)
    perm_score <- calculate_prediction_score(perm_pred, y_subset, is_classification)
    
    # Importance is decrease in score
    perm_importance[feature] <- baseline_score - perm_score
  }
  
  return(perm_importance)
}

#' Predict Extra Trees Ensemble
#'
#' Makes predictions using the entire Extra Trees ensemble
#'
#' @param extra_trees_ensemble Trained Extra Trees ensemble
#' @param X Test data
#' @param is_classification Classification flag
#' @return Ensemble predictions
predict_extra_trees_ensemble <- function(extra_trees_ensemble, X, is_classification) {
  
  n_samples <- nrow(X)
  n_trees <- extra_trees_ensemble$n_trees
  
  if (is_classification) {
    # For classification, use majority voting
    all_predictions <- matrix("", nrow = n_samples, ncol = n_trees)
    
    for (tree_idx in 1:n_trees) {
      tree_pred <- predict_extra_tree(extra_trees_ensemble$trees[[tree_idx]], X, is_classification)
      all_predictions[, tree_idx] <- tree_pred
    }
    
    # Majority vote for each sample
    final_predictions <- character(n_samples)
    for (i in 1:n_samples) {
      votes <- table(all_predictions[i, ])
      final_predictions[i] <- names(which.max(votes))
    }
    
    return(final_predictions)
    
  } else {
    # For regression, average predictions
    predictions <- numeric(n_samples)
    
    for (tree_idx in 1:n_trees) {
      tree_pred <- predict_extra_tree(extra_trees_ensemble$trees[[tree_idx]], X, is_classification)
      predictions <- predictions + tree_pred
    }
    
    return(predictions / n_trees)
  }
}

#' Calculate Prediction Score
#'
#' Calculates prediction score for importance calculation
#'
#' @param predictions Model predictions
#' @param y_true True labels
#' @param is_classification Classification flag
#' @return Prediction score
calculate_prediction_score <- function(predictions, y_true, is_classification) {
  
  if (is_classification) {
    return(mean(predictions == y_true))
  } else {
    # Use R² for regression
    ss_res <- sum((y_true - predictions)^2)
    ss_tot <- sum((y_true - mean(y_true))^2)
    r2 <- 1 - ss_res / ss_tot
    return(max(0, r2))
  }
}

#' Select Features by Extra Trees Importance
#'
#' Selects features based on Extra Trees importance scores
#'
#' @param importance_scores Feature importance scores
#' @param n_features Total number of features
#' @return Selected features and scores
select_features_by_extra_trees_importance <- function(importance_scores, n_features) {
  
  # Set selection threshold
  threshold <- quantile(importance_scores, 0.7)  # Top 30% by default
  
  # Ensure minimum selection
  min_features <- max(1, round(0.1 * n_features))
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

#' Create Empty Result for Extra Trees
#'
#' Creates standardized empty result for failed feature selection
#'
#' @param message Error or warning message
#' @return Empty result list
create_empty_result <- function(message) {
  return(list(
    selected_features = character(0),
    feature_scores = numeric(0),
    method_name = "Extra_Trees",
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
