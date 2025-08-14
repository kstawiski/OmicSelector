#' Gradient Boosting Feature Selection Method
#'
#' Implements feature selection using Gradient Boosting with built-in
#' feature importance measures. This method sequentially builds weak
#' learners and provides feature rankings based on split improvements.
#'
#' @param ramwas.data Expression data matrix (samples as rows, features as columns)  
#' @param target Target variable vector
#' @param smote_sampling Logical, whether to apply SMOTE sampling
#' @return List containing selected features and method metadata
#' @export
gradient_boosting_fs <- function(ramwas.data, target, smote_sampling = FALSE) {
  
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
    
    # Perform Gradient Boosting feature selection
    gb_result <- perform_gradient_boosting_selection(data_clean, target)
    
    # Extract selected features
    if (length(gb_result$selected_features) == 0) {
      return(create_empty_result("No features selected by Gradient Boosting"))
    }
    
    # Calculate timing
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      selected_features = gb_result$selected_features,
      feature_scores = gb_result$feature_scores,
      method_name = "Gradient_Boosting",
      execution_time = execution_time,
      n_features_selected = length(gb_result$selected_features),
      n_features_input = ncol(ramwas.data),
      method_params = list(
        n_estimators = gb_result$n_estimators,
        learning_rate = gb_result$learning_rate,
        max_depth = gb_result$max_depth,
        subsample = gb_result$subsample,
        feature_fraction = gb_result$feature_fraction,
        smote_sampling = smote_sampling
      ),
      preprocessing_info = list(
        constant_features_removed = preprocessed$constant_removed,
        correlated_features_removed = preprocessed$correlated_removed,
        final_features_count = ncol(data_clean)
      ),
      gb_info = gb_result$gb_info
    ))
    
  }, error = function(e) {
    # Return error result
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      selected_features = character(0),
      feature_scores = numeric(0),
      method_name = "Gradient_Boosting",
      execution_time = execution_time,
      n_features_selected = 0,
      n_features_input = ncol(ramwas.data),
      error_message = paste("Gradient Boosting failed:", e$message),
      method_params = list(smote_sampling = smote_sampling)
    ))
  })
}

#' Perform Gradient Boosting Feature Selection
#'
#' Implements Gradient Boosting with feature importance calculation
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with selected features and GB information
perform_gradient_boosting_selection <- function(data_matrix, target_vector) {
  
  n_features <- ncol(data_matrix)
  
  # Determine if classification or regression
  is_classification <- determine_task_type(target_vector)
  
  # Set Gradient Boosting parameters
  n_estimators <- max(50, min(200, round(sqrt(n_features) * 5)))
  learning_rate <- 0.1
  max_depth <- max(2, min(6, round(log2(n_features))))
  subsample <- 0.8
  feature_fraction <- min(1.0, max(0.3, sqrt(n_features) / n_features))
  
  # Standardize features
  data_scaled <- scale(data_matrix)
  
  # Handle any NAs from scaling
  if (any(is.na(data_scaled))) {
    non_na_cols <- apply(data_scaled, 2, function(x) !any(is.na(x)))
    data_scaled <- data_scaled[, non_na_cols, drop = FALSE]
    n_features <- ncol(data_scaled)
  }
  
  # Train Gradient Boosting model
  gb_model <- train_gradient_boosting(data_scaled, target_vector, n_estimators, 
                                     learning_rate, max_depth, subsample, 
                                     feature_fraction, is_classification)
  
  # Calculate feature importance
  feature_importance <- calculate_gb_feature_importance(gb_model, data_scaled, target_vector, is_classification)
  
  # Select features based on importance
  selected_result <- select_features_by_gb_importance(feature_importance, n_features)
  
  return(list(
    selected_features = selected_result$selected_features,
    feature_scores = selected_result$feature_scores,
    n_estimators = n_estimators,
    learning_rate = learning_rate,
    max_depth = max_depth,
    subsample = subsample,
    feature_fraction = feature_fraction,
    gb_info = list(
      task_type = if(is_classification) "classification" else "regression",
      final_loss = gb_model$final_loss,
      n_trees = length(gb_model$trees),
      feature_usage_count = gb_model$feature_usage_count
    )
  ))
}

#' Train Gradient Boosting Model
#'
#' Trains gradient boosting ensemble with weak learners
#'
#' @param X Input data matrix
#' @param y Target vector
#' @param n_estimators Number of boosting rounds
#' @param learning_rate Learning rate for updates
#' @param max_depth Maximum depth of weak learners
#' @param subsample Sample fraction for each round
#' @param feature_fraction Feature fraction for each round
#' @param is_classification Classification flag
#' @return Trained gradient boosting model
train_gradient_boosting <- function(X, y, n_estimators, learning_rate, max_depth, 
                                   subsample, feature_fraction, is_classification) {
  
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  # Initialize model
  trees <- list()
  feature_importance <- numeric(n_features)
  feature_usage_count <- integer(n_features)
  
  # Initialize predictions and residuals
  if (is_classification) {
    # For binary classification, use log-odds
    class_labels <- unique(y)
    if (length(class_labels) > 2) {
      stop("Multi-class classification not implemented")
    }
    
    y_binary <- as.numeric(y == class_labels[2])
    init_pred <- log(mean(y_binary) / (1 - mean(y_binary)))
    predictions <- rep(init_pred, n_samples)
    
    # Calculate initial probabilities and residuals
    probs <- 1 / (1 + exp(-predictions))
    residuals <- y_binary - probs
    
  } else {
    # For regression, use mean
    init_pred <- mean(y)
    predictions <- rep(init_pred, n_samples)
    residuals <- y - predictions
  }
  
  current_loss <- calculate_loss(y, predictions, is_classification)
  
  # Boosting iterations
  for (iter in 1:n_estimators) {
    
    # Subsample data
    sample_idx <- sample(n_samples, round(subsample * n_samples))
    
    # Select features
    n_features_use <- max(1, round(feature_fraction * n_features))
    feature_idx <- sample(n_features, n_features_use)
    
    X_sub <- X[sample_idx, feature_idx, drop = FALSE]
    residuals_sub <- residuals[sample_idx]
    
    # Train weak learner on residuals
    weak_learner <- train_weak_learner(X_sub, residuals_sub, feature_idx, max_depth)
    
    if (is.null(weak_learner)) next
    
    # Make predictions with weak learner
    weak_predictions <- predict_weak_learner(weak_learner, X)
    
    # Update predictions
    predictions <- predictions + learning_rate * weak_predictions
    
    # Update residuals
    if (is_classification) {
      probs <- 1 / (1 + exp(-predictions))
      residuals <- y_binary - probs
    } else {
      residuals <- y - predictions
    }
    
    # Store tree and update importance
    trees[[iter]] <- weak_learner
    feature_importance[weak_learner$feature_indices] <- 
      feature_importance[weak_learner$feature_indices] + weak_learner$importance_gains
    feature_usage_count[weak_learner$feature_indices] <- 
      feature_usage_count[weak_learner$feature_indices] + 1
    
    # Calculate current loss
    new_loss <- calculate_loss(y, predictions, is_classification)
    
    # Early stopping if no improvement
    if (iter > 10 && abs(new_loss - current_loss) < 1e-6) {
      break
    }
    
    current_loss <- new_loss
  }
  
  return(list(
    trees = trees,
    feature_importance = feature_importance,
    feature_usage_count = feature_usage_count,
    final_loss = current_loss,
    learning_rate = learning_rate,
    is_classification = is_classification
  ))
}

#' Train Weak Learner
#'
#' Trains a single weak learner (decision stump or small tree)
#'
#' @param X Training data
#' @param residuals Target residuals
#' @param feature_indices Available feature indices
#' @param max_depth Maximum tree depth
#' @return Weak learner model
train_weak_learner <- function(X, residuals, feature_indices, max_depth) {
  
  if (nrow(X) < 2 || length(unique(residuals)) == 1) {
    return(NULL)
  }
  
  # Build regression tree on residuals
  tree <- build_regression_tree(X, residuals, max_depth)
  
  if (is.null(tree)) {
    return(NULL)
  }
  
  # Calculate feature importance gains
  importance_gains <- calculate_tree_importance(tree, ncol(X))
  
  return(list(
    tree = tree,
    feature_indices = feature_indices,
    importance_gains = importance_gains
  ))
}

#' Build Regression Tree
#'
#' Builds a regression tree for residuals
#'
#' @param X Training data
#' @param y Target values (residuals)
#' @param max_depth Maximum depth
#' @return Regression tree
build_regression_tree <- function(X, y, max_depth) {
  
  tree <- build_regression_node(X, y, 1, max_depth)
  return(tree)
}

#' Build Regression Node
#'
#' Recursively builds regression tree nodes
#'
#' @param X Node data
#' @param y Node targets
#' @param depth Current depth
#' @param max_depth Maximum depth
#' @return Tree node
build_regression_node <- function(X, y, depth, max_depth) {
  
  n_samples <- nrow(X)
  
  # Stopping criteria
  if (depth >= max_depth || n_samples < 5 || var(y) < 1e-6) {
    return(list(
      type = "leaf",
      prediction = mean(y),
      n_samples = n_samples,
      variance = var(y)
    ))
  }
  
  # Find best split
  best_split <- find_best_regression_split(X, y)
  
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
  left_node <- build_regression_node(X[left_idx, , drop = FALSE], y[left_idx], depth + 1, max_depth)
  right_node <- build_regression_node(X[right_idx, , drop = FALSE], y[right_idx], depth + 1, max_depth)
  
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

#' Find Best Regression Split
#'
#' Finds the best split for regression tree
#'
#' @param X Node data
#' @param y Node targets
#' @return Best split information
find_best_regression_split <- function(X, y) {
  
  n_features <- ncol(X)
  best_score <- -Inf
  best_split <- NULL
  
  for (feature in 1:n_features) {
    
    feature_values <- X[, feature]
    unique_values <- sort(unique(feature_values))
    
    if (length(unique_values) < 2) next
    
    # Try different thresholds
    for (i in 1:(length(unique_values) - 1)) {
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

#' Calculate Tree Importance
#'
#' Calculates feature importance for a single tree
#'
#' @param tree Regression tree
#' @param n_features Number of features
#' @return Feature importance vector
calculate_tree_importance <- function(tree, n_features) {
  
  importance <- numeric(n_features)
  
  if (tree$type == "internal") {
    importance[tree$feature] <- tree$importance
    
    # Recursively add importance from child nodes
    left_importance <- calculate_tree_importance(tree$left, n_features)
    right_importance <- calculate_tree_importance(tree$right, n_features)
    
    importance <- importance + left_importance + right_importance
  }
  
  return(importance)
}

#' Predict Weak Learner
#'
#' Makes predictions using a weak learner
#'
#' @param weak_learner Weak learner model
#' @param X Test data
#' @return Predictions
predict_weak_learner <- function(weak_learner, X) {
  
  n_samples <- nrow(X)
  predictions <- numeric(n_samples)
  
  # Map feature indices
  X_mapped <- X[, weak_learner$feature_indices, drop = FALSE]
  
  for (i in 1:n_samples) {
    predictions[i] <- predict_regression_tree(weak_learner$tree, X_mapped[i, ])
  }
  
  return(predictions)
}

#' Predict Regression Tree
#'
#' Predicts using regression tree
#'
#' @param tree Regression tree
#' @param sample Single sample
#' @return Prediction
predict_regression_tree <- function(tree, sample) {
  
  if (tree$type == "leaf") {
    return(tree$prediction)
  }
  
  if (sample[tree$feature] <= tree$threshold) {
    return(predict_regression_tree(tree$left, sample))
  } else {
    return(predict_regression_tree(tree$right, sample))
  }
}

#' Calculate Loss
#'
#' Calculates loss for gradient boosting
#'
#' @param y_true True labels
#' @param y_pred Predictions
#' @param is_classification Classification flag
#' @return Loss value
calculate_loss <- function(y_true, y_pred, is_classification) {
  
  if (is_classification) {
    # Logistic loss for binary classification
    class_labels <- unique(y_true)
    y_binary <- as.numeric(y_true == class_labels[2])
    
    # Convert to probabilities
    probs <- 1 / (1 + exp(-y_pred))
    probs <- pmax(pmin(probs, 1 - 1e-15), 1e-15)  # Numerical stability
    
    loss <- -mean(y_binary * log(probs) + (1 - y_binary) * log(1 - probs))
    
  } else {
    # Mean squared error for regression
    loss <- mean((y_true - y_pred)^2)
  }
  
  return(loss)
}

#' Calculate GB Feature Importance
#'
#' Calculates comprehensive feature importance for gradient boosting
#'
#' @param gb_model Trained GB model
#' @param X Input data
#' @param y Target variable
#' @param is_classification Classification flag
#' @return Feature importance scores
calculate_gb_feature_importance <- function(gb_model, X, y, is_classification) {
  
  n_features <- ncol(X)
  
  # 1. Split-based importance (accumulated from trees)
  split_importance <- gb_model$feature_importance
  
  # 2. Usage-based importance (frequency of use)
  usage_importance <- gb_model$feature_usage_count
  
  # 3. Permutation importance
  perm_importance <- calculate_gb_permutation_importance(gb_model, X, y, is_classification)
  
  # Normalize importance measures
  split_norm <- if(sum(split_importance) > 0) split_importance / sum(split_importance) else rep(0, n_features)
  usage_norm <- if(sum(usage_importance) > 0) usage_importance / sum(usage_importance) else rep(0, n_features)
  perm_norm <- if(sum(abs(perm_importance)) > 0) abs(perm_importance) / sum(abs(perm_importance)) else rep(0, n_features)
  
  # Combine importance measures
  combined_importance <- 0.5 * split_norm + 0.3 * usage_norm + 0.2 * perm_norm
  
  return(combined_importance)
}

#' Calculate GB Permutation Importance
#'
#' Calculates permutation importance for gradient boosting
#'
#' @param gb_model Trained GB model
#' @param X Input data
#' @param y Target variable
#' @param is_classification Classification flag
#' @return Permutation importance scores
calculate_gb_permutation_importance <- function(gb_model, X, y, is_classification) {
  
  n_features <- ncol(X)
  
  # Get baseline predictions
  baseline_pred <- predict_gb_model(gb_model, X)
  baseline_loss <- calculate_loss(y, baseline_pred, is_classification)
  
  perm_importance <- numeric(n_features)
  
  for (feature in 1:n_features) {
    
    # Permute feature
    X_perm <- X
    X_perm[, feature] <- sample(X_perm[, feature])
    
    # Get predictions with permuted feature
    perm_pred <- predict_gb_model(gb_model, X_perm)
    perm_loss <- calculate_loss(y, perm_pred, is_classification)
    
    # Importance is increase in loss
    perm_importance[feature] <- perm_loss - baseline_loss
  }
  
  return(perm_importance)
}

#' Predict GB Model
#'
#' Makes predictions using the full gradient boosting model
#'
#' @param gb_model Trained GB model
#' @param X Test data
#' @return Predictions
predict_gb_model <- function(gb_model, X) {
  
  n_samples <- nrow(X)
  predictions <- numeric(n_samples)
  
  # Initialize with base prediction
  if (gb_model$is_classification) {
    # Initialize with log-odds
    base_prob <- 0.5  # Assume balanced initialization
    predictions <- rep(log(base_prob / (1 - base_prob)), n_samples)
  } else {
    # Initialize with mean (would need to store this in practice)
    predictions <- rep(0, n_samples)
  }
  
  # Add predictions from all trees
  for (tree in gb_model$trees) {
    tree_pred <- predict_weak_learner(tree, X)
    predictions <- predictions + gb_model$learning_rate * tree_pred
  }
  
  return(predictions)
}

#' Select Features by GB Importance
#'
#' Selects features based on gradient boosting importance
#'
#' @param importance_scores Feature importance scores
#' @param n_features Total number of features
#' @return Selected features and scores
select_features_by_gb_importance <- function(importance_scores, n_features) {
  
  # Set selection threshold
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

#' Create Empty Result for Gradient Boosting
#'
#' Creates standardized empty result for failed feature selection
#'
#' @param message Error or warning message
#' @return Empty result list
create_empty_result <- function(message) {
  return(list(
    selected_features = character(0),
    feature_scores = numeric(0),
    method_name = "Gradient_Boosting",
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
