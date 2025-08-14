#' XGBoost Feature Selection
#' 
#' This function implements XGBoost (eXtreme Gradient Boosting) for feature selection.
#' XGBoost is an optimized gradient boosting framework that uses second-order gradients
#' and regularization for better performance.
#' 
#' @param X A matrix or data frame of features (predictors).
#' @param y A vector of target values (binary classification: 0 and 1, or factor/character).
#' @param num_features Number of features to select (default: 10).
#' @param n_estimators Number of boosting rounds (default: 100).
#' @param learning_rate Step size shrinkage (default: 0.1).
#' @param max_depth Maximum depth of trees (default: 6).
#' @param min_child_weight Minimum sum of instance weight in a child (default: 1).
#' @param subsample Subsample ratio for training instances (default: 1.0).
#' @param colsample_bytree Subsample ratio of columns for each tree (default: 1.0).
#' @param reg_alpha L1 regularization parameter (default: 0).
#' @param reg_lambda L2 regularization parameter (default: 1).
#' @param cv_folds Number of cross-validation folds (default: 5).
#' @param random_state Random seed for reproducibility (default: 42).
#' @param verbose Logical, whether to print progress (default: FALSE).
#' @param use_smote Logical, whether to apply SMOTE for class imbalance (default: FALSE).
#' 
#' @return A list containing:
#'   \item{selected_features}{Names of selected features}
#'   \item{importance_scores}{Feature importance scores}
#'   \item{model_params}{Model parameters used}
#'   \item{method_info}{Information about the method}
#'   \item{cv_performance}{Cross-validation performance metrics}
#'   \item{feature_gain}{Feature importance based on gain}
#'   \item{feature_cover}{Feature importance based on cover}
#'   \item{feature_frequency}{Feature usage frequency}
#' 
#' @examples
#' \dontrun{
#' # Generate sample data
#' X <- matrix(rnorm(1000), ncol = 20)
#' y <- sample(0:1, 50, replace = TRUE)
#' 
#' # Apply XGBoost feature selection
#' result <- fs_method_xgboost(X, y, num_features = 5)
#' print(result$selected_features)
#' }
#' 
#' @export
fs_method_xgboost <- function(X, y, num_features = 10, n_estimators = 100, 
                              learning_rate = 0.1, max_depth = 6,
                              min_child_weight = 1, subsample = 1.0,
                              colsample_bytree = 1.0, reg_alpha = 0,
                              reg_lambda = 1, cv_folds = 5, 
                              random_state = 42, verbose = FALSE, 
                              use_smote = FALSE) {
  
  # Set random seed for reproducibility
  set.seed(random_state)
  
  # Validate inputs
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data frame")
  }
  if (length(y) != nrow(X)) {
    stop("Length of y must equal number of rows in X")
  }
  if (num_features <= 0 || num_features > ncol(X)) {
    stop("num_features must be between 1 and ncol(X)")
  }
  if (n_estimators <= 0) {
    stop("n_estimators must be positive")
  }
  if (learning_rate <= 0 || learning_rate > 1) {
    stop("learning_rate must be between 0 and 1")
  }
  
  # Convert inputs to standard format
  X <- as.matrix(X)
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  # Handle feature names
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("Feature_", seq_len(n_features))
  }
  feature_names <- colnames(X)
  
  # Convert target to binary
  if (is.factor(y) || is.character(y)) {
    unique_classes <- unique(y)
    if (length(unique_classes) != 2) {
      stop("XGBoost requires exactly 2 classes")
    }
    y <- as.numeric(as.factor(y)) - 1
  }
  
  # Ensure target is 0/1
  if (!all(y %in% c(0, 1))) {
    unique_vals <- sort(unique(y))
    if (length(unique_vals) == 2) {
      y[y == unique_vals[1]] <- 0
      y[y == unique_vals[2]] <- 1
    } else {
      stop("Target must be binary")
    }
  }
  
  # Apply SMOTE if requested
  original_X <- X
  original_y <- y
  
  if (use_smote) {
    tryCatch({
      smote_result <- apply_simple_smote(X, y)
      X <- smote_result$X
      y <- smote_result$y
      if (verbose) {
        cat("Applied SMOTE: samples increased from", n_samples, "to", nrow(X), "\n")
      }
    }, error = function(e) {
      if (verbose) {
        cat("SMOTE failed, using original data:", e$message, "\n")
      }
    })
  }
  
  # Update dimensions after potential SMOTE
  n_samples <- nrow(X)
  
  # Train XGBoost model
  if (verbose) {
    cat("Training XGBoost with", n_estimators, "rounds...\n")
  }
  
  xgb_model <- train_xgboost_model(X, y, n_estimators, learning_rate, 
                                   max_depth, min_child_weight, subsample,
                                   colsample_bytree, reg_alpha, reg_lambda, 
                                   verbose)
  
  # Calculate feature importance
  importance_metrics <- calculate_xgboost_importance(xgb_model, feature_names)
  
  # Cross-validation performance assessment
  cv_results <- perform_xgboost_cv(original_X, original_y, cv_folds, 
                                   n_estimators, learning_rate, max_depth,
                                   min_child_weight, subsample, colsample_bytree,
                                   reg_alpha, reg_lambda, use_smote, verbose)
  
  # Combine importance scores (weighted average of gain, cover, and frequency)
  final_importance <- 0.5 * importance_metrics$gain + 
                      0.3 * importance_metrics$cover + 
                      0.2 * importance_metrics$frequency
  
  # Select top features
  top_indices <- order(final_importance, decreasing = TRUE)[seq_len(min(num_features, length(final_importance)))]
  selected_features <- feature_names[top_indices]
  selected_importance <- final_importance[top_indices]
  names(selected_importance) <- selected_features
  
  # Prepare results
  result <- list(
    selected_features = selected_features,
    importance_scores = selected_importance,
    model_params = list(
      n_estimators = n_estimators,
      learning_rate = learning_rate,
      max_depth = max_depth,
      min_child_weight = min_child_weight,
      subsample = subsample,
      colsample_bytree = colsample_bytree,
      reg_alpha = reg_alpha,
      reg_lambda = reg_lambda,
      cv_folds = cv_folds,
      random_state = random_state,
      use_smote = use_smote
    ),
    method_info = list(
      method_name = "XGBoost",
      method_type = "Ensemble",
      description = "eXtreme Gradient Boosting with regularization",
      references = "Chen, T. & Guestrin, C. (2016). XGBoost: A scalable tree boosting system."
    ),
    cv_performance = cv_results,
    feature_gain = importance_metrics$gain[top_indices],
    feature_cover = importance_metrics$cover[top_indices],
    feature_frequency = importance_metrics$frequency[top_indices],
    xgb_model = xgb_model
  )
  
  return(result)
}

#' Train XGBoost Model
#' 
#' @param X Feature matrix
#' @param y Target vector
#' @param n_estimators Number of rounds
#' @param learning_rate Learning rate
#' @param max_depth Maximum depth
#' @param min_child_weight Minimum child weight
#' @param subsample Subsample ratio
#' @param colsample_bytree Column sample ratio
#' @param reg_alpha L1 regularization
#' @param reg_lambda L2 regularization
#' @param verbose Whether to print progress
#' @return Trained XGBoost model
train_xgboost_model <- function(X, y, n_estimators, learning_rate, max_depth,
                                min_child_weight, subsample, colsample_bytree,
                                reg_alpha, reg_lambda, verbose) {
  
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  # Initialize gradients and hessians
  predictions <- rep(0.5, n_samples)  # Initialize with 0.5 probability
  trees <- list()
  
  for (round in 1:n_estimators) {
    if (verbose && round %% 20 == 0) {
      cat("Round", round, "of", n_estimators, "\n")
    }
    
    # Calculate gradients and hessians
    probs <- 1 / (1 + exp(-predictions))
    gradients <- probs - y
    hessians <- probs * (1 - probs)
    
    # Subsample data if needed
    if (subsample < 1.0) {
      sample_indices <- sample(seq_len(n_samples), 
                               floor(subsample * n_samples))
      X_subsample <- X[sample_indices, , drop = FALSE]
      gradients_subsample <- gradients[sample_indices]
      hessians_subsample <- hessians[sample_indices]
    } else {
      X_subsample <- X
      gradients_subsample <- gradients
      hessians_subsample <- hessians
      sample_indices <- seq_len(n_samples)
    }
    
    # Sample features for this tree
    if (colsample_bytree < 1.0) {
      feature_indices <- sample(seq_len(n_features), 
                                floor(colsample_bytree * n_features))
    } else {
      feature_indices <- seq_len(n_features)
    }
    
    # Build tree
    tree <- build_xgboost_tree(X_subsample[, feature_indices, drop = FALSE], 
                               gradients_subsample, hessians_subsample,
                               max_depth, min_child_weight, reg_alpha, 
                               reg_lambda, feature_indices)
    
    # Update predictions
    tree_predictions <- predict_xgboost_tree(tree, X)
    predictions <- predictions + learning_rate * tree_predictions
    
    # Store tree
    trees[[round]] <- tree
  }
  
  return(list(
    trees = trees,
    learning_rate = learning_rate,
    n_estimators = length(trees)
  ))
}

#' Build XGBoost Tree
#' 
#' @param X Feature matrix
#' @param gradients Gradient vector
#' @param hessians Hessian vector
#' @param max_depth Maximum depth
#' @param min_child_weight Minimum child weight
#' @param reg_alpha L1 regularization
#' @param reg_lambda L2 regularization
#' @param feature_indices Original feature indices
#' @return Tree structure
build_xgboost_tree <- function(X, gradients, hessians, max_depth, 
                               min_child_weight, reg_alpha, reg_lambda,
                               feature_indices) {
  
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  # Calculate leaf value for current node
  G <- sum(gradients)
  H <- sum(hessians)
  
  if (H < min_child_weight) {
    H <- min_child_weight
  }
  
  leaf_value <- -G / (H + reg_lambda)
  
  # Base case: max depth reached or not enough samples
  if (max_depth <= 0 || n_samples < 2) {
    return(list(
      is_leaf = TRUE,
      leaf_value = leaf_value,
      feature_usage = numeric(length(feature_indices))
    ))
  }
  
  # Find best split
  best_split <- find_best_xgboost_split(X, gradients, hessians, 
                                        min_child_weight, reg_alpha, 
                                        reg_lambda)
  
  # If no good split found, return leaf
  if (is.null(best_split) || best_split$gain <= 0) {
    return(list(
      is_leaf = TRUE,
      leaf_value = leaf_value,
      feature_usage = numeric(length(feature_indices))
    ))
  }
  
  # Split data
  left_mask <- X[, best_split$feature] <= best_split$threshold
  right_mask <- !left_mask
  
  # Initialize feature usage
  feature_usage <- numeric(length(feature_indices))
  feature_usage[best_split$feature] <- 1
  
  # Build child nodes
  left_child <- build_xgboost_tree(X[left_mask, , drop = FALSE],
                                   gradients[left_mask], hessians[left_mask],
                                   max_depth - 1, min_child_weight,
                                   reg_alpha, reg_lambda, feature_indices)
  
  right_child <- build_xgboost_tree(X[right_mask, , drop = FALSE],
                                    gradients[right_mask], hessians[right_mask],
                                    max_depth - 1, min_child_weight,
                                    reg_alpha, reg_lambda, feature_indices)
  
  # Combine feature usage
  feature_usage <- feature_usage + left_child$feature_usage + right_child$feature_usage
  
  return(list(
    is_leaf = FALSE,
    feature = best_split$feature,
    threshold = best_split$threshold,
    gain = best_split$gain,
    left_child = left_child,
    right_child = right_child,
    feature_usage = feature_usage,
    original_feature_idx = feature_indices[best_split$feature]
  ))
}

#' Find Best XGBoost Split
#' 
#' @param X Feature matrix
#' @param gradients Gradient vector
#' @param hessians Hessian vector
#' @param min_child_weight Minimum child weight
#' @param reg_alpha L1 regularization
#' @param reg_lambda L2 regularization
#' @return Best split information
find_best_xgboost_split <- function(X, gradients, hessians, min_child_weight,
                                    reg_alpha, reg_lambda) {
  
  n_features <- ncol(X)
  
  best_gain <- 0
  best_feature <- NULL
  best_threshold <- NULL
  
  for (feature in seq_len(n_features)) {
    feature_values <- X[, feature]
    unique_values <- sort(unique(feature_values))
    
    if (length(unique_values) < 2) {
      next
    }
    
    for (i in seq_len(length(unique_values) - 1)) {
      threshold <- (unique_values[i] + unique_values[i + 1]) / 2
      
      left_mask <- feature_values <= threshold
      right_mask <- !left_mask
      
      # Calculate split statistics
      GL <- sum(gradients[left_mask])
      HL <- sum(hessians[left_mask])
      GR <- sum(gradients[right_mask])
      HR <- sum(hessians[right_mask])
      
      # Check minimum child weight constraint
      if (HL < min_child_weight || HR < min_child_weight) {
        next
      }
      
      # Calculate gain with regularization
      G <- GL + GR
      H <- HL + HR + reg_lambda
      
      left_score <- GL^2 / (HL + reg_lambda)
      right_score <- GR^2 / (HR + reg_lambda)
      parent_score <- G^2 / H
      
      gain <- 0.5 * (left_score + right_score - parent_score) - reg_alpha
      
      if (gain > best_gain) {
        best_gain <- gain
        best_feature <- feature
        best_threshold <- threshold
      }
    }
  }
  
  if (is.null(best_feature)) {
    return(NULL)
  }
  
  return(list(
    feature = best_feature,
    threshold = best_threshold,
    gain = best_gain
  ))
}

#' Predict with XGBoost Tree
#' 
#' @param tree Tree structure
#' @param X Feature matrix
#' @return Predictions
predict_xgboost_tree <- function(tree, X) {
  predictions <- numeric(nrow(X))
  
  for (i in seq_len(nrow(X))) {
    predictions[i] <- predict_single_xgboost_tree(tree, X[i, ])
  }
  
  return(predictions)
}

#' Predict Single Sample with XGBoost Tree
#' 
#' @param tree Tree structure
#' @param x Single sample
#' @return Prediction
predict_single_xgboost_tree <- function(tree, x) {
  if (tree$is_leaf) {
    return(tree$leaf_value)
  }
  
  if (x[tree$feature] <= tree$threshold) {
    return(predict_single_xgboost_tree(tree$left_child, x))
  } else {
    return(predict_single_xgboost_tree(tree$right_child, x))
  }
}

#' Calculate XGBoost Feature Importance
#' 
#' @param model Trained XGBoost model
#' @param feature_names Feature names
#' @return Importance metrics
calculate_xgboost_importance <- function(model, feature_names) {
  n_features <- length(feature_names)
  
  gain_importance <- numeric(n_features)
  cover_importance <- numeric(n_features)
  frequency_importance <- numeric(n_features)
  
  names(gain_importance) <- feature_names
  names(cover_importance) <- feature_names
  names(frequency_importance) <- feature_names
  
  for (tree in model$trees) {
    tree_importance <- calculate_tree_importance(tree, n_features)
    
    gain_importance <- gain_importance + tree_importance$gain
    cover_importance <- cover_importance + tree_importance$cover
    frequency_importance <- frequency_importance + tree_importance$frequency
  }
  
  # Normalize importance scores
  if (sum(gain_importance) > 0) {
    gain_importance <- gain_importance / sum(gain_importance)
  }
  if (sum(cover_importance) > 0) {
    cover_importance <- cover_importance / sum(cover_importance)
  }
  if (sum(frequency_importance) > 0) {
    frequency_importance <- frequency_importance / sum(frequency_importance)
  }
  
  return(list(
    gain = gain_importance,
    cover = cover_importance,
    frequency = frequency_importance
  ))
}

#' Calculate Tree Importance
#' 
#' @param tree Tree structure
#' @param n_features Number of features
#' @return Tree importance metrics
calculate_tree_importance <- function(tree, n_features) {
  gain_importance <- numeric(n_features)
  cover_importance <- numeric(n_features)
  frequency_importance <- numeric(n_features)
  
  if (!tree$is_leaf) {
    # Current node contribution
    if (!is.null(tree$original_feature_idx)) {
      feature_idx <- tree$original_feature_idx
      gain_importance[feature_idx] <- gain_importance[feature_idx] + tree$gain
      cover_importance[feature_idx] <- cover_importance[feature_idx] + 1
      frequency_importance[feature_idx] <- frequency_importance[feature_idx] + 1
    }
    
    # Recursive calculation for children
    left_importance <- calculate_tree_importance(tree$left_child, n_features)
    right_importance <- calculate_tree_importance(tree$right_child, n_features)
    
    gain_importance <- gain_importance + left_importance$gain + right_importance$gain
    cover_importance <- cover_importance + left_importance$cover + right_importance$cover
    frequency_importance <- frequency_importance + left_importance$frequency + right_importance$frequency
  }
  
  return(list(
    gain = gain_importance,
    cover = cover_importance,
    frequency = frequency_importance
  ))
}

#' Perform Cross-Validation for XGBoost
#' 
#' @param X Feature matrix
#' @param y Target vector
#' @param cv_folds Number of CV folds
#' @param n_estimators Number of estimators
#' @param learning_rate Learning rate
#' @param max_depth Maximum depth
#' @param min_child_weight Minimum child weight
#' @param subsample Subsample ratio
#' @param colsample_bytree Column sample ratio
#' @param reg_alpha L1 regularization
#' @param reg_lambda L2 regularization
#' @param use_smote Whether to use SMOTE
#' @param verbose Whether to print progress
#' @return CV performance metrics
perform_xgboost_cv <- function(X, y, cv_folds, n_estimators, learning_rate,
                               max_depth, min_child_weight, subsample,
                               colsample_bytree, reg_alpha, reg_lambda,
                               use_smote, verbose) {
  
  n_samples <- nrow(X)
  fold_size <- floor(n_samples / cv_folds)
  indices <- sample(seq_len(n_samples))
  
  cv_accuracies <- numeric(cv_folds)
  cv_aucs <- numeric(cv_folds)
  
  for (fold in seq_len(cv_folds)) {
    # Create train/test split
    test_start <- (fold - 1) * fold_size + 1
    test_end <- ifelse(fold == cv_folds, n_samples, fold * fold_size)
    test_indices <- indices[test_start:test_end]
    train_indices <- indices[-c(test_start:test_end)]
    
    X_train <- X[train_indices, , drop = FALSE]
    y_train <- y[train_indices]
    X_test <- X[test_indices, , drop = FALSE]
    y_test <- y[test_indices]
    
    # Apply SMOTE to training data if requested
    if (use_smote) {
      tryCatch({
        smote_result <- apply_simple_smote(X_train, y_train)
        X_train <- smote_result$X
        y_train <- smote_result$y
      }, error = function(e) {
        # Continue with original data if SMOTE fails
      })
    }
    
    # Train XGBoost model
    model <- train_xgboost_model(X_train, y_train, n_estimators, 
                                 learning_rate, max_depth, min_child_weight,
                                 subsample, colsample_bytree, reg_alpha, 
                                 reg_lambda, FALSE)
    
    # Make predictions
    predictions <- predict_xgboost_model(model, X_test)
    probabilities <- predict_xgboost_probabilities(model, X_test)
    
    # Calculate metrics
    cv_accuracies[fold] <- mean(predictions == y_test)
    cv_aucs[fold] <- calculate_auc(y_test, probabilities)
  }
  
  return(list(
    mean_accuracy = mean(cv_accuracies),
    std_accuracy = sd(cv_accuracies),
    mean_auc = mean(cv_aucs, na.rm = TRUE),
    std_auc = sd(cv_aucs, na.rm = TRUE),
    fold_accuracies = cv_accuracies,
    fold_aucs = cv_aucs
  ))
}

#' Predict with XGBoost Model
#' 
#' @param model Trained XGBoost model
#' @param X Test features
#' @return Predictions
predict_xgboost_model <- function(model, X) {
  n_samples <- nrow(X)
  predictions <- numeric(n_samples)
  
  for (tree in model$trees) {
    tree_predictions <- predict_xgboost_tree(tree, X)
    predictions <- predictions + model$learning_rate * tree_predictions
  }
  
  # Convert to binary predictions
  probabilities <- 1 / (1 + exp(-predictions))
  binary_predictions <- ifelse(probabilities > 0.5, 1, 0)
  
  return(binary_predictions)
}

#' Predict Probabilities with XGBoost Model
#' 
#' @param model Trained XGBoost model
#' @param X Test features
#' @return Prediction probabilities
predict_xgboost_probabilities <- function(model, X) {
  n_samples <- nrow(X)
  predictions <- numeric(n_samples)
  
  for (tree in model$trees) {
    tree_predictions <- predict_xgboost_tree(tree, X)
    predictions <- predictions + model$learning_rate * tree_predictions
  }
  
  # Convert to probabilities
  probabilities <- 1 / (1 + exp(-predictions))
  return(probabilities)
}

# Helper functions (SMOTE and AUC) would be included here if not already defined
# in other files - using the same implementations as in AdaBoost

#' Simple SMOTE Implementation
#' 
#' @param X Feature matrix
#' @param y Target vector
#' @return List with balanced X and y
apply_simple_smote <- function(X, y) {
  # Count classes
  class_counts <- table(y)
  minority_class <- names(class_counts)[which.min(class_counts)]
  majority_class <- names(class_counts)[which.max(class_counts)]
  
  # If classes are balanced, return original data
  if (min(class_counts) / max(class_counts) > 0.8) {
    return(list(X = X, y = y))
  }
  
  # Get minority and majority samples
  minority_indices <- which(y == minority_class)
  majority_indices <- which(y == majority_class)
  
  # Calculate how many synthetic samples to generate
  n_synthetic <- length(majority_indices) - length(minority_indices)
  
  if (n_synthetic <= 0) {
    return(list(X = X, y = y))
  }
  
  # Generate synthetic samples
  synthetic_X <- matrix(0, nrow = n_synthetic, ncol = ncol(X))
  synthetic_y <- rep(minority_class, n_synthetic)
  
  for (i in seq_len(n_synthetic)) {
    # Random minority sample
    idx1 <- sample(minority_indices, 1)
    idx2 <- sample(minority_indices, 1)
    
    # Generate synthetic sample
    lambda <- runif(1)
    synthetic_X[i, ] <- lambda * X[idx1, ] + (1 - lambda) * X[idx2, ]
  }
  
  # Combine original and synthetic data
  new_X <- rbind(X, synthetic_X)
  new_y <- c(y, synthetic_y)
  
  return(list(X = new_X, y = new_y))
}

#' Calculate AUC (Area Under Curve)
#' 
#' @param y_true True labels
#' @param y_scores Predicted scores
#' @return AUC value
calculate_auc <- function(y_true, y_scores) {
  tryCatch({
    # Simple AUC calculation
    n_pos <- sum(y_true == 1)
    n_neg <- sum(y_true == 0)
    
    if (n_pos == 0 || n_neg == 0) {
      return(0.5)
    }
    
    # Get all pairs
    pos_scores <- y_scores[y_true == 1]
    neg_scores <- y_scores[y_true == 0]
    
    # Count concordant pairs
    concordant <- 0
    for (pos_score in pos_scores) {
      concordant <- concordant + sum(pos_score > neg_scores)
    }
    
    auc <- concordant / (n_pos * n_neg)
    return(auc)
  }, error = function(e) {
    return(0.5)
  })
}
