#' AdaBoost Feature Selection
#' 
#' This function implements AdaBoost (Adaptive Boosting) for feature selection.
#' AdaBoost trains a sequence of weak classifiers, typically decision stumps,
#' with adaptive weighting to focus on misclassified examples.
#' 
#' @param X A matrix or data frame of features (predictors).
#' @param y A vector of target values (binary classification: 0 and 1, or factor/character).
#' @param num_features Number of features to select (default: 10).
#' @param n_estimators Number of weak classifiers (default: 50).
#' @param learning_rate Step size shrinkage used in updates (default: 1.0).
#' @param max_depth Maximum depth for decision stumps (default: 1).
#' @param min_samples_split Minimum samples required to split a node (default: 2).
#' @param min_samples_leaf Minimum samples required at leaf node (default: 1).
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
#'   \item{ensemble_weights}{Weights of each weak classifier}
#'   \item{cv_performance}{Cross-validation performance metrics}
#'   \item{feature_frequency}{Frequency of features used across weak classifiers}
#' 
#' @examples
#' \dontrun{
#' # Generate sample data
#' X <- matrix(rnorm(1000), ncol = 20)
#' y <- sample(0:1, 50, replace = TRUE)
#' 
#' # Apply AdaBoost feature selection
#' result <- fs_method_adaboost(X, y, num_features = 5)
#' print(result$selected_features)
#' }
#' 
#' @export
fs_method_adaboost <- function(X, y, num_features = 10, n_estimators = 50, 
                               learning_rate = 1.0, max_depth = 1,
                               min_samples_split = 2, min_samples_leaf = 1,
                               cv_folds = 5, random_state = 42, 
                               verbose = FALSE, use_smote = FALSE) {
  
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
  if (learning_rate <= 0) {
    stop("learning_rate must be positive")
  }
  
  # Convert inputs to standard format
  X <- as.matrix(X)
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  # Handle feature names
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("Feature_", 1:n_features)
  }
  feature_names <- colnames(X)
  
  # Convert target to binary
  if (is.factor(y) || is.character(y)) {
    unique_classes <- unique(y)
    if (length(unique_classes) != 2) {
      stop("AdaBoost requires exactly 2 classes")
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
      # Simple SMOTE implementation
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
  
  # Initialize variables for AdaBoost
  weak_classifiers <- list()
  classifier_weights <- numeric(n_estimators)
  feature_usage <- integer(n_features)
  feature_importance <- numeric(n_features)
  names(feature_importance) <- feature_names
  
  # Initialize sample weights (uniform)
  sample_weights <- rep(1/n_samples, n_samples)
  
  if (verbose) {
    cat("Training AdaBoost with", n_estimators, "estimators...\n")
  }
  
  # Train weak classifiers
  for (t in 1:n_estimators) {
    if (verbose && t %% 10 == 0) {
      cat("Training classifier", t, "of", n_estimators, "\n")
    }
    
    # Train weak classifier with current sample weights
    weak_classifier <- train_weak_decision_stump(X, y, sample_weights, 
                                                  max_depth, min_samples_split, 
                                                  min_samples_leaf)
    
    # Make predictions
    predictions <- predict_decision_stump(weak_classifier, X)
    
    # Calculate weighted error
    weighted_error <- sum(sample_weights * (predictions != y))
    
    # Avoid division by zero and ensure error < 0.5
    if (weighted_error <= 0) {
      weighted_error <- 1e-10
    }
    if (weighted_error >= 0.5) {
      # If error is too high, stop early
      if (verbose) {
        cat("Stopping early due to high error at iteration", t, "\n")
      }
      break
    }
    
    # Calculate classifier weight
    alpha <- learning_rate * 0.5 * log((1 - weighted_error) / weighted_error)
    
    # Update sample weights
    sample_weights <- sample_weights * exp(-alpha * y * predictions)
    sample_weights <- sample_weights / sum(sample_weights)  # Normalize
    
    # Store classifier and its weight
    weak_classifiers[[t]] <- weak_classifier
    classifier_weights[t] <- alpha
    
    # Update feature usage and importance
    feature_usage[weak_classifier$feature_idx] <- feature_usage[weak_classifier$feature_idx] + 1
    feature_importance[weak_classifier$feature_idx] <- feature_importance[weak_classifier$feature_idx] + alpha
  }
  
  # Remove unused classifier slots
  actual_estimators <- length(weak_classifiers)
  classifier_weights <- classifier_weights[1:actual_estimators]
  
  # Calculate final feature importance scores
  # Combine frequency-based and weight-based importance
  frequency_importance <- feature_usage / sum(feature_usage)
  weight_importance <- feature_importance / sum(feature_importance)
  
  # Final importance as weighted combination
  final_importance <- 0.7 * weight_importance + 0.3 * frequency_importance
  
  # Cross-validation performance assessment
  cv_results <- perform_adaboost_cv(original_X, original_y, cv_folds, 
                                    n_estimators, learning_rate, max_depth,
                                    min_samples_split, min_samples_leaf,
                                    use_smote, verbose)
  
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
      n_estimators = actual_estimators,
      learning_rate = learning_rate,
      max_depth = max_depth,
      min_samples_split = min_samples_split,
      min_samples_leaf = min_samples_leaf,
      cv_folds = cv_folds,
      random_state = random_state,
      use_smote = use_smote
    ),
    method_info = list(
      method_name = "AdaBoost",
      method_type = "Ensemble",
      description = "Adaptive Boosting with decision stumps",
      references = "Freund, Y. & Schapire, R. E. (1997). A decision-theoretic generalization of on-line learning."
    ),
    ensemble_weights = classifier_weights,
    cv_performance = cv_results,
    feature_frequency = feature_usage[top_indices],
    weak_classifiers = weak_classifiers
  )
  
  return(result)
}

#' Train Weak Decision Stump for AdaBoost
#' 
#' @param X Feature matrix
#' @param y Target vector
#' @param sample_weights Sample weights
#' @param max_depth Maximum depth (typically 1 for stumps)
#' @param min_samples_split Minimum samples to split
#' @param min_samples_leaf Minimum samples at leaf
#' @return Trained decision stump
train_weak_decision_stump <- function(X, y, sample_weights, max_depth = 1,
                                      min_samples_split = 2, min_samples_leaf = 1) {
  
  n_features <- ncol(X)
  
  best_feature <- 1
  best_threshold <- 0
  best_weighted_gini <- Inf
  best_left_prediction <- 0
  best_right_prediction <- 1
  
  # Try each feature
  for (feature_idx in 1:n_features) {
    feature_values <- X[, feature_idx]
    unique_values <- sort(unique(feature_values))
    
    # Try different thresholds
    for (i in 1:(length(unique_values) - 1)) {
      threshold <- (unique_values[i] + unique_values[i + 1]) / 2
      
      # Split samples
      left_mask <- feature_values <= threshold
      right_mask <- !left_mask
      
      # Check minimum samples constraints
      if (sum(left_mask) < min_samples_leaf || sum(right_mask) < min_samples_leaf) {
        next
      }
      
      # Calculate weighted Gini impurity
      left_weights <- sample_weights[left_mask]
      right_weights <- sample_weights[right_mask]
      left_y <- y[left_mask]
      right_y <- y[right_mask]
      
      if (length(left_y) == 0 || length(right_y) == 0) {
        next
      }
      
      # Calculate weighted Gini for left and right nodes
      left_gini <- calculate_weighted_gini(left_y, left_weights)
      right_gini <- calculate_weighted_gini(right_y, right_weights)
      
      # Weighted average of Gini impurities
      total_left_weight <- sum(left_weights)
      total_right_weight <- sum(right_weights)
      total_weight <- total_left_weight + total_right_weight
      
      weighted_gini <- (total_left_weight / total_weight) * left_gini +
                       (total_right_weight / total_weight) * right_gini
      
      # Update best split if this is better
      if (weighted_gini < best_weighted_gini) {
        best_weighted_gini <- weighted_gini
        best_feature <- feature_idx
        best_threshold <- threshold
        
        # Determine predictions for left and right nodes
        left_weighted_mean <- sum(left_y * left_weights) / sum(left_weights)
        right_weighted_mean <- sum(right_y * right_weights) / sum(right_weights)
        
        best_left_prediction <- ifelse(left_weighted_mean > 0.5, 1, -1)
        best_right_prediction <- ifelse(right_weighted_mean > 0.5, 1, -1)
      }
    }
  }
  
  # Return decision stump
  return(list(
    feature_idx = best_feature,
    threshold = best_threshold,
    left_prediction = best_left_prediction,
    right_prediction = best_right_prediction,
    weighted_gini = best_weighted_gini
  ))
}

#' Predict using Decision Stump
#' 
#' @param stump Decision stump model
#' @param X Feature matrix
#' @return Predictions vector
predict_decision_stump <- function(stump, X) {
  feature_values <- X[, stump$feature_idx]
  predictions <- ifelse(feature_values <= stump$threshold,
                        stump$left_prediction, stump$right_prediction)
  return(predictions)
}

#' Calculate Weighted Gini Impurity
#' 
#' @param y Target values
#' @param weights Sample weights
#' @return Weighted Gini impurity
calculate_weighted_gini <- function(y, weights) {
  if (length(y) == 0) return(0)
  
  # Calculate weighted class probabilities
  total_weight <- sum(weights)
  if (total_weight == 0) return(0)
  
  class_0_weight <- sum(weights[y == 0])
  class_1_weight <- sum(weights[y == 1])
  
  p0 <- class_0_weight / total_weight
  p1 <- class_1_weight / total_weight
  
  # Gini impurity
  gini <- 1 - p0^2 - p1^2
  return(gini)
}

#' Perform Cross-Validation for AdaBoost
#' 
#' @param X Feature matrix
#' @param y Target vector
#' @param cv_folds Number of CV folds
#' @param n_estimators Number of estimators
#' @param learning_rate Learning rate
#' @param max_depth Maximum depth
#' @param min_samples_split Minimum samples to split
#' @param min_samples_leaf Minimum samples at leaf
#' @param use_smote Whether to use SMOTE
#' @param verbose Whether to print progress
#' @return CV performance metrics
perform_adaboost_cv <- function(X, y, cv_folds, n_estimators, learning_rate,
                                max_depth, min_samples_split, min_samples_leaf,
                                use_smote, verbose) {
  
  n_samples <- nrow(X)
  fold_size <- floor(n_samples / cv_folds)
  indices <- sample(1:n_samples)
  
  cv_accuracies <- numeric(cv_folds)
  cv_aucs <- numeric(cv_folds)
  
  for (fold in 1:cv_folds) {
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
    
    # Train AdaBoost model
    model <- train_adaboost_cv_model(X_train, y_train, n_estimators, 
                                     learning_rate, max_depth, 
                                     min_samples_split, min_samples_leaf)
    
    # Make predictions
    predictions <- predict_adaboost_model(model, X_test)
    probabilities <- predict_adaboost_probabilities(model, X_test)
    
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

#' Train AdaBoost Model for CV
#' 
#' @param X Training features
#' @param y Training targets
#' @param n_estimators Number of estimators
#' @param learning_rate Learning rate
#' @param max_depth Maximum depth
#' @param min_samples_split Minimum samples to split
#' @param min_samples_leaf Minimum samples at leaf
#' @return Trained AdaBoost model
train_adaboost_cv_model <- function(X, y, n_estimators, learning_rate, max_depth,
                                    min_samples_split, min_samples_leaf) {
  
  n_samples <- nrow(X)
  sample_weights <- rep(1/n_samples, n_samples)
  
  weak_classifiers <- list()
  classifier_weights <- numeric(n_estimators)
  
  for (t in 1:n_estimators) {
    # Train weak classifier
    weak_classifier <- train_weak_decision_stump(X, y, sample_weights, 
                                                  max_depth, min_samples_split, 
                                                  min_samples_leaf)
    
    # Make predictions
    predictions <- predict_decision_stump(weak_classifier, X)
    
    # Calculate weighted error
    weighted_error <- sum(sample_weights * (predictions != y))
    
    if (weighted_error <= 0) {
      weighted_error <- 1e-10
    }
    if (weighted_error >= 0.5) {
      break
    }
    
    # Calculate classifier weight
    alpha <- learning_rate * 0.5 * log((1 - weighted_error) / weighted_error)
    
    # Update sample weights
    sample_weights <- sample_weights * exp(-alpha * y * predictions)
    sample_weights <- sample_weights / sum(sample_weights)
    
    # Store classifier
    weak_classifiers[[t]] <- weak_classifier
    classifier_weights[t] <- alpha
  }
  
  return(list(
    weak_classifiers = weak_classifiers,
    classifier_weights = classifier_weights[seq_along(weak_classifiers)]
  ))
}

#' Predict with AdaBoost Model
#' 
#' @param model Trained AdaBoost model
#' @param X Test features
#' @return Predictions
predict_adaboost_model <- function(model, X) {
  n_samples <- nrow(X)
  ensemble_predictions <- numeric(n_samples)
  
  for (i in seq_along(model$weak_classifiers)) {
    weak_predictions <- predict_decision_stump(model$weak_classifiers[[i]], X)
    ensemble_predictions <- ensemble_predictions + 
                           model$classifier_weights[i] * weak_predictions
  }
  
  final_predictions <- ifelse(ensemble_predictions > 0, 1, 0)
  return(final_predictions)
}

#' Predict Probabilities with AdaBoost Model
#' 
#' @param model Trained AdaBoost model
#' @param X Test features
#' @return Prediction probabilities
predict_adaboost_probabilities <- function(model, X) {
  n_samples <- nrow(X)
  ensemble_scores <- numeric(n_samples)
  
  for (i in seq_along(model$weak_classifiers)) {
    weak_predictions <- predict_decision_stump(model$weak_classifiers[[i]], X)
    ensemble_scores <- ensemble_scores + 
                      model$classifier_weights[i] * weak_predictions
  }
  
  # Convert to probabilities using sigmoid
  probabilities <- 1 / (1 + exp(-ensemble_scores))
  return(probabilities)
}

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
  
  for (i in 1:n_synthetic) {
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
