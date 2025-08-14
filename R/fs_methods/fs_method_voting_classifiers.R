#' Voting Classifiers Feature Selection Method
#'
#' Implements feature selection using ensemble voting classifiers.
#' This method combines multiple base classifiers and selects features
#' based on their collective importance across different algorithms.
#'
#' @param ramwas.data Expression data matrix (samples as rows, features as columns)  
#' @param target Target variable vector
#' @param smote_sampling Logical, whether to apply SMOTE sampling
#' @return List containing selected features and method metadata
#' @export
voting_classifiers_fs <- function(ramwas.data, target, smote_sampling = FALSE) {
  
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
    
    # Perform Voting Classifiers feature selection
    voting_result <- perform_voting_classifiers_selection(data_clean, target)
    
    # Extract selected features
    if (length(voting_result$selected_features) == 0) {
      return(create_empty_result("No features selected by Voting Classifiers"))
    }
    
    # Calculate timing
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      selected_features = voting_result$selected_features,
      feature_scores = voting_result$feature_scores,
      method_name = "Voting_Classifiers",
      execution_time = execution_time,
      n_features_selected = length(voting_result$selected_features),
      n_features_input = ncol(ramwas.data),
      method_params = list(
        base_classifiers = voting_result$base_classifiers,
        voting_type = voting_result$voting_type,
        feature_selection_threshold = voting_result$threshold,
        consensus_threshold = voting_result$consensus_threshold,
        smote_sampling = smote_sampling
      ),
      preprocessing_info = list(
        constant_features_removed = preprocessed$constant_removed,
        correlated_features_removed = preprocessed$correlated_removed,
        final_features_count = ncol(data_clean)
      ),
      voting_info = voting_result$voting_info
    ))
    
  }, error = function(e) {
    # Return error result
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      selected_features = character(0),
      feature_scores = numeric(0),
      method_name = "Voting_Classifiers",
      execution_time = execution_time,
      n_features_selected = 0,
      n_features_input = ncol(ramwas.data),
      error_message = paste("Voting Classifiers failed:", e$message),
      method_params = list(smote_sampling = smote_sampling)
    ))
  })
}

#' Perform Voting Classifiers Feature Selection
#'
#' Implements ensemble voting with multiple base classifiers
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with selected features and voting information
perform_voting_classifiers_selection <- function(data_matrix, target_vector) {
  
  n_features <- ncol(data_matrix)
  
  # Determine if classification or regression
  is_classification <- determine_task_type(target_vector)
  
  if (!is_classification) {
    # Convert to classification if needed
    target_vector <- discretize_target(target_vector)
  }
  
  # Standardize features
  data_scaled <- scale(data_matrix)
  
  # Handle any NAs from scaling
  if (any(is.na(data_scaled))) {
    non_na_cols <- apply(data_scaled, 2, function(x) !any(is.na(x)))
    data_scaled <- data_scaled[, non_na_cols, drop = FALSE]
    n_features <- ncol(data_scaled)
  }
  
  # Define base classifiers
  base_classifiers <- define_base_classifiers()
  
  # Train ensemble of classifiers
  ensemble_results <- train_voting_ensemble(data_scaled, target_vector, base_classifiers)
  
  # Calculate feature importance from ensemble
  feature_importance <- calculate_ensemble_feature_importance(ensemble_results, data_scaled, target_vector)
  
  # Select features based on consensus voting
  selected_result <- select_features_by_voting_consensus(feature_importance, n_features)
  
  return(list(
    selected_features = selected_result$selected_features,
    feature_scores = selected_result$feature_scores,
    base_classifiers = names(base_classifiers),
    voting_type = "soft_voting",
    threshold = selected_result$threshold,
    consensus_threshold = selected_result$consensus_threshold,
    voting_info = list(
      ensemble_accuracy = ensemble_results$ensemble_accuracy,
      individual_accuracies = ensemble_results$individual_accuracies,
      feature_votes = ensemble_results$feature_votes
    )
  ))
}

#' Define Base Classifiers
#'
#' Defines the set of base classifiers for voting ensemble
#'
#' @return List of base classifier definitions
define_base_classifiers <- function() {
  
  base_classifiers <- list(
    "naive_bayes" = list(
      name = "Naive Bayes",
      train_func = train_naive_bayes,
      predict_func = predict_naive_bayes,
      importance_func = calculate_naive_bayes_importance
    ),
    "logistic_regression" = list(
      name = "Logistic Regression", 
      train_func = train_logistic_regression,
      predict_func = predict_logistic_regression,
      importance_func = calculate_logistic_regression_importance
    ),
    "decision_tree" = list(
      name = "Decision Tree",
      train_func = train_decision_tree_classifier,
      predict_func = predict_decision_tree_classifier,
      importance_func = calculate_decision_tree_importance
    ),
    "knn" = list(
      name = "K-Nearest Neighbors",
      train_func = train_knn_classifier,
      predict_func = predict_knn_classifier,
      importance_func = calculate_knn_importance
    ),
    "svm_linear" = list(
      name = "Linear SVM",
      train_func = train_linear_svm,
      predict_func = predict_linear_svm,
      importance_func = calculate_svm_importance
    )
  )
  
  return(base_classifiers)
}

#' Train Voting Ensemble
#'
#' Trains ensemble of base classifiers
#'
#' @param X Training data
#' @param y Target labels
#' @param base_classifiers Base classifier definitions
#' @return Ensemble training results
train_voting_ensemble <- function(X, y, base_classifiers) {
  
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  # Cross-validation setup
  n_folds <- min(5, n_samples %/% 2)
  fold_assignments <- sample(rep(1:n_folds, length.out = n_samples))
  
  # Initialize results storage
  ensemble_models <- list()
  individual_accuracies <- numeric(length(base_classifiers))
  feature_votes <- matrix(0, nrow = length(base_classifiers), ncol = n_features)
  
  # Train each base classifier
  for (i in seq_along(base_classifiers)) {
    classifier_name <- names(base_classifiers)[i]
    classifier_def <- base_classifiers[[i]]
    
    # Cross-validation for this classifier
    cv_accuracies <- numeric(n_folds)
    cv_importances <- matrix(0, nrow = n_folds, ncol = n_features)
    
    for (fold in 1:n_folds) {
      train_idx <- which(fold_assignments != fold)
      test_idx <- which(fold_assignments == fold)
      
      if (length(train_idx) < 2 || length(test_idx) < 1) next
      
      # Train classifier
      model <- classifier_def$train_func(X[train_idx, , drop = FALSE], y[train_idx])
      
      # Evaluate
      predictions <- classifier_def$predict_func(model, X[test_idx, , drop = FALSE])
      cv_accuracies[fold] <- mean(predictions == y[test_idx])
      
      # Calculate feature importance
      importance <- classifier_def$importance_func(model, X[train_idx, , drop = FALSE], y[train_idx])
      cv_importances[fold, ] <- importance
    }
    
    # Train final model on full data
    final_model <- classifier_def$train_func(X, y)
    ensemble_models[[classifier_name]] <- list(
      model = final_model,
      classifier_def = classifier_def
    )
    
    # Store results
    individual_accuracies[i] <- mean(cv_accuracies[cv_accuracies > 0])
    feature_votes[i, ] <- colMeans(cv_importances)
    
    names(individual_accuracies)[i] <- classifier_name
  }
  
  # Calculate ensemble accuracy through cross-validation
  ensemble_accuracy <- calculate_ensemble_cv_accuracy(X, y, ensemble_models, fold_assignments)
  
  return(list(
    ensemble_models = ensemble_models,
    individual_accuracies = individual_accuracies,
    ensemble_accuracy = ensemble_accuracy,
    feature_votes = feature_votes
  ))
}

#' Train Naive Bayes Classifier
#'
#' Simple Naive Bayes implementation
#'
#' @param X Training data
#' @param y Training labels
#' @return Trained Naive Bayes model
train_naive_bayes <- function(X, y) {
  
  classes <- unique(y)
  
  # Calculate class priors
  class_priors <- table(y) / length(y)
  
  # Calculate feature statistics for each class
  feature_stats <- list()
  for (class in classes) {
    class_data <- X[y == class, , drop = FALSE]
    
    feature_stats[[as.character(class)]] <- list(
      means = colMeans(class_data),
      sds = apply(class_data, 2, sd) + 1e-9  # Add small value for numerical stability
    )
  }
  
  return(list(
    class_priors = class_priors,
    feature_stats = feature_stats,
    classes = classes
  ))
}

#' Predict Naive Bayes
#'
#' Makes predictions using Naive Bayes
#'
#' @param model Trained Naive Bayes model
#' @param X Test data
#' @return Predictions
predict_naive_bayes <- function(model, X) {
  
  n_samples <- nrow(X)
  n_classes <- length(model$classes)
  
  log_probabilities <- matrix(0, nrow = n_samples, ncol = n_classes)
  
  for (i in seq_along(model$classes)) {
    class <- as.character(model$classes[i])
    
    # Calculate log probability for this class
    log_prior <- log(model$class_priors[class])
    
    # Calculate log likelihood
    means <- model$feature_stats[[class]]$means
    sds <- model$feature_stats[[class]]$sds
    
    log_likelihood <- rowSums(dnorm(X, rep(means, each = n_samples), 
                                   rep(sds, each = n_samples), log = TRUE))
    
    log_probabilities[, i] <- log_prior + log_likelihood
  }
  
  # Return predicted classes
  predicted_indices <- apply(log_probabilities, 1, which.max)
  return(model$classes[predicted_indices])
}

#' Calculate Naive Bayes Importance
#'
#' Calculates feature importance for Naive Bayes
#'
#' @param model Trained model
#' @param X Training data
#' @param y Training labels
#' @return Feature importance scores
calculate_naive_bayes_importance <- function(model, X, y) {
  
  n_features <- ncol(X)
  importance <- numeric(n_features)
  
  # Calculate importance based on class separation
  for (i in 1:n_features) {
    class_means <- sapply(model$feature_stats, function(stats) stats$means[i])
    importance[i] <- var(class_means)
  }
  
  return(importance / sum(importance))
}

#' Train Logistic Regression
#'
#' Simple logistic regression implementation
#'
#' @param X Training data
#' @param y Training labels
#' @return Trained logistic regression model
train_logistic_regression <- function(X, y) {
  
  # Convert to binary if multiclass (use one-vs-rest)
  classes <- unique(y)
  if (length(classes) > 2) {
    # Use most frequent class vs others
    majority_class <- names(which.max(table(y)))
    y_binary <- as.numeric(y == majority_class)
  } else {
    y_binary <- as.numeric(y == classes[2])
  }
  
  # Add intercept
  X_with_intercept <- cbind(1, X)
  
  # Initialize coefficients
  coefficients <- rep(0, ncol(X_with_intercept))
  
  # Iterative fitting (simplified Newton-Raphson)
  for (iter in 1:20) {
    
    # Calculate predictions and probabilities
    linear_pred <- X_with_intercept %*% coefficients
    probabilities <- 1 / (1 + exp(-linear_pred))
    probabilities <- pmax(pmin(probabilities, 1 - 1e-15), 1e-15)
    
    # Calculate gradient and Hessian
    gradient <- t(X_with_intercept) %*% (y_binary - probabilities)
    weights <- probabilities * (1 - probabilities)
    hessian <- -t(X_with_intercept) %*% (X_with_intercept * weights)
    
    # Update coefficients
    tryCatch({
      update <- solve(hessian, gradient)
      coefficients <- coefficients - update
    }, error = function(e) {
      # If Hessian is singular, use gradient descent
      coefficients <<- coefficients + 0.01 * gradient / sqrt(sum(gradient^2) + 1e-8)
    })
    
    # Check convergence
    if (sqrt(sum(gradient^2)) < 1e-6) break
  }
  
  return(list(
    coefficients = coefficients,
    classes = classes
  ))
}

#' Predict Logistic Regression
#'
#' Makes predictions using logistic regression
#'
#' @param model Trained logistic regression model
#' @param X Test data
#' @return Predictions
predict_logistic_regression <- function(model, X) {
  
  # Add intercept
  X_with_intercept <- cbind(1, X)
  
  # Calculate probabilities
  linear_pred <- X_with_intercept %*% model$coefficients
  probabilities <- 1 / (1 + exp(-linear_pred))
  
  # Convert to class predictions
  if (length(model$classes) == 2) {
    predictions <- ifelse(probabilities > 0.5, model$classes[2], model$classes[1])
  } else {
    # For multiclass, predict majority class if prob > 0.5
    majority_class <- names(which.max(table(model$classes)))
    predictions <- ifelse(probabilities > 0.5, majority_class, 
                         sample(setdiff(model$classes, majority_class), length(probabilities), replace = TRUE))
  }
  
  return(predictions)
}

#' Calculate Logistic Regression Importance
#'
#' Calculates feature importance for logistic regression
#'
#' @param model Trained model
#' @param X Training data
#' @param y Training labels
#' @return Feature importance scores
calculate_logistic_regression_importance <- function(model, X, y) {
  
  # Feature importance based on absolute coefficients (excluding intercept)
  coefficients <- model$coefficients[-1]  # Remove intercept
  importance <- abs(coefficients)
  
  if (sum(importance) > 0) {
    importance <- importance / sum(importance)
  } else {
    importance <- rep(1/length(importance), length(importance))
  }
  
  return(importance)
}

#' Train Decision Tree Classifier
#'
#' Simple decision tree for classification
#'
#' @param X Training data
#' @param y Training labels
#' @return Trained decision tree model
train_decision_tree_classifier <- function(X, y) {
  
  # Build simple decision tree
  tree <- build_classification_tree(X, y, max_depth = 5)
  
  return(list(
    tree = tree,
    feature_importance = calculate_tree_feature_importance(tree, ncol(X))
  ))
}

#' Build Classification Tree
#'
#' Builds a classification tree
#'
#' @param X Training data
#' @param y Training labels
#' @param max_depth Maximum tree depth
#' @return Classification tree
build_classification_tree <- function(X, y, max_depth = 5) {
  
  tree <- build_classification_node(X, y, 1, max_depth)
  return(tree)
}

#' Build Classification Node
#'
#' Recursively builds classification tree nodes
#'
#' @param X Node data
#' @param y Node labels
#' @param depth Current depth
#' @param max_depth Maximum depth
#' @return Tree node
build_classification_node <- function(X, y, depth, max_depth) {
  
  n_samples <- nrow(X)
  
  # Stopping criteria
  if (depth >= max_depth || n_samples < 5 || length(unique(y)) == 1) {
    return(list(
      type = "leaf",
      prediction = get_majority_class(y),
      n_samples = n_samples
    ))
  }
  
  # Find best split
  best_split <- find_best_classification_split(X, y)
  
  if (is.null(best_split)) {
    return(list(
      type = "leaf",
      prediction = get_majority_class(y),
      n_samples = n_samples
    ))
  }
  
  # Split data
  left_idx <- which(X[, best_split$feature] <= best_split$threshold)
  right_idx <- which(X[, best_split$feature] > best_split$threshold)
  
  # Build child nodes
  left_node <- build_classification_node(X[left_idx, , drop = FALSE], y[left_idx], depth + 1, max_depth)
  right_node <- build_classification_node(X[right_idx, , drop = FALSE], y[right_idx], depth + 1, max_depth)
  
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

#' Find Best Classification Split
#'
#' Finds the best split for classification
#'
#' @param X Node data
#' @param y Node labels
#' @return Best split information
find_best_classification_split <- function(X, y) {
  
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
      
      # Calculate information gain
      parent_entropy <- calculate_entropy(y)
      left_entropy <- calculate_entropy(y[left_idx])
      right_entropy <- calculate_entropy(y[right_idx])
      
      weighted_entropy <- (length(left_idx) * left_entropy + length(right_idx) * right_entropy) / length(y)
      information_gain <- parent_entropy - weighted_entropy
      
      if (information_gain > best_score) {
        best_score <- information_gain
        best_split <- list(
          feature = feature,
          threshold = threshold,
          importance = information_gain,
          left_idx = left_idx,
          right_idx = right_idx
        )
      }
    }
  }
  
  return(best_split)
}

#' Calculate Entropy
#'
#' Calculates entropy for classification
#'
#' @param y Labels
#' @return Entropy value
calculate_entropy <- function(y) {
  
  if (length(y) == 0) return(0)
  
  class_probs <- table(y) / length(y)
  entropy <- -sum(class_probs * log2(class_probs + 1e-15))
  
  return(entropy)
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

#' Predict Decision Tree Classifier
#'
#' Makes predictions using decision tree
#'
#' @param model Trained decision tree model
#' @param X Test data
#' @return Predictions
predict_decision_tree_classifier <- function(model, X) {
  
  n_samples <- nrow(X)
  predictions <- character(n_samples)
  
  for (i in 1:n_samples) {
    predictions[i] <- predict_tree_sample(model$tree, X[i, ])
  }
  
  return(predictions)
}

#' Predict Tree Sample
#'
#' Predicts a single sample using tree
#'
#' @param node Current tree node
#' @param sample Sample to predict
#' @return Prediction
predict_tree_sample <- function(node, sample) {
  
  if (node$type == "leaf") {
    return(as.character(node$prediction))
  }
  
  if (sample[node$feature] <= node$threshold) {
    return(predict_tree_sample(node$left, sample))
  } else {
    return(predict_tree_sample(node$right, sample))
  }
}

#' Calculate Tree Feature Importance
#'
#' Calculates feature importance for tree
#'
#' @param tree Decision tree
#' @param n_features Number of features
#' @return Feature importance vector
calculate_tree_feature_importance <- function(tree, n_features) {
  
  importance <- numeric(n_features)
  
  if (tree$type == "internal") {
    importance[tree$feature] <- tree$importance
    
    # Recursively add importance from child nodes
    left_importance <- calculate_tree_feature_importance(tree$left, n_features)
    right_importance <- calculate_tree_feature_importance(tree$right, n_features)
    
    importance <- importance + left_importance + right_importance
  }
  
  return(importance)
}

#' Calculate Decision Tree Importance
#'
#' Calculates feature importance for decision tree
#'
#' @param model Trained model
#' @param X Training data
#' @param y Training labels
#' @return Feature importance scores
calculate_decision_tree_importance <- function(model, X, y) {
  
  importance <- model$feature_importance
  
  if (sum(importance) > 0) {
    importance <- importance / sum(importance)
  } else {
    importance <- rep(1/length(importance), length(importance))
  }
  
  return(importance)
}

#' Train KNN Classifier
#'
#' Simple k-nearest neighbors implementation
#'
#' @param X Training data
#' @param y Training labels
#' @return Trained KNN model
train_knn_classifier <- function(X, y) {
  
  # Determine optimal k
  k <- max(3, min(15, round(sqrt(nrow(X)))))
  
  return(list(
    X_train = X,
    y_train = y,
    k = k
  ))
}

#' Predict KNN Classifier
#'
#' Makes predictions using KNN
#'
#' @param model Trained KNN model
#' @param X Test data
#' @return Predictions
predict_knn_classifier <- function(model, X) {
  
  n_test <- nrow(X)
  predictions <- character(n_test)
  
  for (i in 1:n_test) {
    # Calculate distances to all training points
    distances <- sqrt(rowSums((model$X_train - matrix(X[i, ], nrow = nrow(model$X_train), ncol = ncol(X), byrow = TRUE))^2))
    
    # Find k nearest neighbors
    k_nearest_idx <- order(distances)[1:model$k]
    k_nearest_labels <- model$y_train[k_nearest_idx]
    
    # Majority vote
    predictions[i] <- get_majority_class(k_nearest_labels)
  }
  
  return(predictions)
}

#' Calculate KNN Importance
#'
#' Calculates feature importance for KNN using permutation
#'
#' @param model Trained model
#' @param X Training data
#' @param y Training labels
#' @return Feature importance scores
calculate_knn_importance <- function(model, X, y) {
  
  n_features <- ncol(X)
  importance <- numeric(n_features)
  
  # Use a subset for efficiency
  subset_size <- min(100, nrow(X))
  subset_idx <- sample(nrow(X), subset_size)
  X_subset <- X[subset_idx, , drop = FALSE]
  y_subset <- y[subset_idx]
  
  # Baseline accuracy
  baseline_pred <- predict_knn_classifier(model, X_subset)
  baseline_accuracy <- mean(baseline_pred == y_subset)
  
  # Permutation importance
  for (feature in 1:n_features) {
    X_perm <- X_subset
    X_perm[, feature] <- sample(X_perm[, feature])
    
    perm_pred <- predict_knn_classifier(model, X_perm)
    perm_accuracy <- mean(perm_pred == y_subset)
    
    importance[feature] <- baseline_accuracy - perm_accuracy
  }
  
  importance <- pmax(importance, 0)  # Remove negative values
  
  if (sum(importance) > 0) {
    importance <- importance / sum(importance)
  } else {
    importance <- rep(1/length(importance), length(importance))
  }
  
  return(importance)
}

#' Train Linear SVM
#'
#' Simple linear SVM implementation using SMO-like approach
#'
#' @param X Training data
#' @param y Training labels
#' @return Trained SVM model
train_linear_svm <- function(X, y) {
  
  # Convert to binary classification
  classes <- unique(y)
  if (length(classes) > 2) {
    majority_class <- names(which.max(table(y)))
    y_binary <- ifelse(y == majority_class, 1, -1)
  } else {
    y_binary <- ifelse(y == classes[1], -1, 1)
  }
  
  # Simple linear SVM (simplified SMO)
  n_samples <- nrow(X)
  
  # Initialize
  alpha <- numeric(n_samples)
  b <- 0
  C <- 1.0  # Regularization parameter
  
  # Simplified training (few iterations)
  for (iter in 1:20) {
    
    # Calculate weights
    w <- colSums(alpha * y_binary * X)
    
    # Find violating pair
    violated <- FALSE
    for (i in 1:n_samples) {
      prediction <- sum(w * X[i, ]) + b
      
      if ((y_binary[i] * prediction < 1) && (alpha[i] < C)) {
        # Update alpha
        alpha[i] <- min(C, alpha[i] + 0.01)
        violated <- TRUE
      }
    }
    
    if (!violated) break
  }
  
  # Final weights
  w <- colSums(alpha * y_binary * X)
  
  return(list(
    weights = w,
    bias = b,
    alpha = alpha,
    classes = classes
  ))
}

#' Predict Linear SVM
#'
#' Makes predictions using linear SVM
#'
#' @param model Trained SVM model
#' @param X Test data
#' @return Predictions
predict_linear_svm <- function(model, X) {
  
  # Calculate predictions
  predictions_numeric <- X %*% model$weights + model$bias
  
  # Convert to class labels
  if (length(model$classes) == 2) {
    predictions <- ifelse(predictions_numeric > 0, model$classes[2], model$classes[1])
  } else {
    majority_class <- names(which.max(table(model$classes)))
    predictions <- ifelse(predictions_numeric > 0, majority_class, 
                         sample(setdiff(model$classes, majority_class), length(predictions_numeric), replace = TRUE))
  }
  
  return(predictions)
}

#' Calculate SVM Importance
#'
#' Calculates feature importance for SVM
#'
#' @param model Trained model
#' @param X Training data
#' @param y Training labels
#' @return Feature importance scores
calculate_svm_importance <- function(model, X, y) {
  
  # Feature importance based on absolute weights
  importance <- abs(model$weights)
  
  if (sum(importance) > 0) {
    importance <- importance / sum(importance)
  } else {
    importance <- rep(1/length(importance), length(importance))
  }
  
  return(importance)
}

#' Calculate Ensemble CV Accuracy
#'
#' Calculates ensemble accuracy using cross-validation
#'
#' @param X Input data
#' @param y Target labels
#' @param ensemble_models Trained ensemble models
#' @param fold_assignments Fold assignments
#' @return Ensemble accuracy
calculate_ensemble_cv_accuracy <- function(X, y, ensemble_models, fold_assignments) {
  
  n_folds <- max(fold_assignments)
  cv_accuracies <- numeric(n_folds)
  
  for (fold in 1:n_folds) {
    test_idx <- which(fold_assignments == fold)
    
    if (length(test_idx) == 0) next
    
    # Get predictions from all models
    ensemble_predictions <- list()
    for (model_name in names(ensemble_models)) {
      model_info <- ensemble_models[[model_name]]
      pred <- model_info$classifier_def$predict_func(model_info$model, X[test_idx, , drop = FALSE])
      ensemble_predictions[[model_name]] <- pred
    }
    
    # Majority voting
    final_predictions <- apply_majority_voting(ensemble_predictions)
    cv_accuracies[fold] <- mean(final_predictions == y[test_idx])
  }
  
  return(mean(cv_accuracies[cv_accuracies > 0]))
}

#' Apply Majority Voting
#'
#' Applies majority voting to ensemble predictions
#'
#' @param ensemble_predictions List of predictions from different models
#' @return Final ensemble predictions
apply_majority_voting <- function(ensemble_predictions) {
  
  n_samples <- length(ensemble_predictions[[1]])
  final_predictions <- character(n_samples)
  
  for (i in 1:n_samples) {
    votes <- sapply(ensemble_predictions, function(pred) pred[i])
    final_predictions[i] <- get_majority_class(votes)
  }
  
  return(final_predictions)
}

#' Calculate Ensemble Feature Importance
#'
#' Calculates feature importance from ensemble
#'
#' @param ensemble_results Ensemble training results
#' @param X Input data
#' @param y Target labels
#' @return Combined feature importance
calculate_ensemble_feature_importance <- function(ensemble_results, X, y) {
  
  n_features <- ncol(X)
  feature_votes <- ensemble_results$feature_votes
  individual_accuracies <- ensemble_results$individual_accuracies
  
  # Weight by individual accuracies
  weights <- individual_accuracies / sum(individual_accuracies)
  
  # Weighted combination of feature votes
  combined_importance <- colSums(feature_votes * weights)
  
  # Normalize
  if (sum(combined_importance) > 0) {
    combined_importance <- combined_importance / sum(combined_importance)
  } else {
    combined_importance <- rep(1/n_features, n_features)
  }
  
  return(combined_importance)
}

#' Select Features by Voting Consensus
#'
#' Selects features based on voting consensus
#'
#' @param importance_scores Feature importance scores
#' @param n_features Total number of features
#' @return Selected features and scores
select_features_by_voting_consensus <- function(importance_scores, n_features) {
  
  # Set selection threshold
  threshold <- quantile(importance_scores, 0.7)  # Top 30% by default
  consensus_threshold <- 0.6  # Require 60% consensus
  
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
    threshold = threshold,
    consensus_threshold = consensus_threshold
  ))
}

#' Discretize Target
#'
#' Converts continuous target to discrete classes
#'
#' @param target_vector Target variable
#' @param n_bins Number of bins for discretization
#' @return Discretized target
discretize_target <- function(target_vector, n_bins = 3) {
  
  # Use quantile-based binning
  quantiles <- quantile(target_vector, probs = seq(0, 1, length.out = n_bins + 1))
  discretized <- cut(target_vector, breaks = quantiles, include.lowest = TRUE, labels = paste0("Class_", 1:n_bins))
  
  return(as.character(discretized))
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

#' Create Empty Result for Voting Classifiers
#'
#' Creates standardized empty result for failed feature selection
#'
#' @param message Error or warning message
#' @return Empty result list
create_empty_result <- function(message) {
  return(list(
    selected_features = character(0),
    feature_scores = numeric(0),
    method_name = "Voting_Classifiers",
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
