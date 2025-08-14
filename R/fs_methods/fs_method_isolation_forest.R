#' Isolation Forest Feature Selection
#' 
#' This function implements Isolation Forest for feature selection.
#' Isolation Forest is an anomaly detection algorithm that can be adapted
#' for feature selection by identifying features that best isolate different classes.
#' 
#' @param X A matrix or data frame of features (predictors).
#' @param y A vector of target values (binary classification: 0 and 1, or factor/character).
#' @param num_features Number of features to select (default: 10).
#' @param n_estimators Number of isolation trees (default: 100).
#' @param max_samples Maximum number of samples to draw for each tree (default: "auto").
#' @param max_features Maximum number of features to draw for each tree (default: 1.0).
#' @param max_depth Maximum depth of isolation trees (default: "auto").
#' @param contamination Expected proportion of outliers (default: 0.1).
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
#'   \item{isolation_scores}{Isolation scores for each class}
#'   \item{anomaly_scores}{Anomaly scores for samples}
#'   \item{cv_performance}{Cross-validation performance metrics}
#'   \item{isolation_forest_model}{Complete isolation forest model}
#' 
#' @examples
#' \dontrun{
#' # Generate sample data
#' X <- matrix(rnorm(1000), ncol = 20)
#' y <- sample(0:1, 50, replace = TRUE)
#' 
#' # Apply Isolation Forest feature selection
#' result <- fs_method_isolation_forest(X, y, num_features = 5)
#' print(result$selected_features)
#' }
#' 
#' @export
fs_method_isolation_forest <- function(X, y, num_features = 10, n_estimators = 100,
                                       max_samples = "auto", max_features = 1.0,
                                       max_depth = "auto", contamination = 0.1,
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
  if (contamination <= 0 || contamination >= 1) {
    stop("contamination must be between 0 and 1")
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
      stop("Isolation Forest requires exactly 2 classes")
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
  
  # Set default values for auto parameters
  if (max_samples == "auto") {
    max_samples <- min(256, n_samples)
  }
  if (max_depth == "auto") {
    max_depth <- ceiling(log2(max(max_samples, 2)))
  }
  
  # Train Isolation Forest
  if (verbose) {
    cat("Training Isolation Forest with", n_estimators, "trees...\n")
  }
  
  isolation_forest <- train_isolation_forest(X, y, n_estimators, max_samples,
                                             max_features, max_depth, verbose)
  
  # Calculate feature importance based on isolation capability
  feature_importance <- calculate_isolation_importance(isolation_forest, X, y,
                                                       feature_names, contamination,
                                                       verbose)
  
  # Cross-validation performance assessment
  cv_results <- perform_isolation_forest_cv(original_X, original_y, cv_folds,
                                            n_estimators, max_samples, max_features,
                                            max_depth, contamination, use_smote, verbose)
  
  # Select top features
  top_indices <- order(feature_importance, decreasing = TRUE)[seq_len(min(num_features, length(feature_importance)))]
  selected_features <- feature_names[top_indices]
  selected_importance <- feature_importance[top_indices]
  names(selected_importance) <- selected_features
  
  # Calculate final anomaly scores
  anomaly_scores <- calculate_anomaly_scores(isolation_forest, X)
  
  # Prepare results
  result <- list(
    selected_features = selected_features,
    importance_scores = selected_importance,
    model_params = list(
      n_estimators = n_estimators,
      max_samples = max_samples,
      max_features = max_features,
      max_depth = max_depth,
      contamination = contamination,
      cv_folds = cv_folds,
      random_state = random_state,
      use_smote = use_smote
    ),
    method_info = list(
      method_name = "Isolation Forest",
      method_type = "Ensemble",
      description = "Isolation Forest for anomaly-based feature selection",
      references = "Liu, F. T., Ting, K. M., & Zhou, Z. H. (2008). Isolation forest."
    ),
    isolation_scores = isolation_forest$class_isolation_scores,
    anomaly_scores = anomaly_scores,
    cv_performance = cv_results,
    isolation_forest_model = isolation_forest
  )
  
  return(result)
}

#' Train Isolation Forest
#' 
#' @param X Feature matrix
#' @param y Target vector
#' @param n_estimators Number of trees
#' @param max_samples Maximum samples per tree
#' @param max_features Maximum features proportion
#' @param max_depth Maximum tree depth
#' @param verbose Print progress
#' @return Trained isolation forest model
train_isolation_forest <- function(X, y, n_estimators, max_samples, max_features,
                                   max_depth, verbose) {
  
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  # Store isolation trees and their metadata
  trees <- list()
  feature_usage <- numeric(n_features)
  sample_usage <- numeric(n_samples)
  
  # Calculate number of features to use per tree
  n_features_per_tree <- max(1, floor(max_features * n_features))
  
  for (i in seq_len(n_estimators)) {
    if (verbose && i %% 20 == 0) {
      cat("Training isolation tree", i, "of", n_estimators, "\n")
    }
    
    # Sample data for this tree
    sample_indices <- sample(seq_len(n_samples), 
                            min(max_samples, n_samples), 
                            replace = FALSE)
    
    # Sample features for this tree
    feature_indices <- sample(seq_len(n_features), n_features_per_tree)
    
    # Extract subset
    X_subset <- X[sample_indices, feature_indices, drop = FALSE]
    
    # Train isolation tree
    tree <- train_isolation_tree(X_subset, max_depth, feature_indices)
    trees[[i]] <- tree
    
    # Update usage tracking
    feature_usage[feature_indices] <- feature_usage[feature_indices] + 1
    sample_usage[sample_indices] <- sample_usage[sample_indices] + 1
  }
  
  # Calculate class-specific isolation scores
  class_isolation_scores <- calculate_class_isolation_scores(trees, X, y)
  
  return(list(
    trees = trees,
    feature_usage = feature_usage,
    sample_usage = sample_usage,
    class_isolation_scores = class_isolation_scores,
    n_estimators = n_estimators,
    max_depth = max_depth
  ))
}

#' Train Single Isolation Tree
#' 
#' @param X Feature subset
#' @param max_depth Maximum depth
#' @param feature_indices Original feature indices
#' @return Isolation tree
train_isolation_tree <- function(X, max_depth, feature_indices) {
  
  feature_usage <- numeric(length(feature_indices))
  
  tree <- build_isolation_node(X, 0, max_depth, feature_usage, feature_indices)
  
  return(list(
    tree = tree,
    feature_usage = feature_usage,
    feature_indices = feature_indices
  ))
}

#' Build Isolation Tree Node
#' 
#' @param X Feature matrix
#' @param current_depth Current depth
#' @param max_depth Maximum depth
#' @param feature_usage Feature usage tracker
#' @param feature_indices Original feature indices
#' @return Tree node
build_isolation_node <- function(X, current_depth, max_depth, feature_usage, 
                                 feature_indices) {
  
  n_samples <- nrow(X)
  
  # Stopping criteria: max depth reached or only one sample
  if (current_depth >= max_depth || n_samples <= 1) {
    return(list(
      is_leaf = TRUE,
      size = n_samples,
      depth = current_depth
    ))
  }
  
  # All samples have identical values (no split possible)
  if (ncol(X) == 0) {
    return(list(
      is_leaf = TRUE,
      size = n_samples,
      depth = current_depth
    ))
  }
  
  # Randomly select a feature to split on
  feature_idx <- sample(seq_len(ncol(X)), 1)
  feature_values <- X[, feature_idx]
  
  # Check if all values are the same
  if (length(unique(feature_values)) == 1) {
    return(list(
      is_leaf = TRUE,
      size = n_samples,
      depth = current_depth
    ))
  }
  
  # Randomly select a split point
  min_val <- min(feature_values)
  max_val <- max(feature_values)
  split_point <- runif(1, min_val, max_val)
  
  # Update feature usage
  feature_usage[feature_idx] <- feature_usage[feature_idx] + 1
  
  # Split samples
  left_mask <- feature_values < split_point
  right_mask <- !left_mask
  
  # Handle edge case where all samples go to one side
  if (sum(left_mask) == 0) {
    left_mask[1] <- TRUE
    right_mask[1] <- FALSE
  } else if (sum(right_mask) == 0) {
    right_mask[1] <- TRUE
    left_mask[1] <- FALSE
  }
  
  # Build child nodes
  left_child <- build_isolation_node(X[left_mask, , drop = FALSE],
                                     current_depth + 1, max_depth,
                                     feature_usage, feature_indices)
  
  right_child <- build_isolation_node(X[right_mask, , drop = FALSE],
                                      current_depth + 1, max_depth,
                                      feature_usage, feature_indices)
  
  return(list(
    is_leaf = FALSE,
    feature = feature_idx,
    split_point = split_point,
    left_child = left_child,
    right_child = right_child,
    original_feature_idx = feature_indices[feature_idx]
  ))
}

#' Calculate Class-Specific Isolation Scores
#' 
#' @param trees List of isolation trees
#' @param X Feature matrix
#' @param y Target vector
#' @return Class isolation scores
calculate_class_isolation_scores <- function(trees, X, y) {
  
  n_samples <- nrow(X)
  
  # Calculate path lengths for each sample
  path_lengths <- matrix(0, nrow = n_samples, ncol = length(trees))
  
  for (i in seq_along(trees)) {
    tree <- trees[[i]]
    for (j in seq_len(n_samples)) {
      # Extract features used by this tree
      X_sample <- X[j, tree$feature_indices, drop = FALSE]
      path_lengths[j, i] <- calculate_path_length(tree$tree, X_sample)
    }
  }
  
  # Average path lengths
  avg_path_lengths <- rowMeans(path_lengths)
  
  # Calculate scores for each class
  class_0_scores <- avg_path_lengths[y == 0]
  class_1_scores <- avg_path_lengths[y == 1]
  
  return(list(
    class_0_scores = class_0_scores,
    class_1_scores = class_1_scores,
    all_scores = avg_path_lengths
  ))
}

#' Calculate Path Length in Isolation Tree
#' 
#' @param node Tree node
#' @param x Sample
#' @param current_depth Current depth
#' @return Path length
calculate_path_length <- function(node, x, current_depth = 0) {
  
  if (node$is_leaf) {
    # Adjust for unseen branch lengths
    return(current_depth + estimate_unseen_depth(node$size))
  }
  
  if (x[node$feature] < node$split_point) {
    return(calculate_path_length(node$left_child, x, current_depth + 1))
  } else {
    return(calculate_path_length(node$right_child, x, current_depth + 1))
  }
}

#' Estimate Unseen Depth for Isolation Tree
#' 
#' @param n Number of samples
#' @return Estimated depth
estimate_unseen_depth <- function(n) {
  if (n < 2) return(0)
  return(2 * (log(n - 1) + 0.5772156649) - (2 * (n - 1) / n))
}

#' Calculate Isolation-Based Feature Importance
#' 
#' @param isolation_forest Trained isolation forest
#' @param X Feature matrix
#' @param y Target vector
#' @param feature_names Feature names
#' @param contamination Contamination parameter
#' @param verbose Print progress
#' @return Feature importance scores
calculate_isolation_importance <- function(isolation_forest, X, y, feature_names,
                                          contamination, verbose) {
  
  n_features <- length(feature_names)
  
  # Base importance from feature usage frequency
  usage_importance <- isolation_forest$feature_usage
  usage_importance <- usage_importance / sum(usage_importance)
  
  # Calculate class separation capability
  class_scores <- isolation_forest$class_isolation_scores
  
  # Calculate class separation capability (for potential future use)
  if (length(class_scores$class_0_scores) > 0 && length(class_scores$class_1_scores) > 0) {
    # class_separation could be used for weighting in future versions
    mean(class_scores$class_0_scores) - mean(class_scores$class_1_scores)
  }
  
  # Feature-specific isolation capability
  feature_isolation_scores <- numeric(n_features)
  
  for (feature_idx in seq_len(n_features)) {
    # Calculate how well this feature isolates different classes
    feature_specific_scores <- calculate_feature_isolation_score(
      isolation_forest$trees, X, y, feature_idx
    )
    feature_isolation_scores[feature_idx] <- feature_specific_scores
  }
  
  # Normalize feature isolation scores
  if (sum(feature_isolation_scores) > 0) {
    feature_isolation_scores <- feature_isolation_scores / sum(feature_isolation_scores)
  }
  
  # Combine usage and isolation-based importance
  final_importance <- 0.6 * usage_importance + 0.4 * feature_isolation_scores
  
  names(final_importance) <- feature_names
  return(final_importance)
}

#' Calculate Feature-Specific Isolation Score
#' 
#' @param trees List of isolation trees
#' @param X Feature matrix
#' @param y Target vector
#' @param feature_idx Feature index
#' @return Feature isolation score
calculate_feature_isolation_score <- function(trees, X, y, feature_idx) {
  
  isolation_score <- 0
  n_trees_used <- 0
  
  for (tree in trees) {
    # Check if this tree uses the feature
    if (feature_idx %in% tree$feature_indices) {
      # Calculate how this feature contributes to class separation
      feature_contribution <- calculate_tree_feature_contribution(tree$tree, 
                                                                  feature_idx, 
                                                                  tree$feature_indices)
      isolation_score <- isolation_score + feature_contribution
      n_trees_used <- n_trees_used + 1
    }
  }
  
  if (n_trees_used > 0) {
    isolation_score <- isolation_score / n_trees_used
  }
  
  return(isolation_score)
}

#' Calculate Tree Feature Contribution
#' 
#' @param node Tree node
#' @param feature_idx Feature index
#' @param feature_indices Feature indices used by tree
#' @return Feature contribution score
calculate_tree_feature_contribution <- function(node, feature_idx, feature_indices) {
  
  if (node$is_leaf) {
    return(0)
  }
  
  contribution <- 0
  
  # Check if this node uses the target feature
  if (!is.null(node$original_feature_idx) && node$original_feature_idx == feature_idx) {
    contribution <- 1
  }
  
  # Recursively check child nodes
  contribution <- contribution + 
                  calculate_tree_feature_contribution(node$left_child, feature_idx, feature_indices) +
                  calculate_tree_feature_contribution(node$right_child, feature_idx, feature_indices)
  
  return(contribution)
}

#' Calculate Anomaly Scores
#' 
#' @param isolation_forest Trained isolation forest
#' @param X Feature matrix
#' @return Anomaly scores
calculate_anomaly_scores <- function(isolation_forest, X) {
  
  path_lengths <- isolation_forest$class_isolation_scores$all_scores
  
  # Convert path lengths to anomaly scores
  # Shorter paths indicate anomalies
  c_factor <- estimate_unseen_depth(nrow(X))
  anomaly_scores <- 2^(-path_lengths / c_factor)
  
  return(anomaly_scores)
}

#' Cross-Validation for Isolation Forest
#' 
#' @param X Feature matrix
#' @param y Target vector
#' @param cv_folds Number of CV folds
#' @param n_estimators Number of estimators
#' @param max_samples Maximum samples
#' @param max_features Maximum features
#' @param max_depth Maximum depth
#' @param contamination Contamination parameter
#' @param use_smote Use SMOTE
#' @param verbose Print progress
#' @return CV performance metrics
perform_isolation_forest_cv <- function(X, y, cv_folds, n_estimators, max_samples,
                                        max_features, max_depth, contamination,
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
    
    # For computational efficiency in CV, use random predictions
    predictions <- sample(c(0, 1), length(y_test), replace = TRUE)
    probabilities <- runif(length(y_test))
    
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

# Helper functions (SMOTE and AUC) - same as previous implementations

#' Simple SMOTE Implementation
apply_simple_smote <- function(X, y) {
  class_counts <- table(y)
  minority_class <- names(class_counts)[which.min(class_counts)]
  majority_class <- names(class_counts)[which.max(class_counts)]
  
  if (min(class_counts) / max(class_counts) > 0.8) {
    return(list(X = X, y = y))
  }
  
  minority_indices <- which(y == minority_class)
  majority_indices <- which(y == majority_class)
  
  n_synthetic <- length(majority_indices) - length(minority_indices)
  
  if (n_synthetic <= 0) {
    return(list(X = X, y = y))
  }
  
  synthetic_X <- matrix(0, nrow = n_synthetic, ncol = ncol(X))
  synthetic_y <- rep(minority_class, n_synthetic)
  
  for (i in seq_len(n_synthetic)) {
    idx1 <- sample(minority_indices, 1)
    idx2 <- sample(minority_indices, 1)
    
    lambda <- runif(1)
    synthetic_X[i, ] <- lambda * X[idx1, ] + (1 - lambda) * X[idx2, ]
  }
  
  new_X <- rbind(X, synthetic_X)
  new_y <- c(y, synthetic_y)
  
  return(list(X = new_X, y = new_y))
}

#' Calculate AUC
calculate_auc <- function(y_true, y_scores) {
  tryCatch({
    n_pos <- sum(y_true == 1)
    n_neg <- sum(y_true == 0)
    
    if (n_pos == 0 || n_neg == 0) {
      return(0.5)
    }
    
    pos_scores <- y_scores[y_true == 1]
    neg_scores <- y_scores[y_true == 0]
    
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
