#' Ridge Regression Feature Selection Method
#'
#' Implements Ridge regression feature selection using L2 regularization.
#' Features are ranked by the magnitude of their ridge coefficients, and
#' top features are selected based on coefficient stability across different
#' lambda values.
#'
#' @param ramwas.data Expression data matrix (samples as rows, features as columns)  
#' @param target Target variable vector
#' @param smote_sampling Logical, whether to apply SMOTE sampling
#' @return List containing selected features and method metadata
#' @export
ridge_regression_fs <- function(ramwas.data, target, smote_sampling = FALSE) {
  
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
    
    # Convert target to appropriate format
    target_numeric <- handle_target_variable(target)
    
    # Perform Ridge regression feature selection
    ridge_result <- perform_ridge_regression(data_clean, target_numeric)
    
    # Extract selected features
    if (length(ridge_result$selected_features) == 0) {
      return(create_empty_result("No features selected by Ridge regression"))
    }
    
    # Calculate timing
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      selected_features = ridge_result$selected_features,
      feature_scores = ridge_result$feature_scores,
      method_name = "Ridge_Regression",
      execution_time = execution_time,
      n_features_selected = length(ridge_result$selected_features),
      n_features_input = ncol(ramwas.data),
      method_params = list(
        lambda = ridge_result$lambda_optimal,
        cv_folds = ridge_result$cv_folds,
        alpha = 0.0,
        selection_method = "coefficient_magnitude",
        smote_sampling = smote_sampling
      ),
      preprocessing_info = list(
        constant_features_removed = preprocessed$constant_removed,
        correlated_features_removed = preprocessed$correlated_removed,
        final_features_count = ncol(data_clean)
      )
    ))
    
  }, error = function(e) {
    # Return error result
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      selected_features = character(0),
      feature_scores = numeric(0),
      method_name = "Ridge_Regression",
      execution_time = execution_time,
      n_features_selected = 0,
      n_features_input = ncol(ramwas.data),
      error_message = paste("Ridge regression failed:", e$message),
      method_params = list(smote_sampling = smote_sampling)
    ))
  })
}

#' Perform Ridge Regression Feature Selection
#'
#' Implements Ridge regression with cross-validation for lambda selection
#' and feature ranking based on coefficient magnitudes
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with selected features and optimal lambda
perform_ridge_regression <- function(data_matrix, target_vector) {
  
  n_samples <- nrow(data_matrix)
  
  # Standardize features
  data_scaled <- scale(data_matrix)
  feature_means <- attr(data_scaled, "scaled:center")
  feature_sds <- attr(data_scaled, "scaled:scale")
  
  # Handle any NAs from scaling
  if (any(is.na(data_scaled)) || any(feature_sds == 0)) {
    non_zero_sd <- which(feature_sds > 1e-10)
    data_scaled <- data_scaled[, non_zero_sd, drop = FALSE]
    feature_means <- feature_means[non_zero_sd]
    feature_sds <- feature_sds[non_zero_sd]
  }
  
  # Center target
  target_centered <- target_vector - mean(target_vector)
  
  # Generate lambda sequence for Ridge
  lambda_max <- calculate_lambda_max_ridge(data_scaled, target_centered)
  lambda_sequence <- exp(seq(log(lambda_max), log(lambda_max * 0.001), length.out = 50))
  
  # Perform cross-validation to find optimal lambda
  cv_folds <- min(5, n_samples)
  fold_assignments <- sample(rep(1:cv_folds, length.out = n_samples))
  
  cv_errors <- numeric(length(lambda_sequence))
  
  for (i in seq_along(lambda_sequence)) {
    lambda <- lambda_sequence[i]
    fold_errors <- numeric(cv_folds)
    
    for (fold in 1:cv_folds) {
      train_idx <- which(fold_assignments != fold)
      test_idx <- which(fold_assignments == fold)
      
      if (length(train_idx) < 2 || length(test_idx) < 1) next
      
      # Fit Ridge on training set
      coefficients <- ridge_regression_solver(
        data_scaled[train_idx, , drop = FALSE], 
        target_centered[train_idx], 
        lambda
      )
      
      # Predict on test set
      predictions <- data_scaled[test_idx, , drop = FALSE] %*% coefficients
      fold_errors[fold] <- mean((target_centered[test_idx] - predictions)^2)
    }
    
    cv_errors[i] <- mean(fold_errors[fold_errors > 0])
  }
  
  # Select optimal lambda (minimum CV error)
  optimal_idx <- which.min(cv_errors)
  lambda_optimal <- lambda_sequence[optimal_idx]
  
  # Fit final model with optimal lambda
  final_coefficients <- ridge_regression_solver(data_scaled, target_centered, lambda_optimal)
  
  # Calculate feature stability across multiple lambda values
  stability_scores <- calculate_ridge_stability(data_scaled, target_centered, lambda_sequence)
  
  # Combine coefficient magnitude and stability for feature ranking
  feature_importance <- abs(final_coefficients) * stability_scores
  
  # Select top features (default: top 20% or at least 10 features)
  n_features_to_select <- max(10, round(0.2 * length(feature_importance)))
  n_features_to_select <- min(n_features_to_select, length(feature_importance))
  
  top_indices <- order(feature_importance, decreasing = TRUE)[1:n_features_to_select]
  selected_features <- colnames(data_matrix)[top_indices]
  feature_scores <- feature_importance[top_indices]
  names(feature_scores) <- selected_features
  
  return(list(
    selected_features = selected_features,
    feature_scores = feature_scores,
    lambda_optimal = lambda_optimal,
    cv_folds = cv_folds,
    coefficients = final_coefficients,
    stability_scores = stability_scores
  ))
}

#' Ridge Regression Solver
#'
#' Solves Ridge regression using analytical solution
#'
#' @param X Design matrix (standardized)
#' @param y Response vector (centered)
#' @param lambda Regularization parameter
#' @return Coefficient vector
ridge_regression_solver <- function(X, y, lambda) {
  
  n_features <- ncol(X)
  
  # Ridge regression analytical solution: beta = (X'X + lambda*I)^(-1) X'y
  XtX <- crossprod(X)
  Xty <- as.numeric(crossprod(X, y))
  
  # Add ridge penalty to diagonal
  ridge_matrix <- XtX + lambda * diag(n_features)
  
  # Solve for coefficients
  tryCatch({
    coefficients <- solve(ridge_matrix, Xty)
    return(coefficients)
  }, error = function(e) {
    # Fallback to pseudo-inverse if matrix is singular
    ridge_matrix_pseudo <- ridge_matrix + 1e-6 * diag(n_features)
    coefficients <- solve(ridge_matrix_pseudo, Xty)
    return(coefficients)
  })
}

#' Calculate Ridge Feature Stability
#'
#' Calculates feature stability across different lambda values
#'
#' @param X Design matrix
#' @param y Response vector
#' @param lambda_sequence Vector of lambda values
#' @return Vector of stability scores
calculate_ridge_stability <- function(X, y, lambda_sequence) {
  
  n_features <- ncol(X)
  coefficient_matrix <- matrix(0, nrow = length(lambda_sequence), ncol = n_features)
  
  # Fit Ridge regression for each lambda
  for (i in seq_along(lambda_sequence)) {
    coefficients <- ridge_regression_solver(X, y, lambda_sequence[i])
    coefficient_matrix[i, ] <- coefficients
  }
  
  # Calculate stability as inverse of coefficient variance
  coefficient_vars <- apply(coefficient_matrix, 2, var)
  stability_scores <- 1 / (1 + coefficient_vars)
  
  return(stability_scores)
}

#' Calculate Maximum Lambda for Ridge
#'
#' Computes a reasonable maximum lambda value for Ridge regression
#'
#' @param X Design matrix
#' @param y Response vector
#' @return Maximum lambda value
calculate_lambda_max_ridge <- function(X, y) {
  # For Ridge, lambda_max is typically based on the largest eigenvalue
  XtX <- crossprod(X)
  max_eigenvalue <- max(eigen(XtX, symmetric = TRUE, only.values = TRUE)$values)
  return(max_eigenvalue * 0.1)  # Conservative starting point
}

#' Basic Feature Preprocessing for Ridge
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

#' Handle Target Variable Conversion
#'
#' Converts target variable to appropriate numeric format
#'
#' @param target Target variable
#' @return Numeric target variable
handle_target_variable <- function(target) {
  if (is.factor(target)) {
    return(as.numeric(target) - 1)
  } else if (is.character(target)) {
    return(as.numeric(as.factor(target)) - 1)
  } else {
    return(as.numeric(target))
  }
}

#' Create Empty Result for Ridge
#'
#' Creates standardized empty result for failed feature selection
#'
#' @param message Error or warning message
#' @return Empty result list
create_empty_result <- function(message) {
  return(list(
    selected_features = character(0),
    feature_scores = numeric(0),
    method_name = "Ridge_Regression",
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
