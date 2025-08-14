#' Adaptive LASSO Feature Selection Method
#'
#' Implements Adaptive LASSO feature selection which uses weighted L1 penalty
#' where weights are derived from initial estimator (e.g., ridge regression).
#' This method has oracle properties and can achieve variable selection consistency.
#'
#' @param ramwas.data Expression data matrix (samples as rows, features as columns)  
#' @param target Target variable vector
#' @param smote_sampling Logical, whether to apply SMOTE sampling
#' @return List containing selected features and method metadata
#' @export
adaptive_lasso_fs <- function(ramwas.data, target, smote_sampling = FALSE) {
  
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
    
    # Perform Adaptive LASSO feature selection
    adaptive_lasso_result <- perform_adaptive_lasso(data_clean, target_numeric)
    
    # Extract selected features
    if (length(adaptive_lasso_result$selected_features) == 0) {
      return(create_empty_result("No features selected by Adaptive LASSO"))
    }
    
    # Calculate timing
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      selected_features = adaptive_lasso_result$selected_features,
      feature_scores = adaptive_lasso_result$feature_scores,
      method_name = "Adaptive_LASSO",
      execution_time = execution_time,
      n_features_selected = length(adaptive_lasso_result$selected_features),
      n_features_input = ncol(ramwas.data),
      method_params = list(
        lambda = adaptive_lasso_result$lambda_optimal,
        cv_folds = adaptive_lasso_result$cv_folds,
        weight_method = "ridge_based",
        gamma_power = adaptive_lasso_result$gamma_power,
        smote_sampling = smote_sampling
      ),
      preprocessing_info = list(
        constant_features_removed = preprocessed$constant_removed,
        correlated_features_removed = preprocessed$correlated_removed,
        final_features_count = ncol(data_clean)
      ),
      adaptive_weights = adaptive_lasso_result$adaptive_weights
    ))
    
  }, error = function(e) {
    # Return error result
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      selected_features = character(0),
      feature_scores = numeric(0),
      method_name = "Adaptive_LASSO",
      execution_time = execution_time,
      n_features_selected = 0,
      n_features_input = ncol(ramwas.data),
      error_message = paste("Adaptive LASSO failed:", e$message),
      method_params = list(smote_sampling = smote_sampling)
    ))
  })
}

#' Perform Adaptive LASSO Feature Selection
#'
#' Implements Adaptive LASSO with ridge-based weights and cross-validation
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with selected features and optimal parameters
perform_adaptive_lasso <- function(data_matrix, target_vector) {
  
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
  
  # Step 1: Compute initial weights using ridge regression
  ridge_weights <- compute_ridge_weights(data_scaled, target_centered)
  
  # Step 2: Compute adaptive weights
  gamma_power <- 1.0  # Standard choice for gamma
  adaptive_weights <- compute_adaptive_weights(ridge_weights, gamma_power)
  
  # Step 3: Generate lambda sequence for adaptive LASSO
  lambda_max <- calculate_lambda_max_adaptive(data_scaled, target_centered, adaptive_weights)
  lambda_sequence <- exp(seq(log(lambda_max), log(lambda_max * 0.001), length.out = 50))
  
  # Step 4: Cross-validation to find optimal lambda
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
      
      # Fit Adaptive LASSO on training set
      coefficients <- adaptive_lasso_solver(
        data_scaled[train_idx, , drop = FALSE], 
        target_centered[train_idx], 
        adaptive_weights,
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
  
  # Step 5: Fit final model with optimal lambda
  final_coefficients <- adaptive_lasso_solver(
    data_scaled, target_centered, adaptive_weights, lambda_optimal
  )
  
  # Identify selected features (non-zero coefficients)
  selected_idx <- which(abs(final_coefficients) > 1e-6)
  
  if (length(selected_idx) > 0) {
    selected_features <- colnames(data_matrix)[selected_idx]
    feature_scores <- abs(final_coefficients[selected_idx])
    names(feature_scores) <- selected_features
  } else {
    selected_features <- character(0)
    feature_scores <- numeric(0)
  }
  
  return(list(
    selected_features = selected_features,
    feature_scores = feature_scores,
    lambda_optimal = lambda_optimal,
    cv_folds = cv_folds,
    coefficients = final_coefficients,
    adaptive_weights = adaptive_weights,
    gamma_power = gamma_power
  ))
}

#' Compute Ridge Weights for Adaptive LASSO
#'
#' Computes initial weights using ridge regression estimates
#'
#' @param X Design matrix (standardized)
#' @param y Response vector (centered)
#' @return Vector of ridge coefficient magnitudes
compute_ridge_weights <- function(X, y) {
  
  # Use cross-validation to select ridge lambda
  lambda_ridge_max <- max(eigen(crossprod(X), symmetric = TRUE, only.values = TRUE)$values) * 0.1
  lambda_ridge_sequence <- exp(seq(log(lambda_ridge_max), log(lambda_ridge_max * 0.001), length.out = 20))
  
  # Simple cross-validation for ridge
  cv_folds <- min(5, nrow(X))
  fold_assignments <- sample(rep(1:cv_folds, length.out = nrow(X)))
  
  best_error <- Inf
  best_lambda_ridge <- lambda_ridge_sequence[1]
  
  for (lambda_ridge in lambda_ridge_sequence) {
    fold_errors <- numeric(cv_folds)
    
    for (fold in 1:cv_folds) {
      train_idx <- which(fold_assignments != fold)
      test_idx <- which(fold_assignments == fold)
      
      if (length(train_idx) < 2 || length(test_idx) < 1) next
      
      # Fit ridge regression
      X_train <- X[train_idx, , drop = FALSE]
      y_train <- y[train_idx]
      
      ridge_coeffs <- ridge_regression_solver(X_train, y_train, lambda_ridge)
      
      # Predict and calculate error
      predictions <- X[test_idx, , drop = FALSE] %*% ridge_coeffs
      fold_errors[fold] <- mean((y[test_idx] - predictions)^2)
    }
    
    avg_error <- mean(fold_errors[fold_errors > 0])
    if (avg_error < best_error) {
      best_error <- avg_error
      best_lambda_ridge <- lambda_ridge
    }
  }
  
  # Fit final ridge model
  final_ridge_coeffs <- ridge_regression_solver(X, y, best_lambda_ridge)
  
  return(abs(final_ridge_coeffs))
}

#' Ridge Regression Solver
#'
#' Solves ridge regression using analytical solution
#'
#' @param X Design matrix
#' @param y Response vector
#' @param lambda Ridge penalty parameter
#' @return Coefficient vector
ridge_regression_solver <- function(X, y, lambda) {
  
  n_features <- ncol(X)
  
  # Ridge regression: beta = (X'X + lambda*I)^(-1) X'y
  XtX <- crossprod(X)
  Xty <- as.numeric(crossprod(X, y))
  
  ridge_matrix <- XtX + lambda * diag(n_features)
  
  tryCatch({
    coefficients <- solve(ridge_matrix, Xty)
    return(coefficients)
  }, error = function(e) {
    # Fallback with regularization
    ridge_matrix_reg <- ridge_matrix + 1e-6 * diag(n_features)
    coefficients <- solve(ridge_matrix_reg, Xty)
    return(coefficients)
  })
}

#' Compute Adaptive Weights
#'
#' Computes adaptive weights from initial estimates
#'
#' @param initial_coeffs Initial coefficient estimates
#' @param gamma_power Power parameter for weights
#' @return Vector of adaptive weights
compute_adaptive_weights <- function(initial_coeffs, gamma_power = 1.0) {
  
  # Adaptive weights: w_j = 1 / |beta_j|^gamma
  # Add small constant to avoid division by zero
  weights <- 1 / (abs(initial_coeffs) + 1e-6)^gamma_power
  
  return(weights)
}

#' Adaptive LASSO Solver
#'
#' Implements coordinate descent for Adaptive LASSO with weighted penalties
#'
#' @param X Design matrix (standardized)
#' @param y Response vector (centered)
#' @param weights Adaptive weight vector
#' @param lambda Regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return Coefficient vector
adaptive_lasso_solver <- function(X, y, weights, lambda, max_iter = 1000, tol = 1e-6) {
  
  n_features <- ncol(X)
  
  # Initialize coefficients
  beta <- numeric(n_features)
  beta_old <- beta
  
  # Precompute X^T X diagonal and X^T y
  XtX_diag <- colSums(X^2)
  Xty <- as.numeric(crossprod(X, y))
  
  # Coordinate descent iterations
  for (iter in 1:max_iter) {
    
    for (j in 1:n_features) {
      
      # Calculate residual excluding feature j
      residual_j <- Xty[j] - sum(crossprod(X[, j], X[, -j, drop = FALSE]) * beta[-j])
      
      # Adaptive soft thresholding
      if (XtX_diag[j] > 0) {
        weighted_lambda <- lambda * weights[j]
        beta[j] <- soft_threshold_adaptive(residual_j, weighted_lambda) / XtX_diag[j]
      } else {
        beta[j] <- 0
      }
    }
    
    # Check convergence
    if (max(abs(beta - beta_old)) < tol) {
      break
    }
    
    beta_old <- beta
  }
  
  return(beta)
}

#' Soft Thresholding Function for Adaptive LASSO
#'
#' Applies soft thresholding operator with adaptive weights
#'
#' @param x Input value
#' @param lambda Weighted threshold parameter
#' @return Soft thresholded value
soft_threshold_adaptive <- function(x, lambda) {
  return(sign(x) * pmax(0, abs(x) - lambda))
}

#' Calculate Maximum Lambda for Adaptive LASSO
#'
#' Computes the maximum lambda value for Adaptive LASSO regularization path
#'
#' @param X Design matrix
#' @param y Response vector
#' @param weights Adaptive weight vector
#' @return Maximum lambda value
calculate_lambda_max_adaptive <- function(X, y, weights) {
  
  # For adaptive LASSO: lambda_max = max(|X'y| / weights) / n
  Xty <- abs(crossprod(X, y))
  weighted_gradients <- Xty / (weights + 1e-6)  # Avoid division by zero
  
  return(max(weighted_gradients) / nrow(X))
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

#' Create Empty Result for Adaptive LASSO
#'
#' Creates standardized empty result for failed feature selection
#'
#' @param message Error or warning message
#' @return Empty result list
create_empty_result <- function(message) {
  return(list(
    selected_features = character(0),
    feature_scores = numeric(0),
    method_name = "Adaptive_LASSO",
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
