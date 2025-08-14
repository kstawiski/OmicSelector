#' Nuclear Norm Minimization Feature Selection Method
#'
#' Implements Nuclear Norm (trace norm) minimization for feature selection.
#' This method encourages low-rank structure in the coefficient matrix and
#' is particularly useful for multi-task learning or when features have
#' intrinsic low-rank relationships.
#'
#' @param ramwas.data Expression data matrix (samples as rows, features as columns)  
#' @param target Target variable vector
#' @param smote_sampling Logical, whether to apply SMOTE sampling
#' @return List containing selected features and method metadata
#' @export
nuclear_norm_fs <- function(ramwas.data, target, smote_sampling = FALSE) {
  
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
    
    # Perform Nuclear Norm feature selection
    nuclear_result <- perform_nuclear_norm_selection(data_clean, target_numeric)
    
    # Extract selected features
    if (length(nuclear_result$selected_features) == 0) {
      return(create_empty_result("No features selected by Nuclear Norm"))
    }
    
    # Calculate timing
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      selected_features = nuclear_result$selected_features,
      feature_scores = nuclear_result$feature_scores,
      method_name = "Nuclear_Norm",
      execution_time = execution_time,
      n_features_selected = length(nuclear_result$selected_features),
      n_features_input = ncol(ramwas.data),
      method_params = list(
        lambda = nuclear_result$lambda_optimal,
        rank_estimate = nuclear_result$rank_estimate,
        cv_folds = nuclear_result$cv_folds,
        regularization_type = "nuclear_norm",
        smote_sampling = smote_sampling
      ),
      preprocessing_info = list(
        constant_features_removed = preprocessed$constant_removed,
        correlated_features_removed = preprocessed$correlated_removed,
        final_features_count = ncol(data_clean)
      ),
      svd_info = nuclear_result$svd_info
    ))
    
  }, error = function(e) {
    # Return error result
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      selected_features = character(0),
      feature_scores = numeric(0),
      method_name = "Nuclear_Norm",
      execution_time = execution_time,
      n_features_selected = 0,
      n_features_input = ncol(ramwas.data),
      error_message = paste("Nuclear Norm failed:", e$message),
      method_params = list(smote_sampling = smote_sampling)
    ))
  })
}

#' Perform Nuclear Norm Feature Selection
#'
#' Implements nuclear norm minimization with iterative soft thresholding
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with selected features and nuclear norm information
perform_nuclear_norm_selection <- function(data_matrix, target_vector) {
  
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
  
  # Create augmented matrix for nuclear norm formulation
  # We form a matrix where features are related through low-rank structure
  augmented_matrix <- create_augmented_matrix(data_scaled, target_centered)
  
  # Generate lambda sequence for nuclear norm regularization
  lambda_max <- calculate_lambda_max_nuclear(augmented_matrix)
  lambda_sequence <- exp(seq(log(lambda_max), log(lambda_max * 0.01), length.out = 30))
  
  # Cross-validation to find optimal lambda
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
      
      # Fit nuclear norm model on training set
      nuclear_result <- nuclear_norm_solver(
        data_scaled[train_idx, , drop = FALSE], 
        target_centered[train_idx], 
        lambda
      )
      
      # Predict on test set
      if (length(nuclear_result$selected_features) > 0) {
        X_test <- data_scaled[test_idx, nuclear_result$selected_features, drop = FALSE]
        
        tryCatch({
          coefficients <- solve(crossprod(X_test) + 1e-6 * diag(ncol(X_test)),
                               crossprod(X_test, target_centered[test_idx]))
          predictions <- X_test %*% coefficients
          fold_errors[fold] <- mean((target_centered[test_idx] - predictions)^2)
        }, error = function(e) {
          fold_errors[fold] <- mean(target_centered[test_idx]^2)  # Baseline error
        })
      } else {
        fold_errors[fold] <- mean(target_centered[test_idx]^2)
      }
    }
    
    cv_errors[i] <- mean(fold_errors[fold_errors > 0])
  }
  
  # Select optimal lambda
  optimal_idx <- which.min(cv_errors)
  lambda_optimal <- lambda_sequence[optimal_idx]
  
  # Fit final model with optimal lambda
  final_result <- nuclear_norm_solver(data_scaled, target_centered, lambda_optimal)
  
  return(list(
    selected_features = final_result$selected_features,
    feature_scores = final_result$feature_scores,
    lambda_optimal = lambda_optimal,
    cv_folds = cv_folds,
    rank_estimate = final_result$rank_estimate,
    svd_info = final_result$svd_info
  ))
}

#' Create Augmented Matrix for Nuclear Norm
#'
#' Creates augmented matrix for nuclear norm formulation
#'
#' @param X Design matrix
#' @param y Response vector
#' @return Augmented matrix
create_augmented_matrix <- function(X, y) {
  
  # Create correlation-based similarity matrix for features
  cor_matrix <- cor(X, use = "complete.obs")
  cor_matrix[is.na(cor_matrix)] <- 0
  
  # Create augmented matrix with target relationship
  target_correlations <- cor(X, y, use = "complete.obs")
  target_correlations[is.na(target_correlations)] <- 0
  
  # Combine feature correlations and target correlations
  augmented_matrix <- cbind(cor_matrix, target_correlations)
  
  return(augmented_matrix)
}

#' Nuclear Norm Solver
#'
#' Solves nuclear norm minimization using iterative soft thresholding
#'
#' @param X Design matrix
#' @param y Response vector
#' @param lambda Nuclear norm penalty parameter
#' @return List with selected features and SVD information
nuclear_norm_solver <- function(X, y, lambda) {
  
  n_features <- ncol(X)
  
  # Create coefficient matrix (reshape problem for matrix completion)
  # We'll use the correlation structure to define a matrix completion problem
  
  # Start with ordinary least squares solution
  tryCatch({
    beta_ols <- solve(crossprod(X) + 1e-6 * diag(n_features), crossprod(X, y))
  }, error = function(e) {
    beta_ols <<- rep(0, n_features)
  })
  
  # Create matrix to apply nuclear norm to
  # We reshape the coefficient vector into a matrix form
  matrix_dim <- floor(sqrt(n_features))
  if (matrix_dim < 2) matrix_dim <- 2
  
  # Pad or truncate to make square matrix
  if (n_features >= matrix_dim^2) {
    coeff_matrix <- matrix(beta_ols[1:(matrix_dim^2)], nrow = matrix_dim, ncol = matrix_dim)
  } else {
    padded_coeffs <- c(beta_ols, rep(0, matrix_dim^2 - n_features))
    coeff_matrix <- matrix(padded_coeffs, nrow = matrix_dim, ncol = matrix_dim)
  }
  
  # Apply nuclear norm regularization using iterative soft thresholding
  regularized_matrix <- nuclear_norm_soft_threshold(coeff_matrix, lambda)
  
  # Extract coefficients back to vector form
  regularized_coeffs <- as.vector(regularized_matrix)[1:n_features]
  
  # Perform SVD for rank estimation and feature selection
  svd_result <- svd(regularized_matrix)
  rank_estimate <- sum(svd_result$d > 1e-6)
  
  # Select features based on coefficient magnitudes and SVD structure
  feature_importance <- abs(regularized_coeffs)
  
  # Also consider SVD-based importance
  if (rank_estimate > 0) {
    # Weight by leading singular vectors
    max_rank <- min(rank_estimate, ncol(svd_result$u))
    if (max_rank > 0) {
      left_singular <- abs(svd_result$u[, seq_len(max_rank), drop = FALSE])
      
      # Combine singular vector information
      if (nrow(left_singular) >= n_features) {
        svd_importance <- rowMeans(left_singular[seq_len(n_features), , drop = FALSE])
      } else {
        svd_importance <- rep(mean(left_singular), n_features)
      }
      
      # Combine coefficient and SVD importance
      feature_importance <- feature_importance + 0.5 * svd_importance
    }
  }
  
  # Select top features (default: features with non-zero importance)
  selected_idx <- which(feature_importance > 1e-6)
  
  if (length(selected_idx) == 0) {
    # Fallback: select top 10% of features
    n_select <- max(1, round(0.1 * n_features))
    selected_idx <- order(feature_importance, decreasing = TRUE)[1:n_select]
  }
  
  selected_features <- colnames(X)[selected_idx]
  feature_scores <- feature_importance[selected_idx]
  names(feature_scores) <- selected_features
  
  return(list(
    selected_features = selected_features,
    feature_scores = feature_scores,
    rank_estimate = rank_estimate,
    svd_info = list(
      singular_values = svd_result$d,
      rank = rank_estimate
    )
  ))
}

#' Nuclear Norm Soft Thresholding
#'
#' Applies soft thresholding to singular values for nuclear norm regularization
#'
#' @param matrix Input matrix
#' @param lambda Threshold parameter
#' @return Soft thresholded matrix
nuclear_norm_soft_threshold <- function(matrix, lambda) {
  
  # Perform SVD
  svd_result <- svd(matrix)
  
  # Apply soft thresholding to singular values
  thresholded_values <- pmax(0, svd_result$d - lambda)
  
  # Reconstruct matrix
  reconstructed <- svd_result$u %*% diag(thresholded_values) %*% t(svd_result$v)
  
  return(reconstructed)
}

#' Calculate Maximum Lambda for Nuclear Norm
#'
#' Computes the maximum lambda value for nuclear norm regularization
#'
#' @param matrix Input matrix for regularization
#' @return Maximum lambda value
calculate_lambda_max_nuclear <- function(matrix) {
  
  # Maximum lambda is the largest singular value
  svd_result <- svd(matrix)
  return(max(svd_result$d))
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

#' Create Empty Result for Nuclear Norm
#'
#' Creates standardized empty result for failed feature selection
#'
#' @param message Error or warning message
#' @return Empty result list
create_empty_result <- function(message) {
  return(list(
    selected_features = character(0),
    feature_scores = numeric(0),
    method_name = "Nuclear_Norm",
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
