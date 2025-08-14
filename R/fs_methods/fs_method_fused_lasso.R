#' Fused LASSO Feature Selection Method
#'
#' Implements Fused LASSO feature selection which combines standard L1 penalty
#' with an additional penalty that encourages smoothness between adjacent features.
#' Assumes features can be ordered (e.g., genomic positions, time series).
#'
#' @param ramwas.data Expression data matrix (samples as rows, features as columns)  
#' @param target Target variable vector
#' @param smote_sampling Logical, whether to apply SMOTE sampling
#' @return List containing selected features and method metadata
#' @export
fused_lasso_fs <- function(ramwas.data, target, smote_sampling = FALSE) {
  
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
    
    # Perform Fused LASSO feature selection
    fused_lasso_result <- perform_fused_lasso(data_clean, target_numeric)
    
    # Extract selected features
    if (length(fused_lasso_result$selected_features) == 0) {
      return(create_empty_result("No features selected by Fused LASSO"))
    }
    
    # Calculate timing
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      selected_features = fused_lasso_result$selected_features,
      feature_scores = fused_lasso_result$feature_scores,
      method_name = "Fused_LASSO",
      execution_time = execution_time,
      n_features_selected = length(fused_lasso_result$selected_features),
      n_features_input = ncol(ramwas.data),
      method_params = list(
        lambda1 = fused_lasso_result$lambda1_optimal,
        lambda2 = fused_lasso_result$lambda2_optimal,
        cv_folds = fused_lasso_result$cv_folds,
        fusion_penalty = "adjacent_differences",
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
      method_name = "Fused_LASSO",
      execution_time = execution_time,
      n_features_selected = 0,
      n_features_input = ncol(ramwas.data),
      error_message = paste("Fused LASSO failed:", e$message),
      method_params = list(smote_sampling = smote_sampling)
    ))
  })
}

#' Perform Fused LASSO Feature Selection
#'
#' Implements Fused LASSO with dual penalties: L1 for sparsity and fusion for smoothness
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with selected features and optimal parameters
perform_fused_lasso <- function(data_matrix, target_vector) {
  
  n_samples <- nrow(data_matrix)
  n_features <- ncol(data_matrix)
  
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
    n_features <- ncol(data_scaled)
  }
  
  # Center target
  target_centered <- target_vector - mean(target_vector)
  
  # Create fusion matrix for adjacent differences
  fusion_matrix <- create_fusion_matrix(n_features)
  
  # Generate lambda sequences for both penalties
  lambda1_max <- calculate_lambda_max_fused(data_scaled, target_centered)
  lambda2_max <- lambda1_max * 0.5  # Fusion penalty typically smaller
  
  lambda1_sequence <- exp(seq(log(lambda1_max), log(lambda1_max * 0.01), length.out = 20))
  lambda2_sequence <- exp(seq(log(lambda2_max), log(lambda2_max * 0.01), length.out = 15))
  
  # Grid search with cross-validation
  cv_folds <- min(5, n_samples)
  fold_assignments <- sample(rep(1:cv_folds, length.out = n_samples))
  
  best_error <- Inf
  best_lambda1 <- lambda1_sequence[1]
  best_lambda2 <- lambda2_sequence[1]
  
  for (lambda1 in lambda1_sequence) {
    for (lambda2 in lambda2_sequence) {
      
      fold_errors <- numeric(cv_folds)
      
      for (fold in 1:cv_folds) {
        train_idx <- which(fold_assignments != fold)
        test_idx <- which(fold_assignments == fold)
        
        if (length(train_idx) < 2 || length(test_idx) < 1) next
        
        # Fit Fused LASSO on training set
        coefficients <- fused_lasso_solver(
          data_scaled[train_idx, , drop = FALSE], 
          target_centered[train_idx], 
          fusion_matrix,
          lambda1, 
          lambda2
        )
        
        # Predict on test set
        predictions <- data_scaled[test_idx, , drop = FALSE] %*% coefficients
        fold_errors[fold] <- mean((target_centered[test_idx] - predictions)^2)
      }
      
      avg_error <- mean(fold_errors[fold_errors > 0])
      if (avg_error < best_error) {
        best_error <- avg_error
        best_lambda1 <- lambda1
        best_lambda2 <- lambda2
      }
    }
  }
  
  # Fit final model with optimal lambdas
  final_coefficients <- fused_lasso_solver(
    data_scaled, target_centered, fusion_matrix, best_lambda1, best_lambda2
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
    lambda1_optimal = best_lambda1,
    lambda2_optimal = best_lambda2,
    cv_folds = cv_folds,
    coefficients = final_coefficients
  ))
}

#' Create Fusion Matrix for Adjacent Differences
#'
#' Creates difference matrix for Fused LASSO fusion penalty
#'
#' @param n_features Number of features
#' @return Fusion matrix for adjacent differences
create_fusion_matrix <- function(n_features) {
  
  if (n_features < 2) {
    return(matrix(0, nrow = 0, ncol = n_features))
  }
  
  # Create (n-1) x n difference matrix
  fusion_matrix <- matrix(0, nrow = n_features - 1, ncol = n_features)
  
  for (i in 1:(n_features - 1)) {
    fusion_matrix[i, i] <- -1
    fusion_matrix[i, i + 1] <- 1
  }
  
  return(fusion_matrix)
}

#' Fused LASSO Solver
#'
#' Implements coordinate descent for Fused LASSO with dual penalties
#'
#' @param X Design matrix (standardized)
#' @param y Response vector (centered)
#' @param D Fusion matrix for adjacent differences
#' @param lambda1 L1 penalty parameter
#' @param lambda2 Fusion penalty parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return Coefficient vector
fused_lasso_solver <- function(X, y, D, lambda1, lambda2, max_iter = 1000, tol = 1e-6) {
  
  n_features <- ncol(X)
  
  # Initialize coefficients
  beta <- numeric(n_features)
  beta_old <- beta
  
  # Precompute matrices
  XtX <- crossprod(X)
  Xty <- as.numeric(crossprod(X, y))
  DtD <- crossprod(D)
  
  # Coordinate descent iterations
  for (iter in 1:max_iter) {
    
    for (j in 1:n_features) {
      
      # Calculate residual for coordinate j
      residual_j <- Xty[j] - sum(XtX[j, -j] * beta[-j])
      
      # Add fusion penalty contribution
      fusion_contribution <- 0
      if (nrow(D) > 0) {
        fusion_contribution <- lambda2 * sum(DtD[j, -j] * beta[-j])
      }
      
      # Total residual including fusion penalty
      total_residual <- residual_j - fusion_contribution
      
      # Diagonal element with fusion penalty
      diagonal_element <- XtX[j, j] + lambda2 * DtD[j, j]
      
      # Soft thresholding with fused penalty
      if (diagonal_element > 0) {
        beta[j] <- soft_threshold_fused(total_residual, lambda1) / diagonal_element
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

#' Soft Thresholding Function for Fused LASSO
#'
#' Applies soft thresholding operator for L1 penalty
#'
#' @param x Input value
#' @param lambda Threshold parameter
#' @return Soft thresholded value
soft_threshold_fused <- function(x, lambda) {
  return(sign(x) * pmax(0, abs(x) - lambda))
}

#' Calculate Maximum Lambda for Fused LASSO
#'
#' Computes the maximum lambda1 value for Fused LASSO regularization path
#'
#' @param X Design matrix
#' @param y Response vector
#' @return Maximum lambda1 value
calculate_lambda_max_fused <- function(X, y) {
  return(max(abs(crossprod(X, y))) / nrow(X))
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

#' Create Empty Result for Fused LASSO
#'
#' Creates standardized empty result for failed feature selection
#'
#' @param message Error or warning message
#' @return Empty result list
create_empty_result <- function(message) {
  return(list(
    selected_features = character(0),
    feature_scores = numeric(0),
    method_name = "Fused_LASSO",
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
