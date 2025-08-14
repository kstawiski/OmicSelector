#' Least Angle Regression Enhanced Feature Selection Method
#'
#' Implements enhanced Least Angle Regression (LARS) feature selection with
#' improved stopping criteria and cross-validation for optimal feature selection.
#' This is an enhanced version with better regularization path tracking.
#'
#' @param ramwas.data Expression data matrix (samples as rows, features as columns)  
#' @param target Target variable vector
#' @param smote_sampling Logical, whether to apply SMOTE sampling
#' @return List containing selected features and method metadata
#' @export
lars_enhanced_fs <- function(ramwas.data, target, smote_sampling = FALSE) {
  
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
    
    # Perform LARS Enhanced feature selection
    lars_result <- perform_lars_enhanced(data_clean, target_numeric)
    
    # Extract selected features
    if (length(lars_result$selected_features) == 0) {
      return(create_empty_result("No features selected by LARS Enhanced"))
    }
    
    # Calculate timing
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      selected_features = lars_result$selected_features,
      feature_scores = lars_result$feature_scores,
      method_name = "LARS_Enhanced",
      execution_time = execution_time,
      n_features_selected = length(lars_result$selected_features),
      n_features_input = ncol(ramwas.data),
      method_params = list(
        n_steps = lars_result$n_steps,
        optimal_step = lars_result$optimal_step,
        cv_folds = lars_result$cv_folds,
        stopping_criterion = "cross_validation",
        smote_sampling = smote_sampling
      ),
      preprocessing_info = list(
        constant_features_removed = preprocessed$constant_removed,
        correlated_features_removed = preprocessed$correlated_removed,
        final_features_count = ncol(data_clean)
      ),
      lars_path = lars_result$lars_path
    ))
    
  }, error = function(e) {
    # Return error result
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      selected_features = character(0),
      feature_scores = numeric(0),
      method_name = "LARS_Enhanced",
      execution_time = execution_time,
      n_features_selected = 0,
      n_features_input = ncol(ramwas.data),
      error_message = paste("LARS Enhanced failed:", e$message),
      method_params = list(smote_sampling = smote_sampling)
    ))
  })
}

#' Perform Enhanced LARS Feature Selection
#'
#' Implements LARS algorithm with cross-validation for optimal step selection
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with selected features and LARS path information
perform_lars_enhanced <- function(data_matrix, target_vector) {
  
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
  
  # Compute LARS regularization path
  lars_path <- compute_lars_path_enhanced(data_scaled, target_centered)
  
  # Use cross-validation to select optimal step
  cv_folds <- min(5, n_samples)
  optimal_step <- select_optimal_lars_step(data_scaled, target_centered, lars_path, cv_folds)
  
  # Extract features and coefficients at optimal step
  if (optimal_step > 0 && optimal_step <= length(lars_path$steps)) {
    step_info <- lars_path$steps[[optimal_step]]
    active_set <- step_info$active_set
    coefficients <- step_info$coefficients
    
    if (length(active_set) > 0) {
      selected_features <- colnames(data_matrix)[active_set]
      feature_scores <- abs(coefficients[active_set])
      names(feature_scores) <- selected_features
    } else {
      selected_features <- character(0)
      feature_scores <- numeric(0)
    }
  } else {
    selected_features <- character(0)
    feature_scores <- numeric(0)
  }
  
  return(list(
    selected_features = selected_features,
    feature_scores = feature_scores,
    n_steps = length(lars_path$steps),
    optimal_step = optimal_step,
    cv_folds = cv_folds,
    lars_path = lars_path
  ))
}

#' Compute Enhanced LARS Regularization Path
#'
#' Computes the full LARS path with improved numerical stability
#'
#' @param X Design matrix (standardized)
#' @param y Response vector (centered)
#' @return List containing LARS path information
compute_lars_path_enhanced <- function(X, y) {
  
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  # Initialize LARS variables
  active_set <- c()
  inactive_set <- 1:n_features
  residual <- y
  beta <- numeric(n_features)
  
  # Store path information
  lars_steps <- list()
  step <- 0
  
  # Maximum steps (prevent infinite loops)
  max_steps <- min(n_features, n_samples - 1, 100)
  
  # LARS main loop
  while (length(active_set) < max_steps && length(inactive_set) > 0) {
    
    step <- step + 1
    
    # Calculate correlations with residual
    correlations <- numeric(n_features)
    for (j in inactive_set) {
      correlations[j] <- abs(sum(X[, j] * residual))
    }
    
    # Find maximum correlation
    max_correlation <- max(correlations[inactive_set])
    
    # Check for numerical issues
    if (max_correlation < 1e-10) {
      break
    }
    
    # Find features with maximum correlation (equiangular condition)
    equiangular_features <- inactive_set[abs(correlations[inactive_set] - max_correlation) < 1e-8]
    
    # Add feature(s) to active set
    new_feature <- equiangular_features[1]  # Take first if multiple
    active_set <- c(active_set, new_feature)
    inactive_set <- setdiff(inactive_set, new_feature)
    
    # Compute equiangular direction
    if (length(active_set) == 1) {
      # First step: simple direction
      direction <- sign(sum(X[, active_set] * residual))
      X_active <- matrix(X[, active_set], ncol = 1)
      equiangular_dir <- direction * X_active / sqrt(sum(X_active^2))
    } else {
      # Multiple features: compute equiangular direction
      X_active <- X[, active_set, drop = FALSE]
      
      # Calculate correlation signs
      correlation_signs <- sign(colSums(X_active * residual))
      
      # Solve for equiangular direction
      tryCatch({
        G <- crossprod(X_active)
        ones <- correlation_signs
        G_inv_ones <- solve(G, ones)
        normalization <- as.numeric(sqrt(crossprod(ones, G_inv_ones)))
        
        if (normalization > 1e-10) {
          equiangular_dir <- X_active %*% (G_inv_ones / normalization)
        } else {
          # Fallback: use principal component
          svd_result <- svd(X_active)
          equiangular_dir <- svd_result$u[, 1, drop = FALSE]
        }
      }, error = function(e) {
        # Numerical issues: use simpler approach
        equiangular_dir <- X_active[, ncol(X_active), drop = FALSE]
        equiangular_dir <- equiangular_dir / sqrt(sum(equiangular_dir^2))
      })
    }
    
    # Calculate step size
    if (length(inactive_set) > 0) {
      step_size <- calculate_lars_step_size(X, residual, equiangular_dir, active_set, inactive_set)
    } else {
      # Final step: go to the end
      step_size <- max_correlation / sum(equiangular_dir^2)
    }
    
    # Update coefficients and residual
    if (step_size > 0 && is.finite(step_size)) {
      # Update residual
      residual <- residual - step_size * equiangular_dir
      
      # Update coefficients for active features
      if (length(active_set) == 1) {
        beta[active_set] <- beta[active_set] + step_size * sign(sum(X[, active_set] * y))
      } else {
        # More complex update for multiple features
        X_active <- X[, active_set, drop = FALSE]
        correlation_signs <- sign(colSums(X_active * y))
        
        tryCatch({
          G <- crossprod(X_active)
          ones <- correlation_signs
          G_inv_ones <- solve(G, ones)
          normalization <- as.numeric(sqrt(crossprod(ones, G_inv_ones)))
          
          if (normalization > 1e-10) {
            coefficient_update <- step_size * G_inv_ones / normalization
            beta[active_set] <- beta[active_set] + coefficient_update
          }
        }, error = function(e) {
          # Simple fallback update
          beta[active_set[length(active_set)]] <- beta[active_set[length(active_set)]] + step_size
        })
      }
    }
    
    # Store step information
    lars_steps[[step]] <- list(
      active_set = active_set,
      coefficients = beta,
      residual_norm = sqrt(sum(residual^2)),
      max_correlation = max_correlation
    )
    
    # Check stopping criteria
    residual_norm <- sqrt(sum(residual^2))
    if (residual_norm < 1e-10) {
      break
    }
  }
  
  return(list(
    steps = lars_steps,
    final_active_set = active_set,
    final_coefficients = beta
  ))
}

#' Calculate LARS Step Size
#'
#' Computes the step size for LARS algorithm
#'
#' @param X Design matrix
#' @param residual Current residual
#' @param direction Equiangular direction
#' @param active_set Current active set
#' @param inactive_set Current inactive set
#' @return Step size
calculate_lars_step_size <- function(X, residual, direction, active_set, inactive_set) {
  
  current_correlation <- max(abs(colSums(X[, active_set, drop = FALSE] * residual)))
  
  min_step <- Inf
  
  # Check correlations with inactive features
  for (j in inactive_set) {
    X_j <- X[, j]
    
    # Calculate correlation change rate
    correlation_rate <- sum(X_j * direction)
    
    if (abs(correlation_rate) > 1e-10) {
      current_corr_j <- sum(X_j * residual)
      
      # Two cases: correlation increasing or decreasing toward max
      step1 <- (current_correlation - current_corr_j) / correlation_rate
      step2 <- (current_correlation + current_corr_j) / (-correlation_rate)
      
      for (step_candidate in c(step1, step2)) {
        if (step_candidate > 1e-10 && step_candidate < min_step) {
          min_step <- step_candidate
        }
      }
    }
  }
  
  return(if (is.finite(min_step)) min_step else current_correlation)
}

#' Select Optimal LARS Step using Cross-Validation
#'
#' Uses cross-validation to select the optimal step in LARS path
#'
#' @param X Design matrix
#' @param y Response vector
#' @param lars_path LARS path information
#' @param cv_folds Number of CV folds
#' @return Optimal step number
select_optimal_lars_step <- function(X, y, lars_path, cv_folds) {
  
  n_samples <- nrow(X)
  n_steps <- length(lars_path$steps)
  
  if (n_steps == 0) return(0)
  
  # Create fold assignments
  fold_assignments <- sample(rep(1:cv_folds, length.out = n_samples))
  
  # Calculate CV error for each step
  cv_errors <- numeric(n_steps)
  
  for (step in 1:n_steps) {
    step_info <- lars_path$steps[[step]]
    active_set <- step_info$active_set
    
    if (length(active_set) == 0) {
      cv_errors[step] <- Inf
      next
    }
    
    fold_errors <- numeric(cv_folds)
    
    for (fold in 1:cv_folds) {
      train_idx <- which(fold_assignments != fold)
      test_idx <- which(fold_assignments == fold)
      
      if (length(train_idx) < 2 || length(test_idx) < 1) next
      
      # Fit model on training set
      X_train <- X[train_idx, active_set, drop = FALSE]
      y_train <- y[train_idx]
      
      tryCatch({
        coefficients <- solve(crossprod(X_train) + 1e-6 * diag(length(active_set)),
                             crossprod(X_train, y_train))
        
        # Predict on test set
        X_test <- X[test_idx, active_set, drop = FALSE]
        predictions <- X_test %*% coefficients
        fold_errors[fold] <- mean((y[test_idx] - predictions)^2)
      }, error = function(e) {
        fold_errors[fold] <<- Inf
      })
    }
    
    cv_errors[step] <- mean(fold_errors[fold_errors < Inf])
  }
  
  # Select step with minimum CV error
  optimal_step <- which.min(cv_errors)
  
  return(if (length(optimal_step) > 0) optimal_step else 1)
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

#' Create Empty Result for LARS Enhanced
#'
#' Creates standardized empty result for failed feature selection
#'
#' @param message Error or warning message
#' @return Empty result list
create_empty_result <- function(message) {
  return(list(
    selected_features = character(0),
    feature_scores = numeric(0),
    method_name = "LARS_Enhanced",
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
