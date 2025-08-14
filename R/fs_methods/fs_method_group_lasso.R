#' Group LASSO Feature Selection Method
#'
#' Implements Group LASSO feature selection which performs variable selection
#' at the group level. Features are automatically grouped based on correlation
#' structure, and entire groups are selected or excluded together.
#'
#' @param ramwas.data Expression data matrix (samples as rows, features as columns)  
#' @param target Target variable vector
#' @param smote_sampling Logical, whether to apply SMOTE sampling
#' @return List containing selected features and method metadata
#' @export
group_lasso_fs <- function(ramwas.data, target, smote_sampling = FALSE) {
  
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
    
    # Perform Group LASSO feature selection
    group_lasso_result <- perform_group_lasso(data_clean, target_numeric)
    
    # Extract selected features
    if (length(group_lasso_result$selected_features) == 0) {
      return(create_empty_result("No features selected by Group LASSO"))
    }
    
    # Calculate timing
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      selected_features = group_lasso_result$selected_features,
      feature_scores = group_lasso_result$feature_scores,
      method_name = "Group_LASSO",
      execution_time = execution_time,
      n_features_selected = length(group_lasso_result$selected_features),
      n_features_input = ncol(ramwas.data),
      method_params = list(
        lambda = group_lasso_result$lambda_optimal,
        n_groups = group_lasso_result$n_groups,
        cv_folds = group_lasso_result$cv_folds,
        group_method = "correlation_clustering",
        smote_sampling = smote_sampling
      ),
      preprocessing_info = list(
        constant_features_removed = preprocessed$constant_removed,
        correlated_features_removed = preprocessed$correlated_removed,
        final_features_count = ncol(data_clean)
      ),
      group_info = group_lasso_result$group_info
    ))
    
  }, error = function(e) {
    # Return error result
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      selected_features = character(0),
      feature_scores = numeric(0),
      method_name = "Group_LASSO",
      execution_time = execution_time,
      n_features_selected = 0,
      n_features_input = ncol(ramwas.data),
      error_message = paste("Group LASSO failed:", e$message),
      method_params = list(smote_sampling = smote_sampling)
    ))
  })
}

#' Perform Group LASSO Feature Selection
#'
#' Implements Group LASSO with automatic group formation and optimization
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with selected features and group information
perform_group_lasso <- function(data_matrix, target_vector) {
  
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
  
  # Form groups based on correlation structure
  group_info <- form_correlation_groups(data_scaled)
  
  # Generate lambda sequence
  lambda_max <- calculate_lambda_max_group_lasso(data_scaled, target_centered, group_info)
  lambda_sequence <- exp(seq(log(lambda_max), log(lambda_max * 0.001), length.out = 40))
  
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
      
      # Fit Group LASSO on training set
      coefficients <- group_lasso_solver(
        data_scaled[train_idx, , drop = FALSE], 
        target_centered[train_idx], 
        group_info, 
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
  final_coefficients <- group_lasso_solver(data_scaled, target_centered, group_info, lambda_optimal)
  
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
    n_groups = length(unique(group_info$group_assignments)),
    group_info = group_info
  ))
}

#' Form Correlation-Based Groups
#'
#' Creates feature groups based on correlation clustering
#'
#' @param data_matrix Standardized data matrix
#' @param cor_threshold Correlation threshold for grouping
#' @return List with group assignments and sizes
form_correlation_groups <- function(data_matrix, cor_threshold = 0.7) {
  
  n_features <- ncol(data_matrix)
  
  # Compute correlation matrix
  cor_matrix <- cor(data_matrix, use = "complete.obs")
  cor_matrix[is.na(cor_matrix)] <- 0
  
  # Initialize group assignments
  group_assignments <- 1:n_features
  
  # Simple correlation clustering
  for (i in 1:(n_features - 1)) {
    if (group_assignments[i] != i) next  # Already assigned to a group
    
    for (j in (i + 1):n_features) {
      if (abs(cor_matrix[i, j]) > cor_threshold) {
        group_assignments[j] <- group_assignments[i]
      }
    }
  }
  
  # Renumber groups consecutively
  unique_groups <- unique(group_assignments)
  group_mapping <- setNames(seq_along(unique_groups), unique_groups)
  group_assignments <- group_mapping[as.character(group_assignments)]
  
  # Calculate group sizes
  group_sizes <- table(group_assignments)
  
  return(list(
    group_assignments = group_assignments,
    group_sizes = as.numeric(group_sizes),
    n_groups = length(unique_groups)
  ))
}

#' Group LASSO Solver
#'
#' Implements block coordinate descent for Group LASSO
#'
#' @param X Design matrix (standardized)
#' @param y Response vector (centered)
#' @param group_info Group information structure
#' @param lambda Regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return Coefficient vector
group_lasso_solver <- function(X, y, group_info, lambda, max_iter = 1000, tol = 1e-6) {
  
  n_features <- ncol(X)
  group_assignments <- group_info$group_assignments
  unique_groups <- unique(group_assignments)
  
  # Initialize coefficients
  beta <- numeric(n_features)
  beta_old <- beta
  
  # Block coordinate descent iterations
  for (iter in 1:max_iter) {
    
    for (g in unique_groups) {
      group_indices <- which(group_assignments == g)
      group_size <- length(group_indices)
      
      if (group_size == 0) next
      
      # Extract group submatrix
      X_group <- X[, group_indices, drop = FALSE]
      
      # Calculate residual for this group
      residual <- y - X[, -group_indices, drop = FALSE] %*% beta[-group_indices]
      group_gradient <- as.numeric(crossprod(X_group, residual))
      
      # Group soft thresholding
      gradient_norm <- sqrt(sum(group_gradient^2))
      
      if (gradient_norm > 0) {
        # Calculate group Hessian approximation (diagonal)
        XtX_group <- crossprod(X_group)
        
        # Group threshold
        group_threshold <- lambda * sqrt(group_size)
        
        if (gradient_norm > group_threshold) {
          # Update group coefficients
          shrinkage_factor <- (gradient_norm - group_threshold) / gradient_norm
          
          # Solve for group coefficients
          group_coefficients <- tryCatch({
            solve(XtX_group + 1e-6 * diag(group_size), group_gradient * shrinkage_factor)
          }, error = function(e) {
            # Fallback to simple shrinkage
            group_gradient * shrinkage_factor / (sum(diag(XtX_group)) + 1e-6)
          })
          
          beta[group_indices] <- group_coefficients
        } else {
          # Shrink to zero
          beta[group_indices] <- 0
        }
      } else {
        beta[group_indices] <- 0
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

#' Calculate Maximum Lambda for Group LASSO
#'
#' Computes the maximum lambda value for Group LASSO regularization path
#'
#' @param X Design matrix
#' @param y Response vector
#' @param group_info Group information structure
#' @return Maximum lambda value
calculate_lambda_max_group_lasso <- function(X, y, group_info) {
  
  group_assignments <- group_info$group_assignments
  unique_groups <- unique(group_assignments)
  
  max_group_norm <- 0
  
  for (g in unique_groups) {
    group_indices <- which(group_assignments == g)
    group_size <- length(group_indices)
    
    if (group_size > 0) {
      X_group <- X[, group_indices, drop = FALSE]
      group_gradient <- as.numeric(crossprod(X_group, y))
      group_norm <- sqrt(sum(group_gradient^2)) / sqrt(group_size)
      max_group_norm <- max(max_group_norm, group_norm)
    }
  }
  
  return(max_group_norm / nrow(X))
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
  
  # Remove highly correlated features (relaxed threshold for Group LASSO)
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

#' Create Empty Result for Group LASSO
#'
#' Creates standardized empty result for failed feature selection
#'
#' @param message Error or warning message
#' @return Empty result list
create_empty_result <- function(message) {
  return(list(
    selected_features = character(0),
    feature_scores = numeric(0),
    method_name = "Group_LASSO",
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
