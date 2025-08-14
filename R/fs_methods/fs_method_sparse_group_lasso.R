#' Sparse Group LASSO Feature Selection Method
#'
#' Implements Sparse Group LASSO feature selection which combines group sparsity
#' with individual feature sparsity. Features are grouped and the method can
#' select entire groups while also selecting individual features within groups.
#'
#' @param ramwas.data Expression data matrix (samples as rows, features as columns)  
#' @param target Target variable vector
#' @param smote_sampling Logical, whether to apply SMOTE sampling
#' @return List containing selected features and method metadata
#' @export
sparse_group_lasso_fs <- function(ramwas.data, target, smote_sampling = FALSE) {
  
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
    
    # Perform Sparse Group LASSO feature selection
    sgl_result <- perform_sparse_group_lasso(data_clean, target_numeric)
    
    # Extract selected features
    if (length(sgl_result$selected_features) == 0) {
      return(create_empty_result("No features selected by Sparse Group LASSO"))
    }
    
    # Calculate timing
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      selected_features = sgl_result$selected_features,
      feature_scores = sgl_result$feature_scores,
      method_name = "Sparse_Group_LASSO",
      execution_time = execution_time,
      n_features_selected = length(sgl_result$selected_features),
      n_features_input = ncol(ramwas.data),
      method_params = list(
        lambda = sgl_result$lambda_optimal,
        alpha = sgl_result$alpha_optimal,
        cv_folds = sgl_result$cv_folds,
        n_groups = sgl_result$n_groups,
        group_method = "correlation_clustering",
        smote_sampling = smote_sampling
      ),
      preprocessing_info = list(
        constant_features_removed = preprocessed$constant_removed,
        correlated_features_removed = preprocessed$correlated_removed,
        final_features_count = ncol(data_clean)
      ),
      group_info = sgl_result$group_info
    ))
    
  }, error = function(e) {
    # Return error result
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      selected_features = character(0),
      feature_scores = numeric(0),
      method_name = "Sparse_Group_LASSO",
      execution_time = execution_time,
      n_features_selected = 0,
      n_features_input = ncol(ramwas.data),
      error_message = paste("Sparse Group LASSO failed:", e$message),
      method_params = list(smote_sampling = smote_sampling)
    ))
  })
}

#' Perform Sparse Group LASSO Feature Selection
#'
#' Implements Sparse Group LASSO with dual penalties for group and individual sparsity
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with selected features and optimal parameters
perform_sparse_group_lasso <- function(data_matrix, target_vector) {
  
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
  group_info <- form_correlation_groups_sgl(data_scaled)
  
  # Generate lambda and alpha sequences
  lambda_max <- calculate_lambda_max_sgl(data_scaled, target_centered, group_info)
  lambda_sequence <- exp(seq(log(lambda_max), log(lambda_max * 0.001), length.out = 30))
  alpha_sequence <- c(0.1, 0.3, 0.5, 0.7, 0.9)  # Balance between group and individual sparsity
  
  # Grid search with cross-validation
  cv_folds <- min(5, n_samples)
  fold_assignments <- sample(rep(1:cv_folds, length.out = n_samples))
  
  best_error <- Inf
  best_lambda <- lambda_sequence[1]
  best_alpha <- alpha_sequence[1]
  
  for (alpha in alpha_sequence) {
    for (lambda in lambda_sequence) {
      
      fold_errors <- numeric(cv_folds)
      
      for (fold in 1:cv_folds) {
        train_idx <- which(fold_assignments != fold)
        test_idx <- which(fold_assignments == fold)
        
        if (length(train_idx) < 2 || length(test_idx) < 1) next
        
        # Fit Sparse Group LASSO on training set
        coefficients <- sparse_group_lasso_solver(
          data_scaled[train_idx, , drop = FALSE], 
          target_centered[train_idx], 
          group_info,
          lambda,
          alpha
        )
        
        # Predict on test set
        predictions <- data_scaled[test_idx, , drop = FALSE] %*% coefficients
        fold_errors[fold] <- mean((target_centered[test_idx] - predictions)^2)
      }
      
      avg_error <- mean(fold_errors[fold_errors > 0])
      if (avg_error < best_error) {
        best_error <- avg_error
        best_lambda <- lambda
        best_alpha <- alpha
      }
    }
  }
  
  # Fit final model with optimal parameters
  final_coefficients <- sparse_group_lasso_solver(
    data_scaled, target_centered, group_info, best_lambda, best_alpha
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
    lambda_optimal = best_lambda,
    alpha_optimal = best_alpha,
    cv_folds = cv_folds,
    coefficients = final_coefficients,
    n_groups = length(unique(group_info$group_assignments)),
    group_info = group_info
  ))
}

#' Form Correlation-Based Groups for Sparse Group LASSO
#'
#' Creates feature groups based on correlation clustering
#'
#' @param data_matrix Standardized data matrix
#' @param cor_threshold Correlation threshold for grouping
#' @return List with group assignments and sizes
form_correlation_groups_sgl <- function(data_matrix, cor_threshold = 0.6) {
  
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

#' Sparse Group LASSO Solver
#'
#' Implements coordinate descent for Sparse Group LASSO with dual penalties
#'
#' @param X Design matrix (standardized)
#' @param y Response vector (centered)
#' @param group_info Group information structure
#' @param lambda Overall regularization parameter
#' @param alpha Balance between group (alpha=0) and individual (alpha=1) sparsity
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return Coefficient vector
sparse_group_lasso_solver <- function(X, y, group_info, lambda, alpha, max_iter = 1000, tol = 1e-6) {
  
  n_features <- ncol(X)
  group_assignments <- group_info$group_assignments
  unique_groups <- unique(group_assignments)
  
  # Initialize coefficients
  beta <- numeric(n_features)
  beta_old <- beta
  
  # Individual and group penalty weights
  lambda_individual <- lambda * alpha
  lambda_group <- lambda * (1 - alpha)
  
  # Coordinate descent iterations
  for (iter in 1:max_iter) {
    
    # Update each group
    for (g in unique_groups) {
      group_indices <- which(group_assignments == g)
      group_size <- length(group_indices)
      
      if (group_size == 0) next
      
      # Extract group submatrix
      X_group <- X[, group_indices, drop = FALSE]
      
      # Calculate residual for this group
      residual <- y - X[, -group_indices, drop = FALSE] %*% beta[-group_indices]
      
      # Update group coefficients with dual penalties
      if (group_size == 1) {
        # Single feature: use individual penalty only
        j <- group_indices[1]
        residual_j <- as.numeric(crossprod(X_group, residual))
        XtX_diag <- sum(X_group^2)
        
        if (XtX_diag > 0) {
          beta[j] <- soft_threshold_sgl(residual_j, lambda_individual) / XtX_diag
        } else {
          beta[j] <- 0
        }
      } else {
        # Multiple features: use both penalties
        group_coefficients <- update_group_coefficients_sgl(
          X_group, residual, lambda_individual, lambda_group, group_size
        )
        beta[group_indices] <- group_coefficients
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

#' Update Group Coefficients for Sparse Group LASSO
#'
#' Updates coefficients for a group with dual penalties
#'
#' @param X_group Group design matrix
#' @param residual Current residual vector
#' @param lambda_individual Individual penalty parameter
#' @param lambda_group Group penalty parameter
#' @param group_size Number of features in group
#' @return Updated group coefficients
update_group_coefficients_sgl <- function(X_group, residual, lambda_individual, lambda_group, group_size) {
  
  # Calculate group gradient
  group_gradient <- as.numeric(crossprod(X_group, residual))
  
  # Step 1: Apply individual soft thresholding
  beta_individual <- soft_threshold_sgl(group_gradient, lambda_individual)
  
  # Step 2: Apply group soft thresholding
  beta_norm <- sqrt(sum(beta_individual^2))
  group_threshold <- lambda_group * sqrt(group_size)
  
  if (beta_norm > group_threshold) {
    # Group survives: apply group shrinkage
    shrinkage_factor <- (beta_norm - group_threshold) / beta_norm
    beta_group <- beta_individual * shrinkage_factor
    
    # Step 3: Solve local quadratic problem
    XtX_group <- crossprod(X_group)
    
    tryCatch({
      # Solve for final coefficients
      final_coefficients <- solve(XtX_group + 1e-6 * diag(group_size), beta_group)
      return(final_coefficients)
    }, error = function(e) {
      # Fallback to simpler update
      diag_XtX <- diag(XtX_group)
      diag_XtX[diag_XtX == 0] <- 1e-6
      return(beta_group / diag_XtX)
    })
  } else {
    # Group eliminated: set to zero
    return(numeric(group_size))
  }
}

#' Soft Thresholding Function for Sparse Group LASSO
#'
#' Applies soft thresholding operator
#'
#' @param x Input vector or value
#' @param lambda Threshold parameter
#' @return Soft thresholded result
soft_threshold_sgl <- function(x, lambda) {
  return(sign(x) * pmax(0, abs(x) - lambda))
}

#' Calculate Maximum Lambda for Sparse Group LASSO
#'
#' Computes the maximum lambda value for regularization path
#'
#' @param X Design matrix
#' @param y Response vector
#' @param group_info Group information structure
#' @return Maximum lambda value
calculate_lambda_max_sgl <- function(X, y, group_info) {
  
  # Individual penalty contribution
  individual_max <- max(abs(crossprod(X, y))) / nrow(X)
  
  # Group penalty contribution
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
  
  group_max <- max_group_norm / nrow(X)
  
  # Return maximum of both contributions
  return(max(individual_max, group_max))
}

#' Basic Feature Preprocessing
#'
#' Removes constant and highly correlated features
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @param cor_threshold Correlation threshold for removal
#' @return List with processed data and removal info
preprocess_features_basic <- function(data_matrix, target_vector, cor_threshold = 0.98) {
  
  # Remove constant features
  feature_vars <- apply(data_matrix, 2, var, na.rm = TRUE)
  constant_features <- which(feature_vars < 1e-10)
  
  if (length(constant_features) > 0) {
    data_matrix <- data_matrix[, -constant_features, drop = FALSE]
  }
  
  # Remove highly correlated features (slightly relaxed for SGL)
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

#' Create Empty Result for Sparse Group LASSO
#'
#' Creates standardized empty result for failed feature selection
#'
#' @param message Error or warning message
#' @return Empty result list
create_empty_result <- function(message) {
  return(list(
    selected_features = character(0),
    feature_scores = numeric(0),
    method_name = "Sparse_Group_LASSO",
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
