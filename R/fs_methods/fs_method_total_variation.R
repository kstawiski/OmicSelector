#' Total Variation Regularization Feature Selection Method
#'
#' Implements Total Variation (TV) regularization for feature selection.
#' This method promotes piecewise constant solutions and is particularly
#' useful when features have spatial or sequential relationships.
#'
#' @param ramwas.data Expression data matrix (samples as rows, features as columns)  
#' @param target Target variable vector
#' @param smote_sampling Logical, whether to apply SMOTE sampling
#' @return List containing selected features and method metadata
#' @export
total_variation_fs <- function(ramwas.data, target, smote_sampling = FALSE) {
  
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
    
    # Perform Total Variation feature selection
    tv_result <- perform_total_variation_selection(data_clean, target_numeric)
    
    # Extract selected features
    if (length(tv_result$selected_features) == 0) {
      return(create_empty_result("No features selected by Total Variation"))
    }
    
    # Calculate timing
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      selected_features = tv_result$selected_features,
      feature_scores = tv_result$feature_scores,
      method_name = "Total_Variation",
      execution_time = execution_time,
      n_features_selected = length(tv_result$selected_features),
      n_features_input = ncol(ramwas.data),
      method_params = list(
        lambda = tv_result$lambda_optimal,
        cv_folds = tv_result$cv_folds,
        tv_type = "discrete_difference",
        regularization_strength = tv_result$regularization_strength,
        smote_sampling = smote_sampling
      ),
      preprocessing_info = list(
        constant_features_removed = preprocessed$constant_removed,
        correlated_features_removed = preprocessed$correlated_removed,
        final_features_count = ncol(data_clean)
      ),
      tv_info = tv_result$tv_info
    ))
    
  }, error = function(e) {
    # Return error result
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      selected_features = character(0),
      feature_scores = numeric(0),
      method_name = "Total_Variation",
      execution_time = execution_time,
      n_features_selected = 0,
      n_features_input = ncol(ramwas.data),
      error_message = paste("Total Variation failed:", e$message),
      method_params = list(smote_sampling = smote_sampling)
    ))
  })
}

#' Perform Total Variation Feature Selection
#'
#' Implements Total Variation regularization with iterative optimization
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with selected features and TV information
perform_total_variation_selection <- function(data_matrix, target_vector) {
  
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
  
  # Create difference matrix for Total Variation penalty
  D_matrix <- create_difference_matrix(n_features)
  
  # Generate lambda sequence for TV regularization
  lambda_max <- calculate_lambda_max_tv(data_scaled, target_centered, D_matrix)
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
      
      # Fit TV model on training set
      coefficients <- total_variation_solver(
        data_scaled[train_idx, , drop = FALSE], 
        target_centered[train_idx], 
        D_matrix,
        lambda
      )
      
      # Predict on test set
      predictions <- data_scaled[test_idx, , drop = FALSE] %*% coefficients
      fold_errors[fold] <- mean((target_centered[test_idx] - predictions)^2)
    }
    
    cv_errors[i] <- mean(fold_errors[fold_errors > 0])
  }
  
  # Select optimal lambda
  optimal_idx <- which.min(cv_errors)
  lambda_optimal <- lambda_sequence[optimal_idx]
  
  # Fit final model with optimal lambda
  final_coefficients <- total_variation_solver(data_scaled, target_centered, D_matrix, lambda_optimal)
  
  # Calculate total variation and feature importance
  tv_norm <- calculate_tv_norm(final_coefficients, D_matrix)
  feature_importance <- abs(final_coefficients)
  
  # Apply additional TV-based feature ranking
  tv_based_importance <- calculate_tv_feature_importance(final_coefficients, D_matrix)
  combined_importance <- 0.7 * feature_importance + 0.3 * tv_based_importance
  
  # Select features based on combined importance
  # Default: select features with non-zero importance
  selected_idx <- which(combined_importance > 1e-6)
  
  if (length(selected_idx) == 0) {
    # Fallback: select top 15% of features
    n_select <- max(1, round(0.15 * n_features))
    selected_idx <- order(combined_importance, decreasing = TRUE)[1:n_select]
  }
  
  selected_features <- colnames(data_matrix)[selected_idx]
  feature_scores <- combined_importance[selected_idx]
  names(feature_scores) <- selected_features
  
  return(list(
    selected_features = selected_features,
    feature_scores = feature_scores,
    lambda_optimal = lambda_optimal,
    cv_folds = cv_folds,
    regularization_strength = tv_norm,
    tv_info = list(
      total_variation_norm = tv_norm,
      n_iterations = "convergence_based"
    )
  ))
}

#' Create Difference Matrix for Total Variation
#'
#' Creates difference matrix for discrete Total Variation penalty
#'
#' @param n_features Number of features
#' @return Difference matrix
create_difference_matrix <- function(n_features) {
  
  if (n_features < 2) {
    return(matrix(0, nrow = 0, ncol = n_features))
  }
  
  # Create first-order difference matrix
  D_matrix <- matrix(0, nrow = n_features - 1, ncol = n_features)
  
  for (i in 1:(n_features - 1)) {
    D_matrix[i, i] <- -1
    D_matrix[i, i + 1] <- 1
  }
  
  return(D_matrix)
}

#' Total Variation Solver
#'
#' Solves Total Variation regularized regression using iterative methods
#'
#' @param X Design matrix
#' @param y Response vector
#' @param D Difference matrix
#' @param lambda TV penalty parameter
#' @return Coefficient vector
total_variation_solver <- function(X, y, D, lambda) {
  
  n_features <- ncol(X)
  
  # Initial solution using ridge regression
  beta_init <- tryCatch({
    solve(crossprod(X) + 1e-6 * diag(n_features), crossprod(X, y))
  }, error = function(e) {
    rep(0, n_features)
  })
  
  # Iterative TV optimization using split Bregman or ADMM approach
  beta <- tv_admm_solver(X, y, D, lambda, beta_init)
  
  return(beta)
}

#' ADMM Solver for Total Variation
#'
#' Implements ADMM algorithm for TV regularized regression
#'
#' @param X Design matrix
#' @param y Response vector
#' @param D Difference matrix
#' @param lambda TV penalty parameter
#' @param beta_init Initial coefficient vector
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return Optimized coefficient vector
tv_admm_solver <- function(X, y, D, lambda, beta_init, max_iter = 100, tol = 1e-6) {
  
  n_features <- ncol(X)
  n_diff <- nrow(D)
  
  # Initialize ADMM variables
  beta <- beta_init
  z <- D %*% beta  # Auxiliary variable
  u <- numeric(n_diff)  # Dual variable
  
  # ADMM penalty parameter
  rho <- 1.0
  
  # Precompute matrices for efficiency
  XtX <- crossprod(X)
  Xty <- as.numeric(crossprod(X, y))
  DtD <- crossprod(D)
  
  # ADMM iterations
  for (iter in 1:max_iter) {
    
    beta_old <- beta
    
    # Beta update: solve (X'X + rho*D'D) beta = X'y + rho*D'(z - u)
    rhs <- Xty + rho * as.numeric(crossprod(D, z - u))
    lhs_matrix <- XtX + rho * DtD
    
    beta <- tryCatch({
      solve(lhs_matrix + 1e-6 * diag(n_features), rhs)
    }, error = function(e) {
      # Fallback: use diagonal approximation
      diag_lhs <- diag(lhs_matrix)
      diag_lhs[diag_lhs == 0] <- 1e-6
      rhs / diag_lhs
    })
    
    # Z update: soft thresholding for L1 penalty
    Dbeta_u <- D %*% beta + u
    z <- soft_threshold_tv(Dbeta_u, lambda / rho)
    
    # U update: dual variable update
    u <- u + D %*% beta - z
    
    # Check convergence
    primal_residual <- sqrt(sum((D %*% beta - z)^2))
    dual_residual <- sqrt(sum((beta - beta_old)^2))
    
    if (primal_residual < tol && dual_residual < tol) {
      break
    }
    
    # Adaptive rho update
    if (primal_residual > 10 * dual_residual) {
      rho <- 2 * rho
    } else if (dual_residual > 10 * primal_residual) {
      rho <- rho / 2
    }
  }
  
  return(beta)
}

#' Soft Thresholding for Total Variation
#'
#' Applies soft thresholding operator for TV penalty
#'
#' @param x Input vector
#' @param lambda Threshold parameter
#' @return Soft thresholded vector
soft_threshold_tv <- function(x, lambda) {
  return(sign(x) * pmax(0, abs(x) - lambda))
}

#' Calculate TV Norm
#'
#' Computes the total variation norm of coefficient vector
#'
#' @param beta Coefficient vector
#' @param D Difference matrix
#' @return Total variation norm
calculate_tv_norm <- function(beta, D) {
  
  if (nrow(D) == 0) return(0)
  
  differences <- D %*% beta
  tv_norm <- sum(abs(differences))
  
  return(tv_norm)
}

#' Calculate TV-based Feature Importance
#'
#' Computes feature importance based on total variation structure
#'
#' @param beta Coefficient vector
#' @param D Difference matrix
#' @return Feature importance vector
calculate_tv_feature_importance <- function(beta, D) {
  
  n_features <- length(beta)
  
  if (nrow(D) == 0) {
    return(abs(beta))
  }
  
  # Calculate local variation for each feature
  differences <- abs(D %*% beta)
  feature_importance <- numeric(n_features)
  
  # Each feature contributes to differences with its neighbors
  for (i in 1:n_features) {
    if (i == 1) {
      # First feature: only right difference
      if (length(differences) >= 1) {
        feature_importance[i] <- differences[1]
      }
    } else if (i == n_features) {
      # Last feature: only left difference
      if (length(differences) >= i - 1) {
        feature_importance[i] <- differences[i - 1]
      }
    } else {
      # Middle features: average of left and right differences
      if (length(differences) >= i) {
        feature_importance[i] <- mean(c(differences[i - 1], differences[i]))
      }
    }
  }
  
  return(feature_importance)
}

#' Calculate Maximum Lambda for Total Variation
#'
#' Computes the maximum lambda value for TV regularization
#'
#' @param X Design matrix
#' @param y Response vector
#' @param D Difference matrix
#' @return Maximum lambda value
calculate_lambda_max_tv <- function(X, y, D) {
  
  # Calculate gradient of least squares objective
  beta_ols <- if (ncol(X) <= nrow(X)) {
    tryCatch({
      solve(crossprod(X) + 1e-6 * diag(ncol(X)), crossprod(X, y))
    }, error = function(e) {
      rep(0, ncol(X))
    })
  } else {
    rep(0, ncol(X))
  }
  
  # Maximum TV penalty based on differences
  if (nrow(D) > 0) {
    tv_gradient <- D %*% beta_ols
    lambda_max <- max(abs(tv_gradient))
  } else {
    lambda_max <- max(abs(crossprod(X, y))) / nrow(X)
  }
  
  return(lambda_max)
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

#' Create Empty Result for Total Variation
#'
#' Creates standardized empty result for failed feature selection
#'
#' @param message Error or warning message
#' @return Empty result list
create_empty_result <- function(message) {
  return(list(
    selected_features = character(0),
    feature_scores = numeric(0),
    method_name = "Total_Variation",
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
