#' Orthogonal Matching Pursuit Feature Selection Method
#'
#' Implements Orthogonal Matching Pursuit (OMP) feature selection which iteratively
#' selects features that are most correlated with the current residual and performs
#' orthogonal projection to update residuals.
#'
#' @param ramwas.data Expression data matrix (samples as rows, features as columns)  
#' @param target Target variable vector
#' @param smote_sampling Logical, whether to apply SMOTE sampling
#' @return List containing selected features and method metadata
#' @export
orthogonal_matching_pursuit_fs <- function(ramwas.data, target, smote_sampling = FALSE) {
  
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
    
    # Perform Orthogonal Matching Pursuit feature selection
    omp_result <- perform_orthogonal_matching_pursuit(data_clean, target_numeric)
    
    # Extract selected features
    if (length(omp_result$selected_features) == 0) {
      return(create_empty_result("No features selected by OMP"))
    }
    
    # Calculate timing
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      selected_features = omp_result$selected_features,
      feature_scores = omp_result$feature_scores,
      method_name = "Orthogonal_Matching_Pursuit",
      execution_time = execution_time,
      n_features_selected = length(omp_result$selected_features),
      n_features_input = ncol(ramwas.data),
      method_params = list(
        max_features = omp_result$max_features,
        stopping_criterion = omp_result$stopping_criterion,
        final_residual_norm = omp_result$final_residual_norm,
        n_iterations = omp_result$n_iterations,
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
      method_name = "Orthogonal_Matching_Pursuit",
      execution_time = execution_time,
      n_features_selected = 0,
      n_features_input = ncol(ramwas.data),
      error_message = paste("OMP failed:", e$message),
      method_params = list(smote_sampling = smote_sampling)
    ))
  })
}

#' Perform Orthogonal Matching Pursuit Feature Selection
#'
#' Implements OMP algorithm with automatic stopping criteria
#'
#' @param data_matrix Input data matrix
#' @param target_vector Target variable
#' @return List with selected features and algorithm details
perform_orthogonal_matching_pursuit <- function(data_matrix, target_vector) {
  
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
  initial_target_norm <- sqrt(sum(target_centered^2))
  
  # Initialize OMP variables
  selected_indices <- c()
  residual <- target_centered
  max_iterations <- min(n_features, n_samples - 1, 50)  # Reasonable upper bound
  
  # Stopping criteria
  residual_threshold <- 0.01 * initial_target_norm
  correlation_threshold <- 1e-10
  
  feature_correlations <- numeric(n_features)
  
  # OMP main loop
  for (iteration in 1:max_iterations) {
    
    # Calculate correlations with current residual
    for (j in 1:n_features) {
      if (j %in% selected_indices) {
        feature_correlations[j] <- 0  # Already selected
      } else {
        feature_correlations[j] <- abs(sum(data_scaled[, j] * residual))
      }
    }
    
    # Find feature with maximum correlation
    max_correlation <- max(feature_correlations)
    
    # Check stopping criteria
    if (max_correlation < correlation_threshold) {
      break  # No more meaningful correlations
    }
    
    if (sqrt(sum(residual^2)) < residual_threshold) {
      break  # Residual small enough
    }
    
    # Select feature with maximum correlation
    best_feature <- which.max(feature_correlations)
    selected_indices <- c(selected_indices, best_feature)
    
    # Update residual using orthogonal projection
    X_selected <- data_scaled[, selected_indices, drop = FALSE]
    
    # Solve least squares problem: min ||y - X_selected * beta||^2
    tryCatch({
      # Use QR decomposition for numerical stability
      qr_decomp <- qr(X_selected)
      
      if (qr_decomp$rank == ncol(X_selected)) {
        # Full rank: normal least squares
        coefficients <- qr.coef(qr_decomp, target_centered)
        coefficients[is.na(coefficients)] <- 0
        
        # Update residual
        fitted_values <- X_selected %*% coefficients
        residual <- target_centered - fitted_values
      } else {
        # Rank deficient: use pseudo-inverse
        coefficients <- solve(crossprod(X_selected) + 1e-6 * diag(ncol(X_selected)),
                             crossprod(X_selected, target_centered))
        
        fitted_values <- X_selected %*% coefficients
        residual <- target_centered - fitted_values
      }
    }, error = function(e) {
      # Fallback: simple projection
      last_feature <- data_scaled[, best_feature]
      projection_coeff <- sum(residual * last_feature) / sum(last_feature^2)
      residual <<- residual - projection_coeff * last_feature
    })
    
    # Check if we've explained enough variance
    residual_norm <- sqrt(sum(residual^2))
    explained_variance <- 1 - (residual_norm / initial_target_norm)^2
    
    if (explained_variance > 0.95) {
      break  # Explained 95% of variance
    }
  }
  
  # Calculate final feature scores
  if (length(selected_indices) > 0) {
    # Final least squares to get coefficients
    X_final <- data_scaled[, selected_indices, drop = FALSE]
    
    tryCatch({
      final_coefficients <- solve(crossprod(X_final) + 1e-6 * diag(length(selected_indices)),
                                 crossprod(X_final, target_centered))
    }, error = function(e) {
      final_coefficients <<- rep(1, length(selected_indices))
    })
    
    selected_features <- colnames(data_matrix)[selected_indices]
    feature_scores <- abs(final_coefficients)
    names(feature_scores) <- selected_features
  } else {
    selected_features <- character(0)
    feature_scores <- numeric(0)
  }
  
  return(list(
    selected_features = selected_features,
    feature_scores = feature_scores,
    max_features = max_iterations,
    stopping_criterion = "correlation_threshold",
    final_residual_norm = sqrt(sum(residual^2)),
    n_iterations = length(selected_indices)
  ))
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

#' Create Empty Result for OMP
#'
#' Creates standardized empty result for failed feature selection
#'
#' @param message Error or warning message
#' @return Empty result list
create_empty_result <- function(message) {
  return(list(
    selected_features = character(0),
    feature_scores = numeric(0),
    method_name = "Orthogonal_Matching_Pursuit",
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
