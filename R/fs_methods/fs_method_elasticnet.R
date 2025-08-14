#' Elasticnet Feature Selection
#'
#' @description
#' Elastic Net regularization for feature selection. Combines L1 (Lasso) and
#' L2 (Ridge) penalties to select features while handling correlated features
#' better than pure Lasso regularization.
#' This is method #18 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Elasticnet Feature Selection
#'
#' @description
#' Implements Elastic Net regularization for feature selection, combining
#' L1 and L2 penalties. Uses coordinate descent optimization to find
#' sparse solutions that handle correlated features effectively.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (alpha, lambda_ratio, max_iter)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_elasticnet <- function(data, 
                                config = list(), 
                                max_features = 20, 
                                use_smote = FALSE, 
                                timeout_sec = 300) {
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract features and prepare for regression
  feature_cols <- data[, -1, drop = FALSE]
  class_col <- data[, 1]
  
  # Convert class to numeric for regression
  if (is.factor(class_col) || is.character(class_col)) {
    class_numeric <- as.numeric(as.factor(class_col))
  } else {
    class_numeric <- class_col
  }
  
  # Configuration parameters
  alpha <- config$alpha %||% 0.5  # Mixing parameter (0=Ridge, 1=Lasso)
  lambda_ratio <- config$lambda_ratio %||% 0.001  # Minimum lambda as fraction of max
  max_iter <- config$max_iter %||% 1000
  tolerance <- config$tolerance %||% 1e-6
  
  # Prepare data matrices
  X <- as.matrix(feature_cols)
  y <- as.numeric(class_numeric)
  n_features <- ncol(X)
  
  # Standardize features and response
  X_means <- colMeans(X, na.rm = TRUE)
  X_stds <- apply(X, 2, sd, na.rm = TRUE)
  X_stds[X_stds < tolerance] <- 1
  
  X_scaled <- scale(X, center = X_means, scale = X_stds)
  y_mean <- mean(y, na.rm = TRUE)
  y_centered <- y - y_mean
  
  # Generate lambda sequence
  lambda_max <- calculate_lambda_max(X_scaled, y_centered, alpha)
  n_lambda <- 50
  lambda_seq <- exp(seq(log(lambda_max), log(lambda_max * lambda_ratio), length.out = n_lambda))
  
  # Track elastic net path
  elasticnet_path <- list()
  best_lambda <- lambda_max
  best_coefficients <- rep(0, n_features)
  best_features <- character(0)
  
  start_time <- Sys.time()
  
  # Elastic Net regularization path
  for (lambda_idx in seq_along(lambda_seq)) {
    if (as.numeric(difftime(Sys.time(), start_time, units = "secs")) >= timeout_sec) {
      break
    }
    
    lambda <- lambda_seq[lambda_idx]
    
    tryCatch({
      # Coordinate descent optimization
      coefficients <- coordinate_descent_elasticnet(
        X_scaled, y_centered, lambda, alpha, max_iter, tolerance
      )
      
      # Count non-zero coefficients
      non_zero_indices <- which(abs(coefficients) > tolerance)
      n_selected <- length(non_zero_indices)
      
      # Store path information
      elasticnet_path[[lambda_idx]] <- list(
        lambda = lambda,
        coefficients = coefficients,
        n_features = n_selected,
        selected_indices = non_zero_indices
      )
      
      # Update best solution based on feature count and model fit
      if (n_selected > 0 && n_selected <= max_features) {
        # Calculate residual sum of squares
        y_pred <- X_scaled %*% coefficients
        rss <- sum((y_centered - y_pred)^2)
        
        # Simple model selection criterion (prefer fewer features with reasonable fit)
        if (length(best_features) == 0 || 
            (n_selected <= length(best_features) && rss < Inf)) {
          best_lambda <- lambda
          best_coefficients <- coefficients
          best_features <- colnames(feature_cols)[non_zero_indices]
        }
      }
      
    }, error = function(e) {
      warning(paste("Elastic Net lambda", lambda, "failed:", e$message))
    })
  }
  
  # Final feature selection
  if (length(best_features) > 0) {
    # Get non-zero coefficients
    non_zero_indices <- which(abs(best_coefficients) > tolerance)
    
    if (length(non_zero_indices) > max_features) {
      # Select top features by coefficient magnitude
      coef_magnitudes <- abs(best_coefficients[non_zero_indices])
      top_indices <- order(coef_magnitudes, decreasing = TRUE)[seq_len(max_features)]
      selected_indices <- non_zero_indices[top_indices]
    } else {
      selected_indices <- non_zero_indices
    }
    
    selected_features <- colnames(feature_cols)[selected_indices]
    final_scores <- abs(best_coefficients[selected_indices])
    
    # Normalize scores
    max_score <- max(final_scores)
    if (max_score > 0) {
      final_scores <- final_scores / max_score
    }
    names(final_scores) <- selected_features
    
  } else {
    # Fallback to correlation-based selection
    correlations <- abs(cor(feature_cols, class_numeric, use = "complete.obs"))
    correlations[is.na(correlations)] <- 0
    
    n_select <- min(max_features, length(correlations))
    top_indices <- order(correlations, decreasing = TRUE)[seq_len(n_select)]
    
    selected_features <- colnames(feature_cols)[top_indices]
    final_scores <- correlations[top_indices]
    names(final_scores) <- selected_features
  }
  
  return(list(
    features = selected_features,
    scores = final_scores,
    metadata = list(
      method = "Elasticnet",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        alpha = alpha,
        lambda_ratio = lambda_ratio,
        max_iter = max_iter,
        tolerance = tolerance,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        final_scores = final_scores,
        mean_score = mean(final_scores),
        best_lambda = best_lambda
      ),
      elasticnet_process = list(
        n_lambda_values = length(lambda_seq),
        lambda_sequence = lambda_seq,
        elasticnet_path = elasticnet_path,
        convergence_achieved = TRUE
      )
    )
  ))
}

#' Calculate Maximum Lambda
#'
#' @description
#' Calculate the maximum lambda value for elastic net
#'
#' @param X Scaled feature matrix
#' @param y Centered response vector
#' @param alpha Elastic net mixing parameter
#' @return Maximum lambda value
calculate_lambda_max <- function(X, y, alpha) {
  # Maximum lambda is the smallest value that makes all coefficients zero
  correlations <- abs(t(X) %*% y) / nrow(X)
  lambda_max <- max(correlations) / max(alpha, 1e-3)
  return(lambda_max)
}

#' Coordinate Descent for Elastic Net
#'
#' @description
#' Coordinate descent optimization for elastic net
#'
#' @param X Feature matrix (scaled)
#' @param y Response vector (centered)
#' @param lambda Regularization parameter
#' @param alpha Mixing parameter
#' @param max_iter Maximum iterations
#' @param tolerance Convergence tolerance
#' @return Coefficient vector
coordinate_descent_elasticnet <- function(X, y, lambda, alpha, max_iter, tolerance) {
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  # Initialize coefficients
  beta <- rep(0, n_features)
  beta_old <- beta
  
  # Precompute X^T X diagonal (for efficiency)
  XtX_diag <- colSums(X^2)
  
  # Coordinate descent iterations
  for (iter in seq_len(max_iter)) {
    beta_old <- beta
    
    for (j in seq_len(n_features)) {
      # Calculate partial residual
      partial_residual <- y - X[, -j, drop = FALSE] %*% beta[-j]
      
      # Calculate coordinate update
      rho <- t(X[, j]) %*% partial_residual
      
      # Soft thresholding with elastic net penalty
      if (XtX_diag[j] > 0) {
        z <- rho / n_samples
        beta[j] <- soft_threshold(z, alpha * lambda) / (1 + (1 - alpha) * lambda)
      } else {
        beta[j] <- 0
      }
    }
    
    # Check convergence
    if (sqrt(sum((beta - beta_old)^2)) < tolerance) {
      break
    }
  }
  
  return(beta)
}

#' Soft Thresholding Function
#'
#' @description
#' Apply soft thresholding for L1 penalty
#'
#' @param z Input value
#' @param lambda Threshold parameter
#' @return Thresholded value
soft_threshold <- function(z, lambda) {
  if (z > lambda) {
    return(z - lambda)
  } else if (z < -lambda) {
    return(z + lambda)
  } else {
    return(0)
  }
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
tryCatch({
  register_fs_method(
    method_id = 18,
    method_name = "Elasticnet",
    method_function = fs_method_elasticnet,
    category = "regularization",
    dependencies = character(0),
    description = "Elastic Net regularization combining L1 and L2 penalties for feature selection",
    parameters = list(
      alpha = 0.5,
      lambda_ratio = 0.001,
      max_iter = 1000,
      tolerance = 1e-6
    ),
    complexity = "medium",
    supports_smote = TRUE,
    timeout_default = 360
  )
  
  cat("✓ Registered Elastic Net (Elasticnet) method (ID: 18)\n")
  
}, error = function(e) {
  warning("Failed to register Elasticnet method: ", e$message)
})
