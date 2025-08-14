#' LARS Feature Selection
#'
#' @description
#' Least Angle Regression (LARS) for feature selection. Efficiently finds
#' a sequence of feature additions that maximally correlate with the residual,
#' providing a regularization path for feature selection.
#' This is method #17 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' LARS Feature Selection
#'
#' @description
#' Implements Least Angle Regression algorithm for feature selection.
#' Builds a sequence of models by adding features that are most correlated
#' with the current residual, providing an efficient regularization path.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (max_steps, normalize, eps)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_lars <- function(data, 
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
  max_steps <- config$max_steps %||% min(max_features, ncol(feature_cols))
  normalize <- config$normalize %||% TRUE
  eps <- config$eps %||% 1e-10
  
  # Prepare data matrices
  X <- as.matrix(feature_cols)
  y <- as.numeric(class_numeric)
  n_features <- ncol(X)
  
  # Normalize features if requested
  if (normalize) {
    X_means <- colMeans(X, na.rm = TRUE)
    X_stds <- apply(X, 2, sd, na.rm = TRUE)
    X_stds[X_stds < eps] <- 1  # Avoid division by zero
    
    X <- scale(X, center = X_means, scale = X_stds)
    y <- scale(y, center = TRUE, scale = FALSE)
  }
  
  # Initialize LARS algorithm
  active_set <- integer(0)  # Indices of active features
  inactive_set <- seq_len(n_features)  # Indices of inactive features
  
  beta <- rep(0, n_features)  # Coefficient vector
  residual <- y  # Current residual
  
  # Track LARS path
  lars_path <- list()
  step <- 0
  start_time <- Sys.time()
  
  # Main LARS loop
  while (length(active_set) < max_steps && 
         length(inactive_set) > 0 &&
         as.numeric(difftime(Sys.time(), start_time, units = "secs")) < timeout_sec) {
    
    step <- step + 1
    
    tryCatch({
      # Calculate correlations with residual
      correlations <- abs(t(X[, inactive_set, drop = FALSE]) %*% residual)
      
      if (length(correlations) == 0) break
      
      # Find feature with maximum correlation
      max_corr_idx <- which.max(correlations)
      new_feature <- inactive_set[max_corr_idx]
      max_correlation <- correlations[max_corr_idx]
      
      # Add to active set
      active_set <- c(active_set, new_feature)
      inactive_set <- setdiff(inactive_set, new_feature)
      
      # Calculate direction vector using equiangular constraint
      X_active <- X[, active_set, drop = FALSE]
      
      if (length(active_set) == 1) {
        # First step: direction is just the sign of correlation
        direction <- sign(t(X_active) %*% residual)
        gamma <- max_correlation / (t(X_active) %*% X_active)
      } else {
        # Multiple features: solve equiangular system
        G <- t(X_active) %*% X_active  # Gram matrix
        ones <- rep(1, length(active_set))
        
        # Solve for equiangular direction
        tryCatch({
          G_inv <- solve(G + diag(eps, nrow(G)))
          AA <- as.numeric(t(ones) %*% G_inv %*% ones)
          w <- G_inv %*% ones / sqrt(AA)
          u <- X_active %*% w
          
          # Calculate step size
          gamma <- max_correlation / AA
          
          # Check for direction constraints
          if (length(inactive_set) > 0) {
            X_inactive <- X[, inactive_set, drop = FALSE]
            c_inactive <- t(X_inactive) %*% residual
            a_inactive <- t(X_inactive) %*% u
            
            # Find minimum positive step that maintains equiangular property
            valid_steps <- (max_correlation - c_inactive) / (AA - a_inactive)
            valid_steps <- valid_steps[valid_steps > eps]
            
            if (length(valid_steps) > 0) {
              gamma <- min(gamma, min(valid_steps))
            }
          }
          
          direction <- w
          
        }, error = function(e) {
          # Fallback to simple direction
          direction <- rep(0, length(active_set))
          direction[length(direction)] <- 1
        })
      }
      
      # Update coefficients
      beta[active_set] <- beta[active_set] + gamma * direction
      
      # Update residual
      residual <- y - X %*% beta
      
      # Store path information
      lars_path[[step]] <- list(
        step = step,
        feature_added = new_feature,
        active_features = active_set,
        coefficients = beta[beta != 0],
        residual_norm = sqrt(sum(residual^2)),
        correlation = max_correlation
      )
      
    }, error = function(e) {
      warning(paste("LARS step", step, "failed:", e$message))
      break
    })
  }
  
  # Select final features
  if (length(active_set) > 0) {
    # Rank features by order of entry and coefficient magnitude
    feature_importance <- abs(beta[active_set])
    names(feature_importance) <- colnames(feature_cols)[active_set]
    
    # Select top features
    n_select <- min(max_features, length(active_set))
    top_indices <- order(feature_importance, decreasing = TRUE)[seq_len(n_select)]
    
    selected_features <- names(feature_importance)[top_indices]
    final_scores <- feature_importance[top_indices]
    
    # Normalize scores to [0, 1]
    max_score <- max(final_scores)
    if (max_score > 0) {
      final_scores <- final_scores / max_score
    }
    
  } else {
    # Fallback if LARS failed
    correlations_fallback <- abs(cor(feature_cols, class_numeric, use = "complete.obs"))
    correlations_fallback[is.na(correlations_fallback)] <- 0
    
    n_select <- min(max_features, length(correlations_fallback))
    top_indices <- order(correlations_fallback, decreasing = TRUE)[seq_len(n_select)]
    
    selected_features <- colnames(feature_cols)[top_indices]
    final_scores <- correlations_fallback[top_indices]
    names(final_scores) <- selected_features
  }
  
  return(list(
    features = selected_features,
    scores = final_scores,
    metadata = list(
      method = "LARS",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        max_steps = max_steps,
        normalize = normalize,
        eps = eps,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        final_scores = final_scores,
        mean_score = mean(final_scores),
        score_range = range(final_scores)
      ),
      lars_process = list(
        steps_completed = step,
        lars_path = lars_path,
        final_active_set = active_set,
        final_residual_norm = if(step > 0) lars_path[[step]]$residual_norm else NA
      )
    )
  ))
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
tryCatch({
  register_fs_method(
    method_id = 17,
    method_name = "LARS",
    method_function = fs_method_lars,
    category = "regularization",
    dependencies = character(0),
    description = "Least Angle Regression for efficient feature selection with regularization path",
    parameters = list(
      max_steps = 20,
      normalize = TRUE,
      eps = 1e-10
    ),
    complexity = "medium",
    supports_smote = TRUE,
    timeout_default = 300
  )
  
  cat("✓ Registered LARS (LARS) method (ID: 17)\n")
  
}, error = function(e) {
  warning("Failed to register LARS method: ", e$message)
})
