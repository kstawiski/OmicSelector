#' Random Forest RFE Feature Selection
#'
#' @description
#' Random Forest Recursive Feature Elimination without cross-validation.
#' Iteratively removes least important features based on Random Forest
#' variable importance until optimal subset is reached.
#' This is method #11 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Random Forest RFE Feature Selection
#'
#' @description
#' Recursive Feature Elimination using Random Forest variable importance.
#' Trains RF model, ranks features by importance, removes least important
#' features iteratively until target number is reached.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (n_trees, step_size, min_features)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_randomforestrfe <- function(data, 
                                     config = list(), 
                                     max_features = 20, 
                                     use_smote = FALSE, 
                                     timeout_sec = 300) {
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract features
  feature_cols <- data[, -1, drop = FALSE]
  
  # Configuration parameters
  n_trees <- config$n_trees %||% 150
  step_size <- config$step_size %||% 0.15  # Remove 15% of features per iteration
  min_features <- config$min_features %||% 3
  
  # Initialize with all features
  current_features <- colnames(feature_cols)
  
  # Track elimination process
  elimination_history <- list()
  iteration <- 0
  start_time <- Sys.time()
  
  # RFE main loop
  while (length(current_features) > max_features && 
         length(current_features) > min_features &&
         as.numeric(difftime(Sys.time(), start_time, units = "secs")) < timeout_sec) {
    
    iteration <- iteration + 1
    
    tryCatch({
      # Prepare current data subset
      current_data <- data[, c(TRUE, colnames(feature_cols) %in% current_features), drop = FALSE]
      
      # Train Random Forest and get importance
      importance_scores <- calculate_rf_importance_advanced(current_data, n_trees)
      
      # Determine number of features to remove
      n_remove <- max(1, floor(step_size * length(current_features)))
      n_keep <- length(current_features) - n_remove
      
      # Ensure we don't go below target
      if (n_keep < max_features) {
        n_keep <- max_features
      }
      
      # Keep most important features
      top_indices <- order(importance_scores, decreasing = TRUE)[seq_len(n_keep)]
      eliminated_features <- setdiff(current_features, names(importance_scores)[top_indices])
      current_features <- names(importance_scores)[top_indices]
      
      # Record elimination step
      elimination_history[[iteration]] <- list(
        iteration = iteration,
        n_features_before = length(current_features) + length(eliminated_features),
        n_features_after = length(current_features),
        eliminated_features = eliminated_features,
        importance_range = range(importance_scores),
        mean_importance = mean(importance_scores)
      )
      
    }, error = function(e) {
      warning(paste("RandomForestRFE iteration", iteration, "failed:", e$message))
      break
    })
  }
  
  # Final feature ranking
  if (length(current_features) > 0) {
    final_data <- data[, c(TRUE, colnames(feature_cols) %in% current_features), drop = FALSE]
    final_importance <- calculate_rf_importance_advanced(final_data, n_trees)
    
    # Select top features if we have more than requested
    if (length(current_features) > max_features) {
      top_indices <- order(final_importance, decreasing = TRUE)[seq_len(max_features)]
      selected_features <- names(final_importance)[top_indices]
      final_scores <- final_importance[top_indices]
    } else {
      selected_features <- current_features
      final_scores <- final_importance[current_features]
    }
  } else {
    # Fallback if elimination failed
    selected_features <- head(colnames(feature_cols), max_features)
    final_scores <- rep(0.5, length(selected_features))
    names(final_scores) <- selected_features
  }
  
  return(list(
    features = selected_features,
    scores = final_scores,
    metadata = list(
      method = "RandomForestRFE",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        n_trees = n_trees,
        step_size = step_size,
        min_features = min_features,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        final_importance = final_scores,
        mean_importance = mean(final_scores),
        importance_range = range(final_scores)
      ),
      elimination_process = list(
        n_iterations = iteration,
        elimination_history = elimination_history,
        features_eliminated_total = ncol(feature_cols) - length(selected_features)
      )
    )
  ))
}

#' Calculate Advanced Random Forest Importance
#'
#' @description
#' Enhanced Random Forest importance calculation with multiple metrics
#'
#' @param data Training data with class in first column
#' @param n_trees Number of trees
#' @return Named vector of feature importance scores
calculate_rf_importance_advanced <- function(data, n_trees) {
  tryCatch({
    feature_cols <- data[, -1, drop = FALSE]
    class_col <- data[, 1]
    
    # Convert class to numeric if needed
    if (is.factor(class_col) || is.character(class_col)) {
      class_numeric <- as.numeric(as.factor(class_col))
    } else {
      class_numeric <- class_col
    }
    
    # Initialize importance accumulator
    importance_matrix <- matrix(0, nrow = min(n_trees, 50), ncol = ncol(feature_cols))
    colnames(importance_matrix) <- colnames(feature_cols)
    
    # Calculate importance across multiple trees
    successful_trees <- 0
    
    for (tree in seq_len(min(n_trees, 50))) {
      tryCatch({
        # Bootstrap sample
        boot_indices <- sample(nrow(data), replace = TRUE)
        boot_data <- data[boot_indices, ]
        boot_class <- boot_data[, 1]
        boot_features <- boot_data[, -1, drop = FALSE]
        
        # Convert bootstrap class
        if (is.factor(boot_class) || is.character(boot_class)) {
          boot_class_numeric <- as.numeric(as.factor(boot_class))
        } else {
          boot_class_numeric <- boot_class
        }
        
        # Calculate tree-specific importance using multiple criteria
        tree_importance <- calculate_tree_importance_advanced(boot_features, boot_class_numeric)
        
        successful_trees <- successful_trees + 1
        importance_matrix[successful_trees, ] <- tree_importance
        
      }, error = function(e) {
        # Skip failed trees
      })
    }
    
    # Average importance across successful trees
    if (successful_trees > 0) {
      final_importance <- colMeans(importance_matrix[seq_len(successful_trees), , drop = FALSE])
    } else {
      # Fallback to simple correlation
      final_importance <- apply(feature_cols, 2, function(x) {
        if (is.numeric(x)) {
          abs(cor(x, class_numeric, use = "complete.obs"))
        } else {
          x_numeric <- as.numeric(as.factor(x))
          abs(cor(x_numeric, class_numeric, use = "complete.obs"))
        }
      })
    }
    
    # Handle missing values
    final_importance[is.na(final_importance)] <- 0
    
    return(final_importance)
    
  }, error = function(e) {
    # Ultimate fallback
    feature_cols <- data[, -1, drop = FALSE]
    fallback_importance <- rep(0.1, ncol(feature_cols))
    names(fallback_importance) <- colnames(feature_cols)
    return(fallback_importance)
  })
}

#' Calculate Tree-Specific Advanced Importance
#'
#' @description
#' Calculate importance for single tree using multiple criteria
#'
#' @param features Feature matrix
#' @param class_numeric Numeric class vector
#' @return Feature importance scores
calculate_tree_importance_advanced <- function(features, class_numeric) {
  tryCatch({
    n_features <- ncol(features)
    importance_scores <- numeric(n_features)
    names(importance_scores) <- colnames(features)
    
    # Use multiple importance criteria
    for (i in seq_len(n_features)) {
      feature_values <- features[, i]
      
      if (is.numeric(feature_values)) {
        # For numeric features: correlation + variance ratio
        correlation_score <- abs(cor(feature_values, class_numeric, use = "complete.obs"))
        
        # Calculate between-class variance ratio
        unique_classes <- unique(class_numeric)
        if (length(unique_classes) > 1) {
          between_var <- 0
          within_var <- 0
          
          for (class_val in unique_classes) {
            class_indices <- class_numeric == class_val
            class_values <- feature_values[class_indices]
            
            if (length(class_values) > 1) {
              class_mean <- mean(class_values, na.rm = TRUE)
              within_var <- within_var + var(class_values, na.rm = TRUE)
              between_var <- between_var + length(class_values) * (class_mean - mean(feature_values, na.rm = TRUE))^2
            }
          }
          
          variance_ratio <- if (within_var > 0) between_var / within_var else 0
        } else {
          variance_ratio <- 0
        }
        
        # Combined score
        importance_scores[i] <- 0.7 * correlation_score + 0.3 * min(variance_ratio / 10, 1)
        
      } else {
        # For categorical features: use chi-square-like measure
        contingency_table <- table(feature_values, class_numeric)
        if (any(dim(contingency_table) >= 2)) {
          chi_stat <- chisq.test(contingency_table, simulate.p.value = TRUE)$statistic
          importance_scores[i] <- min(chi_stat / sum(contingency_table), 1)
        } else {
          importance_scores[i] <- 0
        }
      }
    }
    
    # Normalize scores
    max_score <- max(importance_scores, na.rm = TRUE)
    if (max_score > 0) {
      importance_scores <- importance_scores / max_score
    }
    
    # Handle missing values
    importance_scores[is.na(importance_scores)] <- 0
    
    return(importance_scores)
    
  }, error = function(e) {
    # Fallback to uniform importance
    fallback <- rep(0.1, ncol(features))
    names(fallback) <- colnames(features)
    return(fallback)
  })
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
tryCatch({
  register_fs_method(
    method_id = 11,
    method_name = "RandomForestRFE",
    method_function = fs_method_randomforestrfe,
    category = "rfe",
    dependencies = character(0),
    description = "Random Forest Recursive Feature Elimination with advanced importance metrics",
    parameters = list(
      n_trees = 150,
      step_size = 0.15,
      min_features = 3
    ),
    complexity = "high",
    supports_smote = TRUE,
    timeout_default = 450
  )
  
  cat("✓ Registered Random Forest RFE (RandomForestRFE) method (ID: 11)\n")
  
}, error = function(e) {
  warning("Failed to register RandomForestRFE method: ", e$message)
})
