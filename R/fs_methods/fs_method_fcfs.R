#' Forward Correlation-based Feature Selection (FCFS)
#'
#' @description
#' Forward stepwise feature selection based on correlation with class variable.
#' Iteratively adds features that improve correlation-based criterion.
#' This is method #5 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Forward Correlation-based Feature Selection
#'
#' @description
#' Forward stepwise selection where features are added based on their
#' correlation with the class variable and low redundancy with already selected features.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (correlation_threshold, redundancy_threshold)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_fcfs <- function(data, 
                          config = list(), 
                          max_features = 20, 
                          use_smote = FALSE, 
                          timeout_sec = 300) {
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract classes and features
  class_col <- data[, 1]
  feature_cols <- data[, -1, drop = FALSE]
  
  # Configuration parameters
  correlation_threshold <- config$correlation_threshold %||% 0.1
  redundancy_threshold <- config$redundancy_threshold %||% 0.8
  
  # Convert class to numeric if factor
  if (is.factor(class_col)) {
    class_numeric <- as.numeric(class_col)
  } else {
    class_numeric <- class_col
  }
  
  # Calculate initial correlations with class
  class_correlations <- apply(feature_cols, 2, function(x) {
    tryCatch({
      if (is.numeric(x)) {
        cor_value <- cor(x, class_numeric, use = "complete.obs")
        return(abs(cor_value))
      } else {
        x_numeric <- as.numeric(as.factor(x))
        cor_value <- cor(x_numeric, class_numeric, use = "complete.obs", method = "spearman")
        return(abs(cor_value))
      }
    }, error = function(e) {
      return(0)
    })
  })
  
  # Handle NAs
  class_correlations[is.na(class_correlations)] <- 0
  class_correlations[is.infinite(class_correlations)] <- 0
  
  # Initialize selection
  selected_features <- character(0)
  remaining_features <- colnames(feature_cols)
  selection_scores <- numeric(0)
  
  # Forward stepwise selection
  for (step in seq_len(max_features)) {
    
    if (length(remaining_features) == 0) {
      break
    }
    
    best_feature <- NULL
    best_score <- -Inf
    
    # Evaluate each remaining feature
    for (feature in remaining_features) {
      
      # Get correlation with class
      class_corr <- class_correlations[feature]
      
      # Skip if correlation is too low
      if (class_corr < correlation_threshold) {
        next
      }
      
      # Calculate redundancy with already selected features
      redundancy_score <- 0
      if (length(selected_features) > 0) {
        feature_values <- feature_cols[[feature]]
        
        for (selected_feature in selected_features) {
          selected_values <- feature_cols[[selected_feature]]
          
          tryCatch({
            if (is.numeric(feature_values) && is.numeric(selected_values)) {
              redundancy <- abs(cor(feature_values, selected_values, use = "complete.obs"))
              redundancy_score <- max(redundancy_score, redundancy, na.rm = TRUE)
            } else {
              # For categorical variables, use Cramér's V approximation
              feature_numeric <- as.numeric(as.factor(feature_values))
              selected_numeric <- as.numeric(as.factor(selected_values))
              redundancy <- abs(cor(feature_numeric, selected_numeric, use = "complete.obs", method = "spearman"))
              redundancy_score <- max(redundancy_score, redundancy, na.rm = TRUE)
            }
          }, error = function(e) {
            # Skip this redundancy calculation
          })
        }
      }
      
      # Skip if too redundant with selected features
      if (redundancy_score > redundancy_threshold) {
        next
      }
      
      # Calculate combined score (high class correlation, low redundancy)
      combined_score <- class_corr - (redundancy_score * 0.5)
      
      if (combined_score > best_score) {
        best_score <- combined_score
        best_feature <- feature
      }
    }
    
    # Add best feature if found
    if (!is.null(best_feature)) {
      selected_features <- c(selected_features, best_feature)
      selection_scores <- c(selection_scores, best_score)
      remaining_features <- setdiff(remaining_features, best_feature)
    } else {
      # No more suitable features
      break
    }
  }
  
  # Return results
  if (length(selected_features) == 0) {
    warning("FCFS selected no features")
    return(NULL)
  }
  
  return(list(
    features = selected_features,
    scores = selection_scores,
    metadata = list(
      method = "fcfs",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        correlation_threshold = correlation_threshold,
        redundancy_threshold = redundancy_threshold,
        max_features = max_features,
        use_smote = use_smote
      ),
      selection_process = list(
        selection_order = selected_features,
        selection_scores = selection_scores,
        n_steps = length(selected_features),
        class_correlations = class_correlations[selected_features]
      ),
      quality_metrics = list(
        mean_class_correlation = mean(class_correlations[selected_features], na.rm = TRUE),
        mean_selection_score = mean(selection_scores, na.rm = TRUE)
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
    method_id = 5,
    method_name = "fcfs",
    method_function = fs_method_fcfs,
    category = "wrapper",
    dependencies = character(0),
    description = "Forward stepwise selection based on correlation with low redundancy",
    parameters = list(
      correlation_threshold = 0.1,
      redundancy_threshold = 0.8
    ),
    complexity = "medium",
    supports_smote = TRUE,
    timeout_default = 240
  )
  
  cat("✓ Registered Forward Correlation-based Feature Selection (fcfs) method (ID: 5)\n")
  
}, error = function(e) {
  warning("Failed to register fcfs method: ", e$message)
})
