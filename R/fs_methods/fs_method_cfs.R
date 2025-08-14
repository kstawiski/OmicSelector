#' Correlation-based Feature Selection (CFS)
#'
#' @description
#' Select features based on correlation with class variable. Features with
#' highest absolute correlation to the class are selected.
#' This is method #3 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Correlation-based Feature Selection
#'
#' @description
#' Select features with highest correlation to the class variable. Handles
#' both numeric and categorical class variables.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (correlation method)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_cfs <- function(data, 
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
  correlation_method <- config$correlation_method %||% "pearson"
  
  # Convert class to numeric if factor
  if (is.factor(class_col)) {
    class_numeric <- as.numeric(class_col)
  } else {
    class_numeric <- class_col
  }
  
  # Calculate correlations with class
  correlations <- apply(feature_cols, 2, function(x) {
    tryCatch({
      if (is.numeric(x)) {
        # Use specified correlation method
        cor_value <- cor(x, class_numeric, use = "complete.obs", method = correlation_method)
        return(abs(cor_value))
      } else if (is.factor(x) || is.character(x)) {
        # For categorical features, use Cramér's V or convert to numeric
        x_numeric <- as.numeric(as.factor(x))
        cor_value <- cor(x_numeric, class_numeric, use = "complete.obs", method = "spearman")
        return(abs(cor_value))
      } else {
        return(0)
      }
    }, error = function(e) {
      return(0)
    })
  })
  
  # Handle NAs
  correlations[is.na(correlations)] <- 0
  correlations[is.infinite(correlations)] <- 0
  
  # Select top correlated features
  if (length(correlations) == 0) {
    return(NULL)
  }
  
  n_select <- min(max_features, length(correlations))
  top_indices <- order(correlations, decreasing = TRUE)[seq_len(n_select)]
  selected_features <- colnames(feature_cols)[top_indices]
  
  return(list(
    features = selected_features,
    scores = correlations[top_indices],
    metadata = list(
      method = "cfs",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        correlation_method = correlation_method,
        max_features = max_features,
        use_smote = use_smote
      ),
      correlation_statistics = list(
        correlations = correlations[top_indices],
        mean_correlation = mean(correlations[top_indices], na.rm = TRUE),
        max_correlation = max(correlations[top_indices], na.rm = TRUE),
        min_correlation = min(correlations[top_indices], na.rm = TRUE)
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
    method_id = 3,
    method_name = "cfs",
    method_function = fs_method_cfs,
    category = "filter",
    dependencies = character(0),
    description = "Select features based on correlation with class variable",
    parameters = list(correlation_method = "pearson"),
    complexity = "low",
    supports_smote = TRUE,
    timeout_default = 120
  )
  
  cat("✓ Registered Correlation-based Feature Selection (cfs) method (ID: 3)\n")
  
}, error = function(e) {
  warning("Failed to register cfs method: ", e$message)
})
