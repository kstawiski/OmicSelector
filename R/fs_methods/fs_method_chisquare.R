#' Chi-Square Test Feature Selection
#'
#' @description
#' Select features using chi-square test for independence between each feature
#' and the class variable. Works best with categorical/discrete features.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Chi-Square Feature Selection
#'
#' @description
#' Select features using chi-square test for categorical data or discretized
#' continuous data.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (n_bins for discretization)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_chisquare <- function(data, 
                               config = list(), 
                               max_features = 20, 
                               use_smote = FALSE, 
                               timeout_sec = 300) {
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract classes and features
  classes <- as.factor(data[, 1])
  features <- data[, -1, drop = FALSE]
  
  # Configuration parameters
  n_bins <- config$n_bins %||% 5  # Number of bins for discretization
  
  # Calculate chi-square statistics
  chi_square_stats <- apply(features, 2, function(x) {
    tryCatch({
      # Discretize continuous variables
      if (is.numeric(x)) {
        x_discrete <- cut(x, breaks = n_bins, include.lowest = TRUE)
      } else {
        x_discrete <- as.factor(x)
      }
      
      # Create contingency table
      cont_table <- table(x_discrete, classes)
      
      # Perform chi-square test
      if (any(dim(cont_table) < 2) || any(cont_table < 5)) {
        # Insufficient data for chi-square test
        return(0)
      }
      
      chi_test <- chisq.test(cont_table)
      return(chi_test$statistic)
      
    }, error = function(e) {
      return(0)
    })
  })
  
  # Remove NAs and handle edge cases
  chi_square_stats[is.na(chi_square_stats)] <- 0
  chi_square_stats[is.infinite(chi_square_stats)] <- 0
  
  # Select top features
  if (length(chi_square_stats) == 0) {
    warning("No features could be tested with chi-square")
    return(NULL)
  }
  
  n_select <- min(max_features, length(chi_square_stats))
  selected_indices <- order(chi_square_stats, decreasing = TRUE)[seq_len(n_select)]
  selected_features <- names(features)[selected_indices]
  
  return(list(
    features = selected_features,
    scores = chi_square_stats[selected_indices],
    metadata = list(
      method = "Chi-Square Test",
      n_features_input = ncol(features),
      n_features_selected = length(selected_features),
      parameters_used = list(
        n_bins = n_bins,
        max_features = max_features,
        use_smote = use_smote
      ),
      test_statistics = chi_square_stats[selected_indices]
    )
  ))
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
tryCatch({
  register_fs_method(
    method_id = 101,
    method_name = "Chi-Square Test",
    method_function = fs_method_chisquare,
    category = "statistical",
    dependencies = character(0),
    description = "Select features using chi-square test for independence",
    parameters = list(n_bins = 5),
    complexity = "low",
    supports_smote = TRUE,
    timeout_default = 120
  )
  
  cat("✓ Registered Chi-Square Test feature selection method (ID: 101)\n")
  
}, error = function(e) {
  warning("Failed to register Chi-Square Test method: ", e$message)
})
