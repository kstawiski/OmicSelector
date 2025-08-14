#' Method Template: [METHOD_NAME]
#'
#' @description
#' Template for implementing new feature selection methods in the modular framework.
#' Copy this file and replace the placeholders with your method implementation.
#'
#' Instructions:
#' 1. Copy this file to R/fs_methods/fs_method_[your_method_name].R
#' 2. Replace [METHOD_NAME] with your method name
#' 3. Replace [METHOD_ID] with a unique integer ID (check existing methods first)
#' 4. Implement the fs_method_[your_method_name] function
#' 5. Update the registration call at the bottom
#' 6. The method will be auto-loaded when the framework initializes
#'
#' @author Your Name

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' [METHOD_NAME] Feature Selection
#'
#' @description
#' [Brief description of what your method does and how it works]
#'
#' @param data Training data with class column as first column
#' @param config Configuration list with method-specific parameters
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with components:
#'   - features: Character vector of selected feature names
#'   - scores: Numeric vector of feature scores/importance (optional)
#'   - metadata: List with method execution information
#'   - performance: Optional performance metrics
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data(iris)
#' iris_data <- data.frame(class = iris$Species, iris[1:4])
#' result <- fs_method_[METHOD_NAME_LOWER](iris_data, max_features = 3)
#' print(result$features)
#' }
fs_method_[METHOD_NAME_LOWER] <- function(data, 
                                         config = list(), 
                                         max_features = 20, 
                                         use_smote = FALSE, 
                                         timeout_sec = 300) {
  
  # Check dependencies (if any)
  required_packages <- c("[PACKAGE1]", "[PACKAGE2]")  # Update with your required packages
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(pkg, " package not available")
      return(NULL)
    }
  }
  
  # Apply SMOTE if requested and supported
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract classes and features
  classes <- data[, 1]
  features <- data[, -1, drop = FALSE]
  
  # Get configuration parameters with defaults
  param1 <- config$param1 %||% "default_value1"  # Replace with your parameters
  param2 <- config$param2 %||% "default_value2"
  
  # ==========================================================================
  # IMPLEMENT YOUR METHOD HERE
  # ==========================================================================
  
  # Example placeholder implementation (replace with your method):
  tryCatch({
    
    # Step 1: Implement your feature selection algorithm
    # Example: Calculate some feature scores
    feature_scores <- apply(features, 2, function(x) {
      # Your scoring logic here
      # This is just a placeholder - replace with your method
      var(x, na.rm = TRUE)
    })
    
    # Step 2: Handle missing/invalid scores
    feature_scores[is.na(feature_scores)] <- 0
    feature_scores[is.infinite(feature_scores)] <- 0
    
    # Step 3: Select top features
    if (length(feature_scores) == 0) {
      warning("No features scored by method")
      return(NULL)
    }
    
    n_select <- min(max_features, length(feature_scores))
    selected_indices <- order(feature_scores, decreasing = TRUE)[seq_len(n_select)]
    selected_features <- names(features)[selected_indices]
    
    # Step 4: Prepare results
    result <- list(
      features = selected_features,
      scores = feature_scores[selected_indices],
      metadata = list(
        method = "[METHOD_NAME]",
        n_features_input = ncol(features),
        n_features_selected = length(selected_features),
        parameters_used = list(
          param1 = param1,
          param2 = param2,
          max_features = max_features,
          use_smote = use_smote
        ),
        execution_info = list(
          # Add any method-specific execution information
        )
      )
    )
    
    return(result)
    
  }, error = function(e) {
    warning("[METHOD_NAME] execution failed: ", e$message)
    return(NULL)
  })
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
# This code runs when the file is sourced
tryCatch({
  register_fs_method(
    method_id = [METHOD_ID],  # Replace with unique integer ID (e.g., 101, 102, etc.)
    method_name = "[METHOD_NAME]",  # Replace with descriptive name
    method_function = fs_method_[METHOD_NAME_LOWER],
    category = "[CATEGORY]",  # One of: statistical, wrapper, embedded, ensemble, filter, hybrid
    dependencies = c("[PACKAGE1]", "[PACKAGE2]"),  # Required packages (empty vector if none)
    description = "[Brief description of the method]",
    parameters = list(
      param1 = "default_value1",  # Default parameter values
      param2 = "default_value2"
    ),
    complexity = "[COMPLEXITY]",  # One of: low, medium, high
    supports_smote = TRUE,  # Whether method works with SMOTE-balanced data
    timeout_default = 300  # Default timeout in seconds
  )
  
  cat("✓ Registered [METHOD_NAME] feature selection method (ID: [METHOD_ID])\n")
  
}, error = function(e) {
  warning("Failed to register [METHOD_NAME] method: ", e$message)
})

# =============================================================================
# USAGE NOTES
# =============================================================================

#' Method Categories:
#' - statistical: Statistical tests (t-test, ANOVA, chi-square, etc.)
#' - filter: Filter methods (correlation, mutual information, etc.)
#' - wrapper: Wrapper methods (RFE, forward/backward selection, etc.)
#' - embedded: Embedded methods (LASSO, Ridge, Elastic Net, etc.)
#' - ensemble: Ensemble methods (Boruta, etc.)
#' - hybrid: Hybrid approaches combining multiple strategies
#'
#' Complexity Levels:
#' - low: Fast methods (< 1 minute on typical datasets)
#' - medium: Moderate methods (1-10 minutes on typical datasets)
#' - high: Slow methods (> 10 minutes on typical datasets)
#'
#' Best Practices:
#' 1. Always handle edge cases (empty features, single class, etc.)
#' 2. Use proper error handling with informative messages
#' 3. Validate inputs and provide meaningful defaults
#' 4. Include comprehensive metadata in results
#' 5. Test with different dataset sizes and characteristics
#' 6. Document any special requirements or limitations
#' 7. Use the %||% operator for default parameter handling
#' 8. Return NULL on failure to allow graceful degradation
