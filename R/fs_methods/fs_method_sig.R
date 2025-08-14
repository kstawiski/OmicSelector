#' Significant Features (t-test/chi-square) Selection
#'
#' @description
#' Select features using statistical significance tests (t-test for continuous
#' features, chi-square for categorical features) with multiple testing correction.
#' This is method #1 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Significant Features Selection
#'
#' @description
#' Select statistically significant features using appropriate tests based on
#' data type and apply Benjamini-Hochberg multiple testing correction.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (p_threshold for significance)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_sig <- function(data, 
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
  p_threshold <- config$p_threshold %||% 0.05
  
  # Perform statistical tests
  p_values <- apply(feature_cols, 2, function(x) {
    tryCatch({
      if (is.numeric(x)) {
        # T-test for continuous features
        if (length(unique(class_col)) == 2) {
          test_result <- t.test(x ~ class_col)
        } else {
          # ANOVA for multi-class
          test_result <- oneway.test(x ~ class_col)
        }
        return(test_result$p.value)
      } else {
        # Chi-square for categorical features
        cont_table <- table(x, class_col)
        if (any(dim(cont_table) < 2) || any(cont_table < 5)) {
          return(1.0)
        }
        test_result <- chisq.test(cont_table)
        return(test_result$p.value)
      }
    }, error = function(e) {
      return(1.0)
    })
  })
  
  # Apply multiple testing correction
  p_adjusted <- p.adjust(p_values, method = "BH")
  
  # Select significant features
  significant_indices <- which(p_adjusted < p_threshold)
  
  if (length(significant_indices) == 0) {
    # If no significant features, return top by p-value
    if (length(p_values) > 0) {
      n_select <- min(max_features, length(p_values))
      significant_indices <- order(p_values)[seq_len(n_select)]
    } else {
      return(NULL)
    }
  }
  
  # Limit to max_features
  if (length(significant_indices) > max_features) {
    top_indices <- order(p_adjusted[significant_indices])[seq_len(max_features)]
    significant_indices <- significant_indices[top_indices]
  }
  
  selected_features <- colnames(feature_cols)[significant_indices]
  
  return(list(
    features = selected_features,
    scores = 1 - p_adjusted[significant_indices],  # Convert p-values to scores
    metadata = list(
      method = "sig",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        p_threshold = p_threshold,
        max_features = max_features,
        use_smote = use_smote,
        correction_method = "BH"
      ),
      test_statistics = list(
        p_values_raw = p_values[significant_indices],
        p_values_adjusted = p_adjusted[significant_indices],
        n_significant = sum(p_adjusted < p_threshold)
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
    method_id = 1,
    method_name = "sig",
    method_function = fs_method_sig,
    category = "statistical",
    dependencies = character(0),
    description = "Select features using statistical significance tests with multiple testing correction",
    parameters = list(p_threshold = 0.05),
    complexity = "low",
    supports_smote = TRUE,
    timeout_default = 120
  )
  
  cat("✓ Registered Significant Features (sig) method (ID: 1)\n")
  
}, error = function(e) {
  warning("Failed to register sig method: ", e$message)
})
