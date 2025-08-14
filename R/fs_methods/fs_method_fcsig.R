#' Fold Change + Significance Feature Selection
#'
#' @description
#' Select features using combination of fold change and statistical significance.
#' First finds significant features, then ranks by fold change between classes.
#' This is method #2 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Fold Change + Significance Selection
#'
#' @description
#' Two-step feature selection: first identify statistically significant features,
#' then rank by fold change magnitude between classes.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (p_threshold, fc_threshold)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_fcsig <- function(data, 
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
  unique_classes <- unique(class_col)
  
  # Configuration parameters
  p_threshold <- config$p_threshold %||% 0.05
  fc_threshold <- config$fc_threshold %||% 1.5
  
  # Step 1: Find significant features
  p_values <- apply(feature_cols, 2, function(x) {
    tryCatch({
      if (is.numeric(x)) {
        if (length(unique_classes) == 2) {
          test_result <- t.test(x ~ class_col)
        } else {
          test_result <- oneway.test(x ~ class_col)
        }
        return(test_result$p.value)
      } else {
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
  significant_features <- colnames(feature_cols)[p_adjusted < p_threshold]
  
  if (length(significant_features) == 0) {
    # If no significant features, use top by p-value
    n_candidates <- min(max_features * 2, length(p_values))
    significant_features <- colnames(feature_cols)[order(p_values)[seq_len(n_candidates)]]
  }
  
  # Step 2: Calculate fold changes for significant features
  if (length(unique_classes) == 2) {
    fold_changes <- sapply(significant_features, function(feature_name) {
      if (feature_name %in% colnames(feature_cols)) {
        feature_values <- feature_cols[[feature_name]]
        
        if (is.numeric(feature_values)) {
          class1_values <- feature_values[class_col == unique_classes[1]]
          class2_values <- feature_values[class_col == unique_classes[2]]
          
          mean1 <- mean(class1_values, na.rm = TRUE)
          mean2 <- mean(class2_values, na.rm = TRUE)
          
          # Calculate fold change (avoiding division by zero)
          if (mean2 != 0 && mean1 > 0 && mean2 > 0) {
            return(abs(log2(mean1 / mean2)))
          } else {
            return(0)
          }
        } else {
          return(0)
        }
      } else {
        return(0)
      }
    })
  } else {
    # For multi-class, use coefficient of variation as proxy
    fold_changes <- sapply(significant_features, function(feature_name) {
      if (feature_name %in% colnames(feature_cols)) {
        feature_values <- feature_cols[[feature_name]]
        if (is.numeric(feature_values)) {
          class_means <- tapply(feature_values, class_col, mean, na.rm = TRUE)
          return(sd(class_means, na.rm = TRUE) / mean(class_means, na.rm = TRUE))
        } else {
          return(0)
        }
      } else {
        return(0)
      }
    })
  }
  
  # Step 3: Select features with highest fold changes
  if (length(fold_changes) == 0) {
    return(NULL)
  }
  
  n_select <- min(max_features, length(fold_changes))
  top_indices <- order(fold_changes, decreasing = TRUE)[seq_len(n_select)]
  selected_features <- significant_features[top_indices]
  
  return(list(
    features = selected_features,
    scores = fold_changes[top_indices],
    metadata = list(
      method = "fcsig",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        p_threshold = p_threshold,
        fc_threshold = fc_threshold,
        max_features = max_features,
        use_smote = use_smote
      ),
      analysis_steps = list(
        n_significant_features = length(significant_features),
        fold_changes = fold_changes[top_indices],
        p_values = p_adjusted[match(selected_features, colnames(feature_cols))]
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
    method_id = 2,
    method_name = "fcsig",
    method_function = fs_method_fcsig,
    category = "statistical",
    dependencies = character(0),
    description = "Select features using combination of fold change and statistical significance",
    parameters = list(
      p_threshold = 0.05,
      fc_threshold = 1.5
    ),
    complexity = "low",
    supports_smote = TRUE,
    timeout_default = 180
  )
  
  cat("✓ Registered Fold Change + Significance (fcsig) method (ID: 2)\n")
  
}, error = function(e) {
  warning("Failed to register fcsig method: ", e$message)
})
