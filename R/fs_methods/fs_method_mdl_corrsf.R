#' MDL with Correlation-based Feature Selection
#'
#' @description
#' Minimum Description Length principle combined with correlation-based feature
#' selection. Balances feature relevance with redundancy control.
#' This is method #8 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' MDL with Correlation-based Feature Selection
#'
#' @description
#' Feature selection using Minimum Description Length principle combined with
#' correlation analysis. Seeks features with high class correlation but low
#' inter-feature correlation while minimizing model complexity.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (corr_threshold, redundancy_threshold, mdl_weight)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_mdl_corrsf <- function(data, 
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
  corr_threshold <- config$corr_threshold %||% 0.3
  redundancy_threshold <- config$redundancy_threshold %||% 0.8
  mdl_weight <- config$mdl_weight %||% 0.2
  
  # Convert class to numeric for correlation
  if (is.factor(class_col) || is.character(class_col)) {
    class_numeric <- as.numeric(as.factor(class_col))
  } else {
    class_numeric <- class_col
  }
  
  # Calculate class correlation for each feature
  class_correlations <- apply(feature_cols, 2, function(x) {
    tryCatch({
      if (is.numeric(x)) {
        abs(cor(x, class_numeric, use = "complete.obs"))
      } else {
        # For categorical features, use point-biserial correlation approximation
        x_numeric <- as.numeric(as.factor(x))
        abs(cor(x_numeric, class_numeric, use = "complete.obs"))
      }
    }, error = function(e) {
      return(0)
    })
  })
  
  # Handle missing correlations
  class_correlations[is.na(class_correlations)] <- 0
  
  # Filter features by correlation threshold
  relevant_features <- which(class_correlations >= corr_threshold)
  
  if (length(relevant_features) == 0) {
    # If no features pass threshold, use top correlated features
    n_select <- min(max_features, length(class_correlations))
    relevant_features <- order(class_correlations, decreasing = TRUE)[seq_len(n_select)]
  }
  
  # Calculate inter-feature correlations for redundancy control
  selected_features <- integer(0)
  remaining_features <- relevant_features
  
  while (length(selected_features) < max_features && length(remaining_features) > 0) {
    if (length(selected_features) == 0) {
      # Select feature with highest class correlation
      best_idx <- which.max(class_correlations[remaining_features])
      selected_features <- c(selected_features, remaining_features[best_idx])
      remaining_features <- remaining_features[-best_idx]
    } else {
      # Calculate MDL-based scores for remaining features
      mdl_scores <- sapply(remaining_features, function(idx) {
        feature_data <- feature_cols[, idx]
        
        # Calculate class correlation benefit
        class_corr_benefit <- class_correlations[idx]
        
        # Calculate redundancy penalty
        selected_data <- feature_cols[, selected_features, drop = FALSE]
        redundancy_penalty <- calculate_feature_redundancy(feature_data, selected_data)
        
        # Calculate MDL complexity penalty
        if (is.numeric(feature_data)) {
          complexity_penalty <- log2(length(unique(feature_data[!is.na(feature_data)]))) / length(feature_data)
        } else {
          complexity_penalty <- log2(length(levels(as.factor(feature_data)))) / length(feature_data)
        }
        
        # Combined MDL score
        mdl_score <- class_corr_benefit - redundancy_penalty - (mdl_weight * complexity_penalty)
        return(mdl_score)
      })
      
      # Select feature with best MDL score
      if (length(mdl_scores) > 0 && max(mdl_scores) > 0) {
        best_idx <- which.max(mdl_scores)
        selected_features <- c(selected_features, remaining_features[best_idx])
        remaining_features <- remaining_features[-best_idx]
      } else {
        break
      }
    }
  }
  
  # Get final selected feature names
  selected_feature_names <- colnames(feature_cols)[selected_features]
  
  # Calculate final scores
  final_scores <- class_correlations[selected_features]
  
  return(list(
    features = selected_feature_names,
    scores = final_scores,
    metadata = list(
      method = "MDL_CorrSF",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_feature_names),
      parameters_used = list(
        corr_threshold = corr_threshold,
        redundancy_threshold = redundancy_threshold,
        mdl_weight = mdl_weight,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        class_correlations = final_scores,
        mean_correlation = mean(final_scores),
        n_above_threshold = sum(class_correlations >= corr_threshold)
      ),
      selection_process = list(
        n_relevant_features = length(relevant_features),
        selection_iterations = length(selected_features)
      )
    )
  ))
}

#' Calculate Feature Redundancy
#'
#' @description
#' Calculate redundancy of a feature with respect to already selected features
#'
#' @param feature_data New feature data
#' @param selected_data Already selected features data
#' @return Redundancy penalty score
calculate_feature_redundancy <- function(feature_data, selected_data) {
  tryCatch({
    if (ncol(selected_data) == 0) {
      return(0)
    }
    
    # Calculate correlation with each selected feature
    correlations <- apply(selected_data, 2, function(x) {
      if (is.numeric(feature_data) && is.numeric(x)) {
        abs(cor(feature_data, x, use = "complete.obs"))
      } else {
        # For categorical features, use contingency table association
        feature_numeric <- as.numeric(as.factor(feature_data))
        x_numeric <- as.numeric(as.factor(x))
        abs(cor(feature_numeric, x_numeric, use = "complete.obs"))
      }
    })
    
    # Handle missing correlations
    correlations[is.na(correlations)] <- 0
    
    # Return maximum correlation as redundancy penalty
    max_redundancy <- max(correlations)
    return(max_redundancy)
    
  }, error = function(e) {
    return(0)
  })
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
tryCatch({
  register_fs_method(
    method_id = 8,
    method_name = "MDL_CorrSF",
    method_function = fs_method_mdl_corrsf,
    category = "statistical",
    dependencies = character(0),
    description = "Feature selection using MDL with correlation-based selection and redundancy control",
    parameters = list(
      corr_threshold = 0.3,
      redundancy_threshold = 0.8,
      mdl_weight = 0.2
    ),
    complexity = "medium",
    supports_smote = TRUE,
    timeout_default = 300
  )
  
  cat("✓ Registered MDL with Correlation SF (MDL_CorrSF) method (ID: 8)\n")
  
}, error = function(e) {
  warning("Failed to register MDL_CorrSF method: ", e$message)
})
