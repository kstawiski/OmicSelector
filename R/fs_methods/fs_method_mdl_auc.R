#' MDL with AUC Feature Selection
#'
#' @description
#' Minimum Description Length principle combined with AUC-based feature evaluation.
#' Selects features that provide good discriminative power with minimal complexity.
#' This is method #6 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' MDL with AUC Feature Selection
#'
#' @description
#' Feature selection using Minimum Description Length principle with Area Under
#' the Curve (AUC) as the evaluation criterion for binary classification.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (auc_threshold, mdl_weight)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_mdl_auc <- function(data, 
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
  auc_threshold <- config$auc_threshold %||% 0.6
  mdl_weight <- config$mdl_weight %||% 0.3
  
  # Convert to binary classification if needed
  if (length(unique_classes) > 2) {
    warning("MDL_AUC: Multi-class detected, using one-vs-rest approach")
    # Use first class vs rest
    class_binary <- ifelse(class_col == unique_classes[1], 1, 0)
  } else {
    class_binary <- as.numeric(as.factor(class_col)) - 1
  }
  
  # Calculate AUC and MDL scores for each feature
  feature_scores <- apply(feature_cols, 2, function(x) {
    tryCatch({
      if (is.numeric(x)) {
        # Calculate AUC
        auc_score <- calculate_auc(x, class_binary)
        
        # Calculate MDL (simplified as information content)
        # MDL approximation: -log2(AUC) + complexity penalty
        complexity_penalty <- log2(length(unique(x))) / length(x)
        mdl_score <- -log2(max(auc_score, 0.001)) + complexity_penalty
        
        # Combined score (higher AUC, lower MDL is better)
        combined_score <- auc_score - (mdl_weight * mdl_score)
        
        return(combined_score)
      } else {
        # For categorical features, use contingency table approach
        cont_table <- table(x, class_binary)
        if (any(dim(cont_table) < 2)) {
          return(0)
        }
        
        # Calculate chi-square as proxy for discriminative power
        chi_test <- chisq.test(cont_table)
        chi_score <- chi_test$statistic / sum(cont_table)
        
        # Simple complexity penalty for categorical variables
        complexity_penalty <- log2(length(unique(x))) / length(x)
        combined_score <- chi_score - (mdl_weight * complexity_penalty)
        
        return(combined_score)
      }
    }, error = function(e) {
      return(0)
    })
  })
  
  # Handle edge cases
  feature_scores[is.na(feature_scores)] <- 0
  feature_scores[is.infinite(feature_scores)] <- 0
  
  # Filter by AUC threshold (approximate)
  valid_features <- which(feature_scores > (auc_threshold - mdl_weight))
  
  if (length(valid_features) == 0) {
    # If no features pass threshold, select top features
    n_select <- min(max_features, length(feature_scores))
    valid_features <- order(feature_scores, decreasing = TRUE)[seq_len(n_select)]
  }
  
  # Select top features
  n_select <- min(max_features, length(valid_features))
  top_indices <- valid_features[order(feature_scores[valid_features], decreasing = TRUE)][seq_len(n_select)]
  selected_features <- colnames(feature_cols)[top_indices]
  
  return(list(
    features = selected_features,
    scores = feature_scores[top_indices],
    metadata = list(
      method = "MDL_AUC",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        auc_threshold = auc_threshold,
        mdl_weight = mdl_weight,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        combined_scores = feature_scores[top_indices],
        mean_score = mean(feature_scores[top_indices]),
        n_above_threshold = sum(feature_scores > (auc_threshold - mdl_weight))
      ),
      classification_info = list(
        n_classes = length(unique_classes),
        binary_conversion = length(unique_classes) > 2
      )
    )
  ))
}

#' Calculate AUC for binary classification
#'
#' @description
#' Simple AUC calculation for a single feature
#'
#' @param feature_values Feature values
#' @param class_binary Binary class labels (0/1)
#' @return AUC value
calculate_auc <- function(feature_values, class_binary) {
  tryCatch({
    # Simple AUC calculation using rank statistics
    n_pos <- sum(class_binary == 1)
    n_neg <- sum(class_binary == 0)
    
    if (n_pos == 0 || n_neg == 0) {
      return(0.5)
    }
    
    # Calculate rank-based AUC
    pos_values <- feature_values[class_binary == 1]
    neg_values <- feature_values[class_binary == 0]
    
    # Count concordant pairs
    concordant <- 0
    tied <- 0
    
    for (pos_val in pos_values) {
      concordant <- concordant + sum(pos_val > neg_values)
      tied <- tied + sum(pos_val == neg_values)
    }
    
    auc <- (concordant + 0.5 * tied) / (n_pos * n_neg)
    return(abs(auc - 0.5) + 0.5)  # Ensure AUC >= 0.5
    
  }, error = function(e) {
    return(0.5)
  })
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
tryCatch({
  register_fs_method(
    method_id = 6,
    method_name = "MDL_AUC",
    method_function = fs_method_mdl_auc,
    category = "statistical",
    dependencies = character(0),
    description = "Feature selection using Minimum Description Length with AUC evaluation",
    parameters = list(
      auc_threshold = 0.6,
      mdl_weight = 0.3
    ),
    complexity = "medium",
    supports_smote = TRUE,
    timeout_default = 240
  )
  
  cat("✓ Registered MDL with AUC (MDL_AUC) method (ID: 6)\n")
  
}, error = function(e) {
  warning("Failed to register MDL_AUC method: ", e$message)
})
