#' MDL with SU Feature Selection
#'
#' @description
#' Minimum Description Length principle combined with Symmetrical Uncertainty.
#' Uses information-theoretic measures for feature evaluation and selection.
#' This is method #7 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' MDL with SU Feature Selection
#'
#' @description
#' Feature selection using Minimum Description Length principle with Symmetrical
#' Uncertainty as the evaluation criterion. Balances feature informativeness
#' with model complexity.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (su_threshold, mdl_weight, bins)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_mdl_su <- function(data, 
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
  su_threshold <- config$su_threshold %||% 0.1
  mdl_weight <- config$mdl_weight %||% 0.2
  n_bins <- config$bins %||% 10
  
  # Discretize continuous features for entropy calculation
  discretized_features <- apply(feature_cols, 2, function(x) {
    if (is.numeric(x)) {
      # Discretize using quantile-based binning
      tryCatch({
        cut(x, breaks = quantile(x, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE),
            include.lowest = TRUE, labels = FALSE)
      }, error = function(e) {
        # Fallback to simple binning
        cut(x, breaks = n_bins, labels = FALSE)
      })
    } else {
      # For categorical variables, convert to numeric
      as.numeric(as.factor(x))
    }
  })
  
  # Discretize class variable
  class_discrete <- as.numeric(as.factor(class_col))
  
  # Calculate SU and MDL scores for each feature
  feature_scores <- apply(discretized_features, 2, function(x) {
    tryCatch({
      # Calculate Symmetrical Uncertainty
      su_score <- calculate_symmetrical_uncertainty(x, class_discrete)
      
      # Calculate MDL (complexity penalty)
      n_unique_feature <- length(unique(x[!is.na(x)]))
      n_unique_class <- length(unique(class_discrete))
      n_samples <- length(x)
      
      # MDL approximation: encoding cost for feature-class relationship
      # Complexity penalty based on description length
      complexity_penalty <- (log2(n_unique_feature) + log2(n_unique_class)) / n_samples
      
      # Combined score: SU minus weighted complexity
      combined_score <- su_score - (mdl_weight * complexity_penalty)
      
      return(combined_score)
      
    }, error = function(e) {
      return(0)
    })
  })
  
  # Handle edge cases
  feature_scores[is.na(feature_scores)] <- 0
  feature_scores[is.infinite(feature_scores)] <- 0
  
  # Filter by SU threshold
  valid_features <- which(feature_scores > su_threshold)
  
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
      method = "MDL_SU",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        su_threshold = su_threshold,
        mdl_weight = mdl_weight,
        bins = n_bins,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        combined_scores = feature_scores[top_indices],
        mean_score = mean(feature_scores[top_indices]),
        n_above_threshold = sum(feature_scores > su_threshold)
      ),
      discretization_info = list(
        n_bins_used = n_bins,
        n_classes = length(unique(class_col))
      )
    )
  ))
}

#' Calculate Symmetrical Uncertainty
#'
#' @description
#' Calculate symmetrical uncertainty between two discrete variables
#'
#' @param x First variable (feature)
#' @param y Second variable (class)
#' @return Symmetrical uncertainty value
calculate_symmetrical_uncertainty <- function(x, y) {
  tryCatch({
    # Calculate entropies
    entropy_x <- calculate_entropy(x)
    entropy_y <- calculate_entropy(y)
    joint_entropy <- calculate_joint_entropy(x, y)
    
    # Calculate mutual information
    mutual_info <- entropy_x + entropy_y - joint_entropy
    
    # Calculate symmetrical uncertainty
    if (entropy_x + entropy_y == 0) {
      return(0)
    }
    
    su <- 2 * mutual_info / (entropy_x + entropy_y)
    return(su)
    
  }, error = function(e) {
    return(0)
  })
}

#' Calculate entropy
#'
#' @description
#' Calculate Shannon entropy of a discrete variable
#'
#' @param x Discrete variable
#' @return Entropy value
calculate_entropy <- function(x) {
  tryCatch({
    # Remove missing values
    x_clean <- x[!is.na(x)]
    
    if (length(x_clean) == 0) {
      return(0)
    }
    
    # Calculate probabilities
    prob_table <- table(x_clean)
    probs <- prob_table / sum(prob_table)
    
    # Calculate entropy
    entropy <- -sum(probs * log2(probs + 1e-10))  # Add small constant to avoid log(0)
    return(entropy)
    
  }, error = function(e) {
    return(0)
  })
}

#' Calculate joint entropy
#'
#' @description
#' Calculate joint entropy of two discrete variables
#'
#' @param x First variable
#' @param y Second variable
#' @return Joint entropy value
calculate_joint_entropy <- function(x, y) {
  tryCatch({
    # Remove missing values
    valid_indices <- !is.na(x) & !is.na(y)
    x_clean <- x[valid_indices]
    y_clean <- y[valid_indices]
    
    if (length(x_clean) == 0) {
      return(0)
    }
    
    # Create joint variable
    joint_var <- paste(x_clean, y_clean, sep = "_")
    
    # Calculate joint entropy
    entropy <- calculate_entropy(joint_var)
    return(entropy)
    
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
    method_id = 7,
    method_name = "MDL_SU",
    method_function = fs_method_mdl_su,
    category = "information_theory",
    dependencies = character(0),
    description = "Feature selection using Minimum Description Length with Symmetrical Uncertainty",
    parameters = list(
      su_threshold = 0.1,
      mdl_weight = 0.2,
      bins = 10
    ),
    complexity = "medium",
    supports_smote = TRUE,
    timeout_default = 240
  )
  
  cat("✓ Registered MDL with SU (MDL_SU) method (ID: 7)\n")
  
}, error = function(e) {
  warning("Failed to register MDL_SU method: ", e$message)
})
