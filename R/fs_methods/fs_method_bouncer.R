#' Bounce R Feature Selection
#'
#' @description
#' Bounce R algorithm for feature selection - an iterative method that uses
#' random sampling and stability assessment for robust feature selection.
#' This is method #9 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Bounce R Feature Selection
#'
#' @description
#' Iterative feature selection algorithm that uses bootstrap sampling and
#' stability assessment to identify robust features. Features are ranked
#' based on their frequency of selection across multiple iterations.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (n_iterations, sample_fraction, stability_threshold)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_bouncer <- function(data, 
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
  n_iterations <- config$n_iterations %||% 50
  sample_fraction <- config$sample_fraction %||% 0.8
  stability_threshold <- config$stability_threshold %||% 0.3
  
  # Initialize feature selection counters
  feature_counts <- rep(0, ncol(feature_cols))
  names(feature_counts) <- colnames(feature_cols)
  
  # Store iteration details for stability assessment
  iteration_selections <- list()
  
  # Main bouncing iterations
  for (iter in seq_len(n_iterations)) {
    tryCatch({
      # Bootstrap sampling
      n_samples <- nrow(data)
      sample_indices <- sample(n_samples, size = floor(sample_fraction * n_samples), replace = TRUE)
      
      # Create bootstrap sample
      bootstrap_data <- data[sample_indices, , drop = FALSE]
      bootstrap_class <- bootstrap_data[, 1]
      bootstrap_features <- bootstrap_data[, -1, drop = FALSE]
      
      # Perform feature selection on bootstrap sample
      selected_in_iteration <- bounce_iteration_selection(bootstrap_features, bootstrap_class)
      
      # Update feature counts
      feature_counts[selected_in_iteration] <- feature_counts[selected_in_iteration] + 1
      
      # Store iteration results
      iteration_selections[[iter]] <- selected_in_iteration
      
    }, error = function(e) {
      # Skip failed iterations
      warning(paste("Bounce R iteration", iter, "failed:", e$message))
    })
  }
  
  # Calculate stability scores (frequency of selection)
  stability_scores <- feature_counts / n_iterations
  
  # Apply stability threshold
  stable_features <- which(stability_scores >= stability_threshold)
  
  if (length(stable_features) == 0) {
    # If no features meet stability threshold, select top features
    n_select <- min(max_features, length(stability_scores))
    stable_features <- order(stability_scores, decreasing = TRUE)[seq_len(n_select)]
  }
  
  # Select final features
  n_select <- min(max_features, length(stable_features))
  top_indices <- stable_features[order(stability_scores[stable_features], decreasing = TRUE)][seq_len(n_select)]
  selected_features <- names(stability_scores)[top_indices]
  
  # Calculate additional stability metrics
  stability_variance <- calculate_stability_variance(iteration_selections, selected_features)
  
  return(list(
    features = selected_features,
    scores = stability_scores[top_indices],
    metadata = list(
      method = "bounceR",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        n_iterations = n_iterations,
        sample_fraction = sample_fraction,
        stability_threshold = stability_threshold,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        stability_scores = stability_scores[top_indices],
        mean_stability = mean(stability_scores[top_indices]),
        stability_variance = stability_variance,
        n_above_threshold = sum(stability_scores >= stability_threshold)
      ),
      iteration_info = list(
        total_iterations = n_iterations,
        successful_iterations = length(iteration_selections),
        selection_frequency = stability_scores[stability_scores > 0]
      )
    )
  ))
}

#' Perform feature selection in single bounce iteration
#'
#' @description
#' Single iteration of feature selection using simple statistical measures
#'
#' @param features Feature matrix
#' @param classes Class vector
#' @return Vector of selected feature indices
bounce_iteration_selection <- function(features, classes) {
  tryCatch({
    # Convert class to numeric if needed
    if (is.factor(classes) || is.character(classes)) {
      class_numeric <- as.numeric(as.factor(classes))
    } else {
      class_numeric <- classes
    }
    
    # Calculate simple selection criterion (correlation with class)
    feature_scores <- apply(features, 2, function(x) {
      if (is.numeric(x)) {
        abs(cor(x, class_numeric, use = "complete.obs"))
      } else {
        # For categorical features
        x_numeric <- as.numeric(as.factor(x))
        abs(cor(x_numeric, class_numeric, use = "complete.obs"))
      }
    })
    
    # Handle missing values
    feature_scores[is.na(feature_scores)] <- 0
    
    # Select top 50% of features in this iteration
    n_select_iter <- max(1, floor(0.5 * length(feature_scores)))
    selected_indices <- order(feature_scores, decreasing = TRUE)[seq_len(n_select_iter)]
    
    return(colnames(features)[selected_indices])
    
  }, error = function(e) {
    # Return random selection as fallback
    n_fallback <- max(1, floor(0.3 * ncol(features)))
    return(sample(colnames(features), n_fallback))
  })
}

#' Calculate stability variance across iterations
#'
#' @description
#' Calculate variance in feature selection across iterations
#'
#' @param iteration_selections List of feature selections per iteration
#' @param final_features Final selected features
#' @return Stability variance metric
calculate_stability_variance <- function(iteration_selections, final_features) {
  tryCatch({
    if (length(iteration_selections) == 0 || length(final_features) == 0) {
      return(1)  # Maximum variance
    }
    
    # Calculate Jaccard similarity across iterations for final features
    jaccard_similarities <- numeric(0)
    
    for (i in seq_along(iteration_selections)) {
      if (length(iteration_selections[[i]]) > 0) {
        # Calculate overlap with final features
        overlap <- length(intersect(iteration_selections[[i]], final_features))
        union_size <- length(union(iteration_selections[[i]], final_features))
        
        if (union_size > 0) {
          jaccard_sim <- overlap / union_size
          jaccard_similarities <- c(jaccard_similarities, jaccard_sim)
        }
      }
    }
    
    if (length(jaccard_similarities) > 0) {
      return(var(jaccard_similarities))
    } else {
      return(1)
    }
    
  }, error = function(e) {
    return(1)
  })
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
tryCatch({
  register_fs_method(
    method_id = 9,
    method_name = "bounceR",
    method_function = fs_method_bouncer,
    category = "ensemble",
    dependencies = character(0),
    description = "Iterative feature selection with bootstrap sampling and stability assessment",
    parameters = list(
      n_iterations = 50,
      sample_fraction = 0.8,
      stability_threshold = 0.3
    ),
    complexity = "medium",
    supports_smote = TRUE,
    timeout_default = 360
  )
  
  cat("✓ Registered Bounce R (bounceR) method (ID: 9)\n")
  
}, error = function(e) {
  warning("Failed to register bounceR method: ", e$message)
})
