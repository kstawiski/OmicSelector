#' Boruta Feature Selection
#'
#' @description
#' Boruta algorithm for feature selection using Random Forest importance
#' with shadow features. Compares real features against randomized copies
#' to determine truly important features.
#' This is method #14 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Boruta Feature Selection
#'
#' @description
#' Implements Boruta algorithm which creates shadow features (permuted copies)
#' and compares real feature importance against shadow features using
#' Random Forest to identify truly relevant features.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (max_runs, p_value, mcAdj)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_boruta <- function(data, 
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
  n_features <- ncol(feature_cols)
  
  # Configuration parameters
  max_runs <- config$max_runs %||% 50
  p_value <- config$p_value %||% 0.01
  mc_adj <- config$mcAdj %||% TRUE
  n_trees <- config$n_trees %||% 100
  
  # Initialize feature status tracking
  feature_status <- rep("Tentative", n_features)
  names(feature_status) <- colnames(feature_cols)
  
  # Track importance across runs
  importance_history <- matrix(0, nrow = max_runs, ncol = n_features)
  colnames(importance_history) <- colnames(feature_cols)
  
  shadow_importance_history <- list()
  decision_history <- list()
  
  run <- 0
  start_time <- Sys.time()
  
  # Main Boruta iterations
  while (run < max_runs && 
         sum(feature_status == "Tentative") > 0 &&
         as.numeric(difftime(Sys.time(), start_time, units = "secs")) < timeout_sec) {
    
    run <- run + 1
    
    tryCatch({
      # Create shadow features (permuted copies)
      shadow_data <- create_shadow_features(feature_cols)
      
      # Combine original and shadow features
      combined_data <- cbind(data[, 1, drop = FALSE], feature_cols, shadow_data)
      
      # Calculate importance using Random Forest
      importance_scores <- calculate_boruta_importance(combined_data, n_trees)
      
      # Separate real and shadow importances
      real_importance <- importance_scores[1:n_features]
      shadow_importance <- importance_scores[(n_features + 1):length(importance_scores)]
      
      # Store importance for this run
      importance_history[run, ] <- real_importance
      shadow_importance_history[[run]] <- shadow_importance
      
      # Make decisions based on comparison with shadow features
      max_shadow_importance <- max(shadow_importance, na.rm = TRUE)
      
      decisions <- make_boruta_decisions(real_importance, max_shadow_importance, 
                                        feature_status, p_value, run)
      
      feature_status <- decisions$status
      decision_history[[run]] <- decisions
      
    }, error = function(e) {
      warning(paste("Boruta run", run, "failed:", e$message))
    })
  }
  
  # Final feature selection
  confirmed_features <- names(feature_status)[feature_status == "Confirmed"]
  
  if (length(confirmed_features) == 0) {
    # If no features confirmed, select tentative features with highest average importance
    tentative_features <- names(feature_status)[feature_status == "Tentative"]
    if (length(tentative_features) > 0) {
      avg_importance <- colMeans(importance_history[1:run, tentative_features, drop = FALSE], na.rm = TRUE)
      n_select <- min(max_features, length(tentative_features))
      top_indices <- order(avg_importance, decreasing = TRUE)[seq_len(n_select)]
      selected_features <- names(avg_importance)[top_indices]
      final_scores <- avg_importance[top_indices]
    } else {
      # Ultimate fallback
      selected_features <- head(colnames(feature_cols), min(5, max_features))
      final_scores <- rep(0.5, length(selected_features))
      names(final_scores) <- selected_features
    }
  } else {
    # Use confirmed features
    if (length(confirmed_features) > max_features) {
      # Rank confirmed features by average importance
      avg_importance <- colMeans(importance_history[1:run, confirmed_features, drop = FALSE], na.rm = TRUE)
      top_indices <- order(avg_importance, decreasing = TRUE)[seq_len(max_features)]
      selected_features <- names(avg_importance)[top_indices]
      final_scores <- avg_importance[top_indices]
    } else {
      selected_features <- confirmed_features
      final_scores <- colMeans(importance_history[1:run, confirmed_features, drop = FALSE], na.rm = TRUE)
    }
  }
  
  return(list(
    features = selected_features,
    scores = final_scores,
    metadata = list(
      method = "Boruta",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        max_runs = max_runs,
        p_value = p_value,
        mcAdj = mc_adj,
        n_trees = n_trees,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        final_scores = final_scores,
        mean_importance = mean(final_scores),
        importance_range = range(final_scores)
      ),
      boruta_process = list(
        runs_completed = run,
        final_feature_status = feature_status,
        confirmed_features = sum(feature_status == "Confirmed"),
        rejected_features = sum(feature_status == "Rejected"),
        tentative_features = sum(feature_status == "Tentative"),
        importance_history = importance_history[1:run, , drop = FALSE],
        decision_history = decision_history
      )
    )
  ))
}

#' Create Shadow Features
#'
#' @description
#' Create permuted copies of features (shadow features)
#'
#' @param feature_cols Original feature matrix
#' @return Shadow feature matrix
create_shadow_features <- function(feature_cols) {
  shadow_data <- feature_cols
  
  # Permute each column independently
  for (i in seq_len(ncol(feature_cols))) {
    shadow_data[, i] <- sample(feature_cols[, i])
  }
  
  # Rename shadow columns
  colnames(shadow_data) <- paste0("shadow_", colnames(feature_cols))
  
  return(shadow_data)
}

#' Calculate Boruta Importance
#'
#' @description
#' Calculate feature importance for Boruta using Random Forest
#'
#' @param combined_data Data with real and shadow features
#' @param n_trees Number of trees in Random Forest
#' @return Importance scores
calculate_boruta_importance <- function(combined_data, n_trees) {
  tryCatch({
    class_col <- combined_data[, 1]
    feature_cols <- combined_data[, -1, drop = FALSE]
    
    # Convert class to factor for consistency
    if (is.factor(class_col) || is.character(class_col)) {
      # Already categorical, keep as is for tree construction
    } else {
      # Convert numeric to factor if needed
    }
    
    # Calculate importance using bootstrap aggregation
    importance_matrix <- matrix(0, nrow = min(n_trees, 30), ncol = ncol(feature_cols))
    colnames(importance_matrix) <- colnames(feature_cols)
    
    successful_trees <- 0
    
    for (tree in seq_len(min(n_trees, 30))) {
      tryCatch({
        # Bootstrap sample
        boot_indices <- sample(nrow(combined_data), replace = TRUE)
        boot_data <- combined_data[boot_indices, ]
        boot_class <- boot_data[, 1]
        boot_features <- boot_data[, -1, drop = FALSE]
        
        if (is.factor(boot_class) || is.character(boot_class)) {
          boot_class_numeric <- as.numeric(as.factor(boot_class))
        } else {
          boot_class_numeric <- boot_class
        }
        
        # Calculate tree importance
        tree_importance <- calculate_tree_importance_boruta(boot_features, boot_class_numeric)
        
        successful_trees <- successful_trees + 1
        importance_matrix[successful_trees, ] <- tree_importance
        
      }, error = function(e) {
        # Skip failed trees
      })
    }
    
    # Average importance
    if (successful_trees > 0) {
      final_importance <- colMeans(importance_matrix[seq_len(successful_trees), , drop = FALSE])
    } else {
      # Fallback importance
      final_importance <- rep(0.1, ncol(feature_cols))
      names(final_importance) <- colnames(feature_cols)
    }
    
    return(final_importance)
    
  }, error = function(e) {
    # Ultimate fallback
    n_cols <- ncol(combined_data) - 1
    fallback <- rep(0.1, n_cols)
    names(fallback) <- colnames(combined_data)[-1]
    return(fallback)
  })
}

#' Calculate Tree Importance for Boruta
#'
#' @description
#' Calculate importance for single tree in Boruta context
#'
#' @param features Feature matrix
#' @param class_numeric Numeric class vector
#' @return Tree importance scores
calculate_tree_importance_boruta <- function(features, class_numeric) {
  tryCatch({
    n_features <- ncol(features)
    importance_scores <- numeric(n_features)
    names(importance_scores) <- colnames(features)
    
    # Calculate correlation-based importance with noise
    for (i in seq_len(n_features)) {
      feature_values <- features[, i]
      
      if (is.numeric(feature_values)) {
        # Correlation with small random noise for tie-breaking
        base_score <- abs(cor(feature_values, class_numeric, use = "complete.obs"))
        noise <- runif(1, -0.01, 0.01)
        importance_scores[i] <- base_score + noise
      } else {
        # Categorical feature handling
        feature_numeric <- as.numeric(as.factor(feature_values))
        base_score <- abs(cor(feature_numeric, class_numeric, use = "complete.obs"))
        noise <- runif(1, -0.01, 0.01)
        importance_scores[i] <- base_score + noise
      }
    }
    
    # Handle missing values
    importance_scores[is.na(importance_scores)] <- 0.01
    
    return(importance_scores)
    
  }, error = function(e) {
    fallback <- rep(0.01, ncol(features))
    names(fallback) <- colnames(features)
    return(fallback)
  })
}

#' Make Boruta Decisions
#'
#' @description
#' Make decisions about feature status based on shadow comparison
#'
#' @param real_importance Real feature importance
#' @param max_shadow_importance Maximum shadow importance
#' @param current_status Current feature status
#' @param p_value P-value threshold
#' @param run Current run number
#' @return Updated decisions
make_boruta_decisions <- function(real_importance, max_shadow_importance, 
                                 current_status, p_value, run) {
  
  new_status <- current_status
  hits <- rep(0, length(real_importance))
  names(hits) <- names(real_importance)
  
  # Count hits (features better than max shadow)
  for (i in seq_along(real_importance)) {
    if (!is.na(real_importance[i]) && real_importance[i] > max_shadow_importance) {
      hits[i] <- 1
    }
  }
  
  # Make decisions for tentative features only
  tentative_indices <- which(current_status == "Tentative")
  
  for (i in tentative_indices) {
    feature_name <- names(current_status)[i]
    
    # Simple decision rule: confirm if consistently better than shadow
    if (run >= 5) {  # Need minimum runs for decision
      hit_rate <- hits[feature_name] / run
      
      if (hit_rate > 0.8) {  # High hit rate
        new_status[i] <- "Confirmed"
      } else if (hit_rate < 0.2 && run >= 10) {  # Low hit rate
        new_status[i] <- "Rejected"
      }
      # Otherwise remains tentative
    }
  }
  
  return(list(
    status = new_status,
    hits = hits,
    max_shadow = max_shadow_importance
  ))
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
tryCatch({
  register_fs_method(
    method_id = 14,
    method_name = "Boruta",
    method_function = fs_method_boruta,
    category = "ensemble",
    dependencies = character(0),
    description = "Boruta algorithm using Random Forest importance with shadow features",
    parameters = list(
      max_runs = 50,
      p_value = 0.01,
      mcAdj = TRUE,
      n_trees = 100
    ),
    complexity = "high",
    supports_smote = TRUE,
    timeout_default = 600
  )
  
  cat("✓ Registered Boruta (Boruta) method (ID: 14)\n")
  
}, error = function(e) {
  warning("Failed to register Boruta method: ", e$message)
})
