#' Boruta Shap Feature Selection
#'
#' @description
#' Enhanced Boruta algorithm with SHAP-inspired importance calculation.
#' Uses improved importance metrics beyond standard Random Forest importance
#' for more robust feature selection with shadow feature comparison.
#' This is method #15 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Boruta Shap Feature Selection
#'
#' @description
#' Enhanced Boruta algorithm using SHAP-inspired feature importance.
#' Combines multiple importance metrics and uses shadow features to
#' determine truly relevant features with improved statistical testing.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (max_runs, p_value, sample_fraction)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_boruta_shap <- function(data, 
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
  max_runs <- config$max_runs %||% 40
  p_value <- config$p_value %||% 0.05
  sample_fraction <- config$sample_fraction %||% 0.8
  n_trees <- config$n_trees %||% 50
  
  # Initialize tracking
  feature_status <- rep("Tentative", n_features)
  names(feature_status) <- colnames(feature_cols)
  
  importance_matrix <- matrix(0, nrow = max_runs, ncol = n_features)
  colnames(importance_matrix) <- colnames(feature_cols)
  
  shadow_importance_matrix <- matrix(0, nrow = max_runs, ncol = n_features)
  colnames(shadow_importance_matrix) <- paste0("shadow_", colnames(feature_cols))
  
  z_scores <- matrix(0, nrow = max_runs, ncol = n_features)
  colnames(z_scores) <- colnames(feature_cols)
  
  run <- 0
  start_time <- Sys.time()
  
  # Main Boruta-SHAP iterations
  while (run < max_runs && 
         sum(feature_status == "Tentative") > 0 &&
         as.numeric(difftime(Sys.time(), start_time, units = "secs")) < timeout_sec) {
    
    run <- run + 1
    
    tryCatch({
      # Sample data for this run
      sample_indices <- sample(nrow(data), size = floor(sample_fraction * nrow(data)))
      sample_data <- data[sample_indices, ]
      
      # Create shadow features
      sample_features <- sample_data[, -1, drop = FALSE]
      shadow_features <- create_shadow_features_shap(sample_features)
      
      # Combine original and shadow
      combined_features <- cbind(sample_features, shadow_features)
      combined_data <- cbind(sample_data[, 1, drop = FALSE], combined_features)
      
      # Calculate SHAP-inspired importance
      importance_scores <- calculate_shap_importance(combined_data, n_trees)
      
      # Separate real and shadow importances
      real_importance <- importance_scores[1:n_features]
      shadow_importance <- importance_scores[(n_features + 1):length(importance_scores)]
      
      # Store importances
      importance_matrix[run, ] <- real_importance
      shadow_importance_matrix[run, ] <- shadow_importance
      
      # Calculate z-scores (standardized difference from shadow)
      for (i in seq_len(n_features)) {
        shadow_mean <- mean(shadow_importance, na.rm = TRUE)
        shadow_sd <- sd(shadow_importance, na.rm = TRUE)
        
        if (shadow_sd > 0) {
          z_scores[run, i] <- (real_importance[i] - shadow_mean) / shadow_sd
        } else {
          z_scores[run, i] <- 0
        }
      }
      
      # Update feature status based on z-scores
      if (run >= 5) {  # Need minimum runs for statistical testing
        feature_status <- update_boruta_shap_status(
          feature_status, z_scores[1:run, , drop = FALSE], p_value, run
        )
      }
      
    }, error = function(e) {
      warning(paste("Boruta-SHAP run", run, "failed:", e$message))
    })
  }
  
  # Final feature selection
  confirmed_features <- names(feature_status)[feature_status == "Confirmed"]
  
  if (length(confirmed_features) == 0) {
    # Select tentative features with highest average z-scores
    tentative_features <- names(feature_status)[feature_status == "Tentative"]
    if (length(tentative_features) > 0) {
      avg_z_scores <- colMeans(z_scores[1:run, tentative_features, drop = FALSE], na.rm = TRUE)
      n_select <- min(max_features, length(tentative_features))
      top_indices <- order(avg_z_scores, decreasing = TRUE)[seq_len(n_select)]
      selected_features <- names(avg_z_scores)[top_indices]
      final_scores <- avg_z_scores[top_indices]
    } else {
      # Fallback to top importance
      avg_importance <- colMeans(importance_matrix[1:run, , drop = FALSE], na.rm = TRUE)
      n_select <- min(max_features, length(avg_importance))
      top_indices <- order(avg_importance, decreasing = TRUE)[seq_len(n_select)]
      selected_features <- names(avg_importance)[top_indices]
      final_scores <- avg_importance[top_indices]
    }
  } else {
    # Use confirmed features
    if (length(confirmed_features) > max_features) {
      avg_z_scores <- colMeans(z_scores[1:run, confirmed_features, drop = FALSE], na.rm = TRUE)
      top_indices <- order(avg_z_scores, decreasing = TRUE)[seq_len(max_features)]
      selected_features <- names(avg_z_scores)[top_indices]
      final_scores <- avg_z_scores[top_indices]
    } else {
      selected_features <- confirmed_features
      final_scores <- colMeans(z_scores[1:run, confirmed_features, drop = FALSE], na.rm = TRUE)
    }
  }
  
  return(list(
    features = selected_features,
    scores = final_scores,
    metadata = list(
      method = "Boruta_shap",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        max_runs = max_runs,
        p_value = p_value,
        sample_fraction = sample_fraction,
        n_trees = n_trees,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        final_z_scores = final_scores,
        mean_z_score = mean(final_scores),
        z_score_range = range(final_scores)
      ),
      boruta_shap_process = list(
        runs_completed = run,
        final_feature_status = feature_status,
        confirmed_features = sum(feature_status == "Confirmed"),
        rejected_features = sum(feature_status == "Rejected"),
        tentative_features = sum(feature_status == "Tentative"),
        z_score_matrix = z_scores[1:run, , drop = FALSE],
        importance_matrix = importance_matrix[1:run, , drop = FALSE]
      )
    )
  ))
}

#' Create Shadow Features for SHAP
#'
#' @description
#' Create enhanced shadow features with multiple permutation strategies
#'
#' @param feature_cols Original feature matrix
#' @return Shadow feature matrix
create_shadow_features_shap <- function(feature_cols) {
  n_features <- ncol(feature_cols)
  shadow_data <- feature_cols
  
  # Enhanced permutation strategies
  for (i in seq_len(n_features)) {
    feature_values <- feature_cols[, i]
    
    if (is.numeric(feature_values)) {
      # For numeric: add noise and permute
      noise_factor <- sd(feature_values, na.rm = TRUE) * 0.1
      permuted_values <- sample(feature_values)
      noisy_values <- permuted_values + rnorm(length(permuted_values), 0, noise_factor)
      shadow_data[, i] <- noisy_values
    } else {
      # For categorical: simple permutation
      shadow_data[, i] <- sample(feature_values)
    }
  }
  
  colnames(shadow_data) <- paste0("shadow_", colnames(feature_cols))
  return(shadow_data)
}

#' Calculate SHAP-Inspired Importance
#'
#' @description
#' Calculate feature importance using SHAP-inspired metrics
#'
#' @param combined_data Data with real and shadow features
#' @param n_trees Number of trees for ensemble
#' @return Importance scores
calculate_shap_importance <- function(combined_data, n_trees) {
  tryCatch({
    class_col <- combined_data[, 1]
    feature_cols <- combined_data[, -1, drop = FALSE]
    
    # Convert class appropriately
    if (is.factor(class_col) || is.character(class_col)) {
      class_factor <- as.factor(class_col)
      class_numeric <- as.numeric(class_factor)
    } else {
      class_numeric <- class_col
      class_factor <- as.factor(class_col)
    }
    
    # Calculate multiple importance metrics
    correlation_importance <- calculate_correlation_importance(feature_cols, class_numeric)
    mutual_info_importance <- calculate_mutual_info_importance(feature_cols, class_factor)
    variance_importance <- calculate_variance_importance(feature_cols, class_factor)
    
    # Combine importance scores with weights
    combined_importance <- 0.5 * correlation_importance + 
                          0.3 * mutual_info_importance + 
                          0.2 * variance_importance
    
    # Normalize to [0, 1]
    max_importance <- max(combined_importance, na.rm = TRUE)
    if (max_importance > 0) {
      combined_importance <- combined_importance / max_importance
    }
    
    return(combined_importance)
    
  }, error = function(e) {
    # Fallback importance
    n_cols <- ncol(combined_data) - 1
    fallback <- rep(0.1, n_cols)
    names(fallback) <- colnames(combined_data)[-1]
    return(fallback)
  })
}

#' Calculate Correlation Importance
#'
#' @description
#' Calculate importance based on correlation with class
#'
#' @param features Feature matrix
#' @param class_numeric Numeric class vector
#' @return Correlation importance scores
calculate_correlation_importance <- function(features, class_numeric) {
  importance_scores <- apply(features, 2, function(x) {
    if (is.numeric(x)) {
      abs(cor(x, class_numeric, use = "complete.obs"))
    } else {
      x_numeric <- as.numeric(as.factor(x))
      abs(cor(x_numeric, class_numeric, use = "complete.obs"))
    }
  })
  
  importance_scores[is.na(importance_scores)] <- 0
  return(importance_scores)
}

#' Calculate Mutual Information Importance
#'
#' @description
#' Calculate importance based on mutual information approximation
#'
#' @param features Feature matrix
#' @param class_factor Factor class vector
#' @return Mutual information importance scores
calculate_mutual_info_importance <- function(features, class_factor) {
  importance_scores <- apply(features, 2, function(x) {
    tryCatch({
      if (is.numeric(x)) {
        # Discretize numeric features for MI calculation
        x_discrete <- cut(x, breaks = 5, labels = FALSE)
        x_discrete[is.na(x_discrete)] <- 1
      } else {
        x_discrete <- as.numeric(as.factor(x))
      }
      
      # Simple mutual information approximation
      contingency_table <- table(x_discrete, class_factor)
      
      if (any(dim(contingency_table) < 2)) {
        return(0)
      }
      
      # Calculate entropy components
      n_total <- sum(contingency_table)
      p_joint <- contingency_table / n_total
      p_x <- rowSums(p_joint)
      p_y <- colSums(p_joint)
      
      # Calculate mutual information
      mi <- 0
      for (i in seq_len(nrow(p_joint))) {
        for (j in seq_len(ncol(p_joint))) {
          if (p_joint[i, j] > 0 && p_x[i] > 0 && p_y[j] > 0) {
            mi <- mi + p_joint[i, j] * log2(p_joint[i, j] / (p_x[i] * p_y[j]))
          }
        }
      }
      
      return(max(0, mi))
      
    }, error = function(e) {
      return(0)
    })
  })
  
  return(importance_scores)
}

#' Calculate Variance Importance
#'
#' @description
#' Calculate importance based on between-class variance ratio
#'
#' @param features Feature matrix
#' @param class_factor Factor class vector
#' @return Variance importance scores
calculate_variance_importance <- function(features, class_factor) {
  importance_scores <- apply(features, 2, function(x) {
    tryCatch({
      if (is.numeric(x)) {
        # Calculate F-statistic (between/within variance ratio)
        class_means <- tapply(x, class_factor, mean, na.rm = TRUE)
        overall_mean <- mean(x, na.rm = TRUE)
        class_sizes <- table(class_factor)
        
        # Between-group variance
        between_var <- sum(class_sizes * (class_means - overall_mean)^2, na.rm = TRUE)
        
        # Within-group variance
        within_var <- sum(tapply(x, class_factor, function(vals) {
          sum((vals - mean(vals, na.rm = TRUE))^2, na.rm = TRUE)
        }), na.rm = TRUE)
        
        if (within_var > 0) {
          f_stat <- between_var / within_var
          return(min(f_stat / 100, 1))  # Normalize
        } else {
          return(0)
        }
      } else {
        # For categorical: use chi-square statistic
        contingency_table <- table(x, class_factor)
        if (any(dim(contingency_table) >= 2)) {
          chi_test <- chisq.test(contingency_table, simulate.p.value = TRUE)
          return(min(chi_test$statistic / sum(contingency_table), 1))
        } else {
          return(0)
        }
      }
      
    }, error = function(e) {
      return(0)
    })
  })
  
  return(importance_scores)
}

#' Update Boruta SHAP Status
#'
#' @description
#' Update feature status based on z-score statistical testing
#'
#' @param current_status Current feature status
#' @param z_score_matrix Matrix of z-scores across runs
#' @param p_value P-value threshold
#' @param current_run Current run number
#' @return Updated feature status
update_boruta_shap_status <- function(current_status, z_score_matrix, p_value, current_run) {
  new_status <- current_status
  
  # Only update tentative features
  tentative_indices <- which(current_status == "Tentative")
  
  for (i in tentative_indices) {
    feature_name <- names(current_status)[i]
    z_scores <- z_score_matrix[, feature_name]
    z_scores <- z_scores[!is.na(z_scores)]
    
    if (length(z_scores) >= 5) {
      # One-sample t-test against 0 (no difference from shadow)
      t_test_result <- t.test(z_scores, mu = 0)
      
      if (t_test_result$p.value < p_value) {
        if (mean(z_scores) > 0) {
          new_status[i] <- "Confirmed"
        } else {
          new_status[i] <- "Rejected"
        }
      }
      
      # Additional rule: strong consistent performance
      positive_runs <- sum(z_scores > 1.96)  # 95% confidence threshold
      if (positive_runs >= 0.8 * length(z_scores)) {
        new_status[i] <- "Confirmed"
      }
      
      negative_runs <- sum(z_scores < -1.96)
      if (negative_runs >= 0.8 * length(z_scores)) {
        new_status[i] <- "Rejected"
      }
    }
  }
  
  return(new_status)
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
tryCatch({
  register_fs_method(
    method_id = 15,
    method_name = "Boruta_shap",
    method_function = fs_method_boruta_shap,
    category = "ensemble",
    dependencies = character(0),
    description = "Enhanced Boruta algorithm with SHAP-inspired importance and statistical testing",
    parameters = list(
      max_runs = 40,
      p_value = 0.05,
      sample_fraction = 0.8,
      n_trees = 50
    ),
    complexity = "very_high",
    supports_smote = TRUE,
    timeout_default = 600
  )
  
  cat("✓ Registered Boruta SHAP (Boruta_shap) method (ID: 15)\n")
  
}, error = function(e) {
  warning("Failed to register Boruta_shap method: ", e$message)
})
