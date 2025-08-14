#' Optimized Helper Functions for OmicSelector
#'
#' @description
#' Helper functions specifically for the optimized OmicSelector implementation
#' for miRNA Case-Control analysis.

#' Setup Logging for OmicSelector
#'
#' @param output_dir Character. Output directory for logs
#' @param prefix Character. Log file prefix
#' @keywords internal
setup_omicselector_logging <- function(output_dir, prefix = "omicselector") {
  if (!requireNamespace("logger", quietly = TRUE)) {
    # Fallback logging without logger package
    return(invisible(NULL))
  }
  
  log_file <- file.path(output_dir, paste0(prefix, "_", Sys.Date(), ".log"))
  
  # Setup logger with both console and file output
  logger::log_appender(logger::appender_tee(log_file))
  logger::log_layout(logger::layout_glue_colors)
  logger::log_threshold(logger::INFO)
  
  return(log_file)
}

#' Log Info Message
#' @param message Character message to log
#' @keywords internal
log_info <- function(message) {
  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info(message)
  } else {
    cat("[INFO]", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message, "\n")
  }
}

#' Log Warning Message
#' @param message Character message to log
#' @keywords internal
log_warn <- function(message) {
  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_warn(message)
  } else {
    cat("[WARN]", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message, "\n")
  }
}

#' Log Error Message
#' @param message Character message to log
#' @keywords internal
log_error <- function(message) {
  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_error(message)
  } else {
    cat("[ERROR]", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message, "\n")
  }
}

#' Create Progress Bar
#' @param total Integer total number of items
#' @param title Character title for progress bar
#' @return Progress bar object
#' @keywords internal
create_progress_bar <- function(total, title = "Progress") {
  if (requireNamespace("cli", quietly = TRUE)) {
    # cli_progress_bar returns an ID, we need to wrap it
    pb_id <- cli::cli_progress_bar(title, total = total)
    list(
      tick = function() cli::cli_progress_update(id = pb_id),
      terminate = function() cli::cli_progress_done(id = pb_id),
      id = pb_id
    )
  } else {
    # Fallback progress bar
    count <- 0
    list(
      tick = function() {
        count <<- count + 1
        cat(".")
        if (count %% 50 == 0) cat("\n")
      },
      terminate = function() cat("\n")
    )
  }
}

#' Safe Timeout Function
#'
#' @param expr Expression to execute
#' @param timeout Timeout in seconds
#' @return Result of expression or timeout error
#' @keywords internal
withTimeout <- function(expr, timeout) {
  if (requireNamespace("R.utils", quietly = TRUE)) {
    R.utils::withTimeout(expr, timeout = timeout)
  } else {
    # Fallback without timeout (just execute)
    eval(expr, envir = parent.frame())
  }
}

#' Null Coalescing Operator
#' @param x Left operand
#' @param y Right operand (default if x is NULL)
#' @return x if not NULL, otherwise y
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Fast Implementation of Correlation-based Feature Selection
#'
#' @param data Data frame with features and class
#' @param class_column Character name of class column
#' @param max_features Integer maximum features to return
#' @return Character vector of selected features
#' @keywords internal
fast_cfs <- function(data, class_column, max_features = 20) {
  
  # Get feature columns (excluding class)
  feature_cols <- setdiff(colnames(data), class_column)
  class_values <- data[[class_column]]
  
  if (length(feature_cols) == 0) {
    return(character(0))
  }
  
  # Convert class to numeric for correlation
  class_numeric <- as.numeric(as.factor(class_values))
  
  # Calculate feature-class correlations
  feature_class_cors <- sapply(feature_cols, function(col) {
    feature_values <- data[[col]]
    if (all(is.na(feature_values)) || var(feature_values, na.rm = TRUE) == 0) {
      return(0)
    }
    abs(cor(feature_values, class_numeric, use = "complete.obs"))
  })
  
  # Select top features by correlation with class
  top_features <- names(sort(feature_class_cors, decreasing = TRUE))[seq_len(min(max_features, length(feature_class_cors)))]
  
  return(top_features)
}

#' Fast Random Forest Feature Importance
#'
#' @param data Data frame with features and class
#' @param class_column Character name of class column
#' @param max_features Integer maximum features to return
#' @return Character vector of selected features
#' @keywords internal
fast_rf_importance <- function(data, class_column, max_features = 20) {
  
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    # Fallback to correlation-based selection
    return(fast_cfs(data, class_column, max_features))
  }
  
  tryCatch({
    # Prepare data
    formula_str <- paste(class_column, "~ .")
    rf_formula <- as.formula(formula_str)
    
    # Train random forest with reduced trees for speed
    rf_model <- randomForest::randomForest(
      formula = rf_formula,
      data = data,
      ntree = 100,  # Reduced for speed
      importance = TRUE,
      na.action = na.omit
    )
    
    # Get feature importance
    importance_scores <- randomForest::importance(rf_model)[, "MeanDecreaseGini"]
    
    # Select top features
    top_features <- names(sort(importance_scores, decreasing = TRUE))[seq_len(min(max_features, length(importance_scores)))]
    
    return(top_features)
    
  }, error = function(e) {
    # Fallback on error
    return(fast_cfs(data, class_column, max_features))
  })
}

#' Optimized t-test for Multiple Features
#'
#' @param feature_matrix Matrix of features (samples x features)
#' @param class_labels Factor of class labels
#' @param case_label Character positive case label
#' @param control_label Character negative control label
#' @return Data frame with statistical results
#' @keywords internal
fast_ttest <- function(feature_matrix, class_labels, case_label, control_label) {
  
  case_indices <- which(class_labels == case_label)
  control_indices <- which(class_labels == control_label)
  
  # Pre-allocate results
  n_features <- ncol(feature_matrix)
  feature_names <- colnames(feature_matrix)
  
  results <- data.frame(
    feature = feature_names,
    mean_case = numeric(n_features),
    mean_control = numeric(n_features), 
    log2fc = numeric(n_features),
    p_value = numeric(n_features),
    stringsAsFactors = FALSE
  )
  
  # Vectorized calculations where possible
  case_data <- feature_matrix[case_indices, , drop = FALSE]
  control_data <- feature_matrix[control_indices, , drop = FALSE]
  
  # Calculate means
  results$mean_case <- colMeans(case_data, na.rm = TRUE)
  results$mean_control <- colMeans(control_data, na.rm = TRUE)
  
  # Calculate log2 fold change
  pseudocount <- 1e-6
  results$log2fc <- log2((results$mean_case + pseudocount) / (results$mean_control + pseudocount))
  
  # Calculate p-values using vectorized t-test
  for (i in seq_len(n_features)) {
    case_vals <- case_data[, i]
    control_vals <- control_data[, i]
    
    # Remove NA values
    case_vals <- case_vals[!is.na(case_vals)]
    control_vals <- control_vals[!is.na(control_vals)]
    
    if (length(case_vals) > 1 && length(control_vals) > 1) {
      tryCatch({
        t_result <- t.test(case_vals, control_vals)
        results$p_value[i] <- t_result$p.value
      }, error = function(e) {
        results$p_value[i] <- 1.0
      })
    } else {
      results$p_value[i] <- 1.0
    }
  }
  
  return(results)
}

#' Create Feature Formula for Modeling
#'
#' @param features Character vector of feature names
#' @param class_column Character name of class column (default: "Class")
#' @return Formula object
#' @keywords internal
create_feature_formula <- function(features, class_column = "Class") {
  if (length(features) == 0) {
    return(as.formula(paste(class_column, "~ 1")))  # Intercept only
  }
  
  # Escape feature names that might have special characters
  escaped_features <- paste0("`", features, "`")
  feature_string <- paste(escaped_features, collapse = " + ")
  formula_string <- paste(class_column, "~", feature_string)
  
  return(as.formula(formula_string))
}

#' Validate miRNA Feature Names
#'
#' @param features Character vector of feature names
#' @param prefix Character expected prefix (default: "hsa_")
#' @return Logical vector indicating valid features
#' @keywords internal
validate_mirna_features <- function(features, prefix = "hsa_") {
  
  # Check if features start with expected prefix
  has_prefix <- startsWith(features, prefix)
  
  # Check for valid miRNA naming pattern (basic validation)
  # Pattern: hsa_miR_xxx or hsa_let_xxx
  valid_pattern <- grepl(paste0("^", prefix, "(miR|let)"), features)
  
  return(has_prefix & valid_pattern)
}

#' Filter Data for miRNA Features Only
#'
#' @param data Data frame containing mixed columns
#' @param class_column Character name of class column
#' @param feature_prefix Character prefix for miRNA features
#' @return Data frame with class column and miRNA features only
#' @keywords internal
filter_mirna_data <- function(data, class_column = "Class", feature_prefix = "hsa_") {
  
  # Get miRNA feature columns
  mirna_cols <- colnames(data)[startsWith(colnames(data), feature_prefix)]
  
  # Validate miRNA features
  valid_mirna <- validate_mirna_features(mirna_cols, feature_prefix)
  mirna_cols <- mirna_cols[valid_mirna]
  
  if (length(mirna_cols) == 0) {
    stop("No valid miRNA features found with prefix: ", feature_prefix)
  }
  
  # Return filtered data
  selected_cols <- c(class_column, mirna_cols)
  return(data[, selected_cols, drop = FALSE])
}

#' Calculate Feature Selection Metrics
#'
#' @param selected_features List of feature sets from different methods
#' @param de_results Data frame of differential expression results
#' @return Data frame with selection metrics
#' @keywords internal
calculate_selection_metrics <- function(selected_features, de_results) {
  
  method_names <- names(selected_features)
  n_methods <- length(method_names)
  
  metrics <- data.frame(
    method = method_names,
    n_features = integer(n_methods),
    n_significant = integer(n_methods),
    prop_significant = numeric(n_methods),
    mean_pvalue = numeric(n_methods),
    mean_abs_log2fc = numeric(n_methods),
    stringsAsFactors = FALSE
  )
  
  significant_features <- de_results$feature[de_results$significant]
  
  for (i in seq_len(n_methods)) {
    features <- selected_features[[i]]
    
    if (length(features) > 0) {
      # Match with DE results
      feature_stats <- de_results[de_results$feature %in% features, ]
      
      metrics$n_features[i] <- length(features)
      metrics$n_significant[i] <- sum(features %in% significant_features)
      metrics$prop_significant[i] <- metrics$n_significant[i] / metrics$n_features[i]
      
      if (nrow(feature_stats) > 0) {
        metrics$mean_pvalue[i] <- mean(feature_stats$p_value_adj, na.rm = TRUE)
        metrics$mean_abs_log2fc[i] <- mean(abs(feature_stats$log2fc), na.rm = TRUE)
      }
    }
  }
  
  return(metrics)
}

#' Export Results to Multiple Formats
#'
#' @param results OmicSelector results object
#' @param output_dir Character output directory
#' @param prefix Character file prefix
#' @keywords internal
export_omicselector_results <- function(results, output_dir, prefix = "omicselector") {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Export selected features to CSV
  features_df <- data.frame(
    method = rep(names(results$selected_features), 
                sapply(results$selected_features, length)),
    feature = unlist(results$selected_features),
    stringsAsFactors = FALSE
  )
  
  if (nrow(features_df) > 0) {
    features_file <- file.path(output_dir, paste0(prefix, "_selected_features_", timestamp, ".csv"))
    write.csv(features_df, features_file, row.names = FALSE)
  }
  
  # Export differential expression results
  de_file <- file.path(output_dir, paste0(prefix, "_differential_expression_", timestamp, ".csv"))
  write.csv(results$differential_expression, de_file, row.names = FALSE)
  
  # Export summary statistics
  if (!is.null(results$summary)) {
    summary_file <- file.path(output_dir, paste0(prefix, "_summary_", timestamp, ".txt"))
    sink(summary_file)
    print(results)
    sink()
  }
  
  # Export complete results as RDS
  rds_file <- file.path(output_dir, paste0(prefix, "_complete_results_", timestamp, ".rds"))
  saveRDS(results, rds_file)
  
  return(list(
    features_csv = if(exists("features_file")) features_file else NULL,
    de_csv = de_file,
    summary_txt = if(exists("summary_file")) summary_file else NULL,
    complete_rds = rds_file
  ))
}
