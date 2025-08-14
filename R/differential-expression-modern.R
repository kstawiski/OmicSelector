#' Modern Differential Expression Analysis
#'
#' @description
#' Modernized differential expression analysis using t-test with multiple testing correction.
#' This function replaces the excessive suppressMessages() usage and provides better
#' error handling and logging.
#'
#' @param features Matrix or data frame with features as columns and samples as rows
#' @param classes Factor or character vector with sample classes ("Case", "Control")
#' @param mode Character string: "logtpm", "deltact", or "auto"
#' @param correction_method Character string: "BH" (default), "bonferroni", "holm", etc.
#' @param alpha Numeric: significance threshold (default 0.05)
#' @param log_fold_change Logical: whether to calculate log fold change (default TRUE)
#'
#' @return Data frame with differential expression results
#'
#' @examples
#' \dontrun{
#' # Create example data
#' set.seed(123)
#' features <- matrix(rnorm(1000), nrow = 50, ncol = 20)
#' colnames(features) <- paste0("feature_", 1:20)
#' classes <- factor(rep(c("Case", "Control"), each = 25))
#' 
#' # Run differential expression
#' results <- differential_expression_modern(features, classes)
#' }
#'
#' @export
differential_expression_modern <- function(features, 
                                         classes, 
                                         mode = "auto",
                                         correction_method = "BH",
                                         alpha = 0.05,
                                         log_fold_change = TRUE) {
  
  log_info("Starting differential expression analysis")
  
  # Input validation
  if (!is.matrix(features) && !is.data.frame(features)) {
    stop("Features must be a matrix or data frame", call. = FALSE)
  }
  
  if (length(classes) != nrow(features)) {
    stop("Length of classes must equal number of rows in features", call. = FALSE)
  }
  
  # Handle auto mode
  if (mode == "auto") {
    if (file.exists("var_type.txt")) {
      mode <- readLines("var_type.txt", warn = FALSE, n = 1)
      log_debug("Auto mode detected type: {mode}")
    } else {
      mode <- "logtpm"
      log_debug("Auto mode defaulting to: {mode}")
    }
  }
  
  # Ensure classes are factor with correct levels
  classes <- factor(classes)
  unique_classes <- levels(classes)
  
  if (length(unique_classes) != 2) {
    stop("Exactly two classes required, found: ", 
         paste(unique_classes, collapse = ", "), call. = FALSE)
  }
  
  log_info("Analyzing {ncol(features)} features across {nrow(features)} samples")
  log_info("Classes: {paste(unique_classes, collapse = ' vs ')}")
  log_info("Class distribution: {table(classes)}")
  
  # Convert to data frame if matrix
  if (is.matrix(features)) {
    features <- as.data.frame(features)
  }
  
  # Initialize results
  results <- data.frame(
    feature = colnames(features),
    stringsAsFactors = FALSE
  )
  
  # Perform t-tests with error handling
  log_info("Performing t-tests...")
  pb <- create_progress_bar(ncol(features), "Running t-tests")
  
  p_values <- numeric(ncol(features))
  fold_changes <- numeric(ncol(features))
  mean_case <- numeric(ncol(features))
  mean_control <- numeric(ncol(features))
  
  for (i in seq_len(ncol(features))) {
    
    tryCatch({
      feature_values <- features[, i]
      
      # Remove missing values
      valid_indices <- !is.na(feature_values)
      if (sum(valid_indices) < 4) {
        log_warn("Feature {colnames(features)[i]} has too few valid values, skipping")
        p_values[i] <- 1.0
        next
      }
      
      feature_clean <- feature_values[valid_indices]
      classes_clean <- classes[valid_indices]
      
      # Perform t-test
      test_result <- t.test(feature_clean ~ classes_clean)
      p_values[i] <- test_result$p.value
      
      # Calculate means and fold changes
      case_values <- feature_clean[classes_clean == unique_classes[1]]
      control_values <- feature_clean[classes_clean == unique_classes[2]]
      
      mean_case[i] <- mean(case_values, na.rm = TRUE)
      mean_control[i] <- mean(control_values, na.rm = TRUE)
      
      if (log_fold_change && mode %in% c("logtpm", "auto")) {
        # For log-transformed data, difference is log fold change
        fold_changes[i] <- mean_case[i] - mean_control[i]
      } else if (mode == "deltact") {
        # For deltaCt data, reverse the calculation
        fold_changes[i] <- mean_control[i] - mean_case[i]
      } else {
        # For linear data, calculate actual fold change
        fold_changes[i] <- log2((mean_case[i] + 1e-6) / (mean_control[i] + 1e-6))
      }
      
    }, error = function(e) {
      log_warn("Error in t-test for feature {colnames(features)[i]}: {e$message}")
      p_values[i] <- 1.0
      fold_changes[i] <- 0
    })
    
    update_progress_bar(pb, i)
  }
  
  finish_progress_bar(pb)
  
  # Multiple testing correction
  log_info("Applying multiple testing correction: {correction_method}")
  p_adjusted <- p.adjust(p_values, method = correction_method)
  
  # Compile results
  results$p_value <- p_values
  results$p_adjusted <- p_adjusted
  results$log_fold_change <- fold_changes
  results$mean_case <- mean_case
  results$mean_control <- mean_control
  results$significant <- p_adjusted < alpha
  
  # Add effect size (if possible)
  if (requireNamespace("effsize", quietly = TRUE)) {
    log_debug("Calculating effect sizes...")
    effect_sizes <- numeric(ncol(features))
    
    for (i in seq_len(ncol(features))) {
      tryCatch({
        feature_values <- features[, i]
        valid_indices <- !is.na(feature_values)
        
        if (sum(valid_indices) >= 4) {
          cohen_d <- effsize::cohen.d(
            feature_values[valid_indices], 
            classes[valid_indices]
          )
          effect_sizes[i] <- cohen_d$estimate
        }
      }, error = function(e) {
        effect_sizes[i] <- NA
      })
    }
    
    results$cohens_d <- effect_sizes
  }
  
  # Sort by adjusted p-value
  results <- results[order(results$p_adjusted), ]
  rownames(results) <- NULL
  
  # Summary statistics
  n_significant <- sum(results$significant, na.rm = TRUE)
  log_info("Analysis complete: {n_significant}/{nrow(results)} features significant (p < {alpha})")
  
  if (n_significant > 0) {
    log_info("Top significant features:")
    top_features <- head(results[results$significant, ], 5)
    for (i in seq_len(nrow(top_features))) {
      log_info("  {top_features$feature[i]}: p = {signif(top_features$p_adjusted[i], 3)}, FC = {signif(top_features$log_fold_change[i], 3)}")
    }
  }
  
  return(results)
}

#' Wrapper for Backward Compatibility
#'
#' @description
#' Wrapper function that maintains compatibility with the original function
#' while using the modernized implementation internally.
#'
#' @param ttpm_features Legacy parameter name for features
#' @param classes Sample classes
#' @param mode Analysis mode
#'
#' @return Data frame with results (same format as original)
#' @export
OmicSelector_differential_expression_ttest_modern <- function(ttpm_features, classes, mode = "auto") {
  
  log_warn("Using legacy function name. Consider switching to differential_expression_modern()")
  
  # Convert legacy format to modern format
  results <- differential_expression_modern(
    features = ttpm_features,
    classes = classes,
    mode = mode
  )
  
  # Convert back to legacy format for compatibility
  legacy_results <- data.frame(
    miR = results$feature,
    p.value = results$p_value,
    p.adjust.BH = results$p_adjusted,
    log2FoldChange = results$log_fold_change,
    mean.Case = results$mean_case,
    mean.Control = results$mean_control,
    stringsAsFactors = FALSE
  )
  
  return(legacy_results)
}
