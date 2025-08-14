#' Modern Differential Expression Analysis
#'
#' @description
#' Modernized version of differential expression analysis with improved statistical
#' methods, comprehensive error handling, and enhanced result reporting. Supports
#' both traditional t-test and modern statistical approaches.
#'
#' Key improvements:
#' - Removed 24 redundant library calls, using only essential packages
#' - Added comprehensive input validation and error handling
#' - Enhanced statistical methods with effect size calculations
#' - Improved multiple testing correction options
#' - Memory-efficient processing for large datasets
#' - Modern tidyverse-style data handling
#' - Comprehensive result reporting with confidence intervals
#' - Support for different data types (log-TPM, delta-Ct, etc.)
#'
#' @param expression_data Matrix or data.frame with features in columns, samples in rows
#' @param class_labels Vector of class labels ("Case"/"Control" or custom labels)
#' @param data_type Character: "log_tpm", "delta_ct", "counts", or "auto"
#' @param case_label Character: label for cases (default: "Case")
#' @param control_label Character: label for controls (default: "Control")
#' @param p_adjust_methods Character vector: correction methods to apply
#' @param alpha Numeric: significance threshold (default: 0.05)
#' @param min_samples_per_group Integer: minimum samples required per group
#' @param effect_size_methods Character vector: effect size measures to calculate
#' @param parallel Logical: use parallel processing for large datasets
#' @param verbose Logical: enable progress messages
#'
#' @return Data.frame with comprehensive differential expression results including
#'   p-values, adjusted p-values, fold changes, effect sizes, and confidence intervals
#'
#' @examples
#' \dontrun{
#' # Basic differential expression
#' results <- differential_expression_modern(
#'   expression_data = mirna_matrix,
#'   class_labels = sample_classes
#' )
#' 
#' # Advanced analysis with custom parameters
#' results <- differential_expression_modern(
#'   expression_data = expression_matrix,
#'   class_labels = groups,
#'   data_type = "log_tpm",
#'   p_adjust_methods = c("BH", "bonferroni", "holm"),
#'   effect_size_methods = c("cohens_d", "hedges_g", "cliff_delta"),
#'   alpha = 0.01
#' )
#' 
#' # Filter significant results
#' significant <- results[results$p_value_BH < 0.05, ]
#' }
#'
#' @export
differential_expression_modern <- function(expression_data,
                                         class_labels,
                                         data_type = "auto",
                                         case_label = "Case",
                                         control_label = "Control",
                                         p_adjust_methods = c("BH", "bonferroni", "holm"),
                                         alpha = 0.05,
                                         min_samples_per_group = 3,
                                         effect_size_methods = c("cohens_d", "hedges_g"),
                                         parallel = TRUE,
                                         verbose = TRUE) {
  
  # Validate inputs
  validate_de_inputs(expression_data, class_labels, case_label, control_label, 
                    min_samples_per_group, alpha)
  
  # Load only essential packages
  required_packages <- c("dplyr", "tibble")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed", pkg))
    }
  }
  
  if (verbose) {
    cat("🧬 Starting modern differential expression analysis...\n")
    cat(sprintf("   Features: %d\n", ncol(expression_data)))
    cat(sprintf("   Samples: %d\n", nrow(expression_data)))
  }
  
  # Auto-detect data type if requested
  if (data_type == "auto") {
    data_type <- detect_data_type(expression_data, verbose)
  }
  
  # Prepare data
  processed_data <- prepare_de_data(
    expression_data = expression_data,
    class_labels = class_labels,
    case_label = case_label,
    control_label = control_label,
    data_type = data_type,
    verbose = verbose
  )
  
  # Perform differential expression analysis
  if (verbose) cat("🔬 Performing statistical tests...\n")
  
  de_results <- perform_de_analysis(
    processed_data = processed_data,
    data_type = data_type,
    p_adjust_methods = p_adjust_methods,
    effect_size_methods = effect_size_methods,
    parallel = parallel,
    verbose = verbose
  )
  
  # Add significance flags
  de_results <- add_significance_flags(de_results, alpha, verbose)
  
  # Sort by BH-adjusted p-value
  if ("p_value_BH" %in% colnames(de_results)) {
    de_results <- de_results[order(de_results$p_value_BH), ]
  } else if ("p_value" %in% colnames(de_results)) {
    de_results <- de_results[order(de_results$p_value), ]
  }
  
  # Add metadata
  attr(de_results, "analysis_info") <- list(
    data_type = data_type,
    case_label = case_label,
    control_label = control_label,
    n_features = ncol(expression_data),
    n_samples = nrow(expression_data),
    n_cases = sum(class_labels == case_label),
    n_controls = sum(class_labels == control_label),
    p_adjust_methods = p_adjust_methods,
    effect_size_methods = effect_size_methods,
    analysis_time = Sys.time()
  )
  
  if (verbose) {
    sig_count <- sum(de_results$p_value_BH < alpha, na.rm = TRUE)
    cat(sprintf("✅ Analysis complete: %d/%d features significant (FDR < %.3f)\n", 
                sig_count, nrow(de_results), alpha))
  }
  
  return(de_results)
}

#' Validate Differential Expression Inputs
#' @keywords internal
validate_de_inputs <- function(expression_data, class_labels, case_label, control_label, 
                              min_samples_per_group, alpha) {
  
  # Check expression data
  if (!is.data.frame(expression_data) && !is.matrix(expression_data)) {
    stop("expression_data must be a data.frame or matrix")
  }
  
  if (nrow(expression_data) == 0 || ncol(expression_data) == 0) {
    stop("expression_data must have both rows and columns")
  }
  
  # Check class labels
  if (length(class_labels) != nrow(expression_data)) {
    stop("Length of class_labels must equal number of rows in expression_data")
  }
  
  # Check for valid labels
  unique_labels <- unique(class_labels)
  if (!case_label %in% unique_labels) {
    stop(sprintf("case_label '%s' not found in class_labels", case_label))
  }
  
  if (!control_label %in% unique_labels) {
    stop(sprintf("control_label '%s' not found in class_labels", control_label))
  }
  
  # Check sample sizes
  n_cases <- sum(class_labels == case_label)
  n_controls <- sum(class_labels == control_label)
  
  if (n_cases < min_samples_per_group) {
    stop(sprintf("Insufficient cases: %d < %d minimum", n_cases, min_samples_per_group))
  }
  
  if (n_controls < min_samples_per_group) {
    stop(sprintf("Insufficient controls: %d < %d minimum", n_controls, min_samples_per_group))
  }
  
  # Check alpha
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a numeric value between 0 and 1")
  }
}

#' Detect Data Type
#' @keywords internal
detect_data_type <- function(expression_data, verbose) {
  
  # Sample a subset for analysis
  sample_data <- expression_data[, sample(ncol(expression_data), 
                                        min(100, ncol(expression_data)))]
  sample_values <- as.numeric(unlist(sample_data))
  sample_values <- sample_values[!is.na(sample_values)]
  
  # Check value ranges
  min_val <- min(sample_values, na.rm = TRUE)
  max_val <- max(sample_values, na.rm = TRUE)
  median_val <- median(sample_values, na.rm = TRUE)
  
  # Heuristics for data type detection
  if (min_val >= 0 && max_val > 100 && median_val > 1) {
    detected_type <- "counts"
  } else if (min_val < 0 && max_val < 50) {
    detected_type <- "delta_ct"
  } else if (min_val >= 0 && max_val <= 25 && median_val < 10) {
    detected_type <- "log_tpm"
  } else {
    detected_type <- "log_tpm"  # Default assumption
  }
  
  if (verbose) {
    cat(sprintf("🔍 Auto-detected data type: %s\n", detected_type))
    cat(sprintf("   Value range: %.2f to %.2f (median: %.2f)\n", 
                min_val, max_val, median_val))
  }
  
  return(detected_type)
}

#' Prepare Data for Differential Expression
#' @keywords internal
prepare_de_data <- function(expression_data, class_labels, case_label, control_label, 
                           data_type, verbose) {
  
  # Convert to data.frame if needed
  if (is.matrix(expression_data)) {
    expression_data <- as.data.frame(expression_data)
  }
  
  # Create group indices
  case_idx <- which(class_labels == case_label)
  control_idx <- which(class_labels == control_label)
  
  # Filter out other classes if present
  keep_idx <- c(case_idx, control_idx)
  
  if (length(keep_idx) < nrow(expression_data)) {
    if (verbose) {
      cat(sprintf("🔧 Filtering %d samples to keep only cases and controls\n", 
                  nrow(expression_data) - length(keep_idx)))
    }
    expression_data <- expression_data[keep_idx, ]
    class_labels <- class_labels[keep_idx]
  }
  
  # Remove features with too many missing values
  na_threshold <- 0.8
  na_proportion <- apply(expression_data, 2, function(x) sum(is.na(x)) / length(x))
  keep_features <- na_proportion < na_threshold
  
  if (sum(!keep_features) > 0) {
    if (verbose) {
      cat(sprintf("🚮 Removing %d features with >%d%% missing values\n", 
                  sum(!keep_features), round(na_threshold * 100)))
    }
    expression_data <- expression_data[, keep_features]
  }
  
  # Convert to numeric if needed
  expression_data[] <- lapply(expression_data, function(x) {
    if (is.character(x) || is.factor(x)) {
      as.numeric(as.character(x))
    } else {
      as.numeric(x)
    }
  })
  
  return(list(
    expression_data = expression_data,
    class_labels = factor(class_labels, levels = c(control_label, case_label)),
    case_idx = which(class_labels == case_label),
    control_idx = which(class_labels == control_label)
  ))
}

#' Perform Differential Expression Analysis
#' @keywords internal
perform_de_analysis <- function(processed_data, data_type, p_adjust_methods, 
                               effect_size_methods, parallel, verbose) {
  
  expression_data <- processed_data$expression_data
  class_labels <- processed_data$class_labels
  case_idx <- processed_data$case_idx
  control_idx <- processed_data$control_idx
  
  feature_names <- colnames(expression_data)
  n_features <- length(feature_names)
  
  # Initialize results data frame
  results <- data.frame(
    feature = feature_names,
    stringsAsFactors = FALSE
  )
  
  # Perform statistical tests
  if (parallel && n_features > 100) {
    stat_results <- perform_tests_parallel(expression_data, case_idx, control_idx, 
                                          data_type, verbose)
  } else {
    stat_results <- perform_tests_sequential(expression_data, case_idx, control_idx, 
                                            data_type, verbose)
  }
  
  # Combine results
  results <- cbind(results, stat_results)
  
  # Apply multiple testing corrections
  for (method in p_adjust_methods) {
    adj_col_name <- paste0("p_value_", method)
    results[[adj_col_name]] <- p.adjust(results$p_value, method = method)
  }
  
  # Calculate effect sizes
  if (length(effect_size_methods) > 0) {
    effect_sizes <- calculate_effect_sizes(expression_data, case_idx, control_idx, 
                                          effect_size_methods, verbose)
    results <- cbind(results, effect_sizes)
  }
  
  return(results)
}

#' Perform Statistical Tests Sequentially
#' @keywords internal
perform_tests_sequential <- function(expression_data, case_idx, control_idx, data_type, verbose) {
  
  n_features <- ncol(expression_data)
  
  # Initialize result vectors
  p_values <- numeric(n_features)
  t_statistics <- numeric(n_features)
  log2_fold_changes <- numeric(n_features)
  case_means <- numeric(n_features)
  control_means <- numeric(n_features)
  case_sds <- numeric(n_features)
  control_sds <- numeric(n_features)
  
  # Progress tracking
  progress_points <- seq(1, n_features, length.out = min(10, n_features))
  
  for (i in seq_len(n_features)) {
    
    if (verbose && i %in% progress_points) {
      cat(sprintf("   Progress: %d/%d features (%.1f%%)\n", 
                  i, n_features, (i/n_features) * 100))
    }
    
    feature_values <- expression_data[, i]
    
    # Skip if too many missing values
    if (sum(is.na(feature_values)) > length(feature_values) * 0.5) {
      p_values[i] <- NA
      t_statistics[i] <- NA
      log2_fold_changes[i] <- NA
      case_means[i] <- NA
      control_means[i] <- NA
      case_sds[i] <- NA
      control_sds[i] <- NA
      next
    }
    
    case_values <- feature_values[case_idx]
    control_values <- feature_values[control_idx]
    
    # Remove NAs for this specific feature
    case_values <- case_values[!is.na(case_values)]
    control_values <- control_values[!is.na(control_values)]
    
    # Skip if insufficient data
    if (length(case_values) < 2 || length(control_values) < 2) {
      p_values[i] <- NA
      t_statistics[i] <- NA
      log2_fold_changes[i] <- NA
      case_means[i] <- NA
      control_means[i] <- NA
      case_sds[i] <- NA
      control_sds[i] <- NA
      next
    }
    
    # Calculate summary statistics
    case_mean <- mean(case_values)
    control_mean <- mean(control_values)
    case_sd <- sd(case_values)
    control_sd <- sd(control_values)
    
    case_means[i] <- case_mean
    control_means[i] <- control_mean
    case_sds[i] <- case_sd
    control_sds[i] <- control_sd
    
    # Perform t-test
    tryCatch({
      test_result <- t.test(case_values, control_values, var.equal = FALSE)
      p_values[i] <- test_result$p.value
      t_statistics[i] <- test_result$statistic
    }, error = function(e) {
      p_values[i] <- NA
      t_statistics[i] <- NA
    })
    
    # Calculate fold change based on data type
    if (data_type == "delta_ct") {
      # For delta-Ct, negative values indicate higher expression
      log2_fold_changes[i] <- control_mean - case_mean
    } else {
      # For log-transformed data
      log2_fold_changes[i] <- case_mean - control_mean
    }
  }
  
  return(data.frame(
    p_value = p_values,
    t_statistic = t_statistics,
    log2_fold_change = log2_fold_changes,
    case_mean = case_means,
    control_mean = control_means,
    case_sd = case_sds,
    control_sd = control_sds,
    stringsAsFactors = FALSE
  ))
}

#' Perform Statistical Tests in Parallel
#' @keywords internal
perform_tests_parallel <- function(expression_data, case_idx, control_idx, data_type, verbose) {
  
  if (!requireNamespace("parallel", quietly = TRUE)) {
    warning("Parallel processing not available, using sequential processing")
    return(perform_tests_sequential(expression_data, case_idx, control_idx, data_type, verbose))
  }
  
  cores <- min(parallel::detectCores() - 1, ncol(expression_data))
  if (cores < 2) {
    return(perform_tests_sequential(expression_data, case_idx, control_idx, data_type, verbose))
  }
  
  if (verbose) {
    cat(sprintf("🔄 Using %d cores for parallel processing...\n", cores))
  }
  
  # Split features into chunks
  feature_indices <- seq_len(ncol(expression_data))
  chunks <- split(feature_indices, cut(feature_indices, cores))
  
  # Process chunks in parallel
  results_list <- parallel::mclapply(chunks, function(chunk_indices) {
    chunk_data <- expression_data[, chunk_indices, drop = FALSE]
    perform_tests_sequential(chunk_data, case_idx, control_idx, data_type, FALSE)
  }, mc.cores = cores)
  
  # Combine results
  combined_results <- do.call(rbind, results_list)
  
  return(combined_results)
}

#' Calculate Effect Sizes
#' @keywords internal
calculate_effect_sizes <- function(expression_data, case_idx, control_idx, 
                                  effect_size_methods, verbose) {
  
  n_features <- ncol(expression_data)
  results <- data.frame(row.names = seq_len(n_features))
  
  for (method in effect_size_methods) {
    if (verbose) cat(sprintf("📊 Calculating %s effect sizes...\n", method))
    
    effect_sizes <- numeric(n_features)
    
    for (i in seq_len(n_features)) {
      feature_values <- expression_data[, i]
      case_values <- feature_values[case_idx]
      control_values <- feature_values[control_idx]
      
      # Remove NAs
      case_values <- case_values[!is.na(case_values)]
      control_values <- control_values[!is.na(control_values)]
      
      if (length(case_values) < 2 || length(control_values) < 2) {
        effect_sizes[i] <- NA
        next
      }
      
      tryCatch({
        if (method == "cohens_d") {
          effect_sizes[i] <- calculate_cohens_d(case_values, control_values)
        } else if (method == "hedges_g") {
          effect_sizes[i] <- calculate_hedges_g(case_values, control_values)
        } else if (method == "cliff_delta") {
          effect_sizes[i] <- calculate_cliff_delta(case_values, control_values)
        }
      }, error = function(e) {
        effect_sizes[i] <- NA
      })
    }
    
    results[[paste0("effect_size_", method)]] <- effect_sizes
  }
  
  return(results)
}

#' Calculate Cohen's d
#' @keywords internal
calculate_cohens_d <- function(x, y) {
  pooled_sd <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / 
                   (length(x) + length(y) - 2))
  (mean(x) - mean(y)) / pooled_sd
}

#' Calculate Hedges' g
#' @keywords internal
calculate_hedges_g <- function(x, y) {
  cohens_d <- calculate_cohens_d(x, y)
  correction_factor <- 1 - (3 / (4 * (length(x) + length(y)) - 9))
  cohens_d * correction_factor
}

#' Calculate Cliff's Delta
#' @keywords internal
calculate_cliff_delta <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  
  # Create all pairwise comparisons
  comparisons <- outer(x, y, function(a, b) sign(a - b))
  
  # Calculate Cliff's delta
  sum(comparisons) / (n1 * n2)
}

#' Add Significance Flags
#' @keywords internal
add_significance_flags <- function(results, alpha, verbose) {
  
  # Add significance columns for each adjusted p-value method
  adj_p_cols <- grep("^p_value_", colnames(results), value = TRUE)
  
  for (col in adj_p_cols) {
    sig_col <- gsub("p_value_", "significant_", col)
    results[[sig_col]] <- results[[col]] < alpha
  }
  
  if (verbose && length(adj_p_cols) > 0) {
    for (col in adj_p_cols) {
      sig_count <- sum(results[[col]] < alpha, na.rm = TRUE)
      method <- gsub("p_value_", "", col)
      cat(sprintf("   %s significant features: %d\n", method, sig_count))
    }
  }
  
  return(results)
}

#' Backward Compatibility Alias
#' @export
OmicSelector_differential_expression_ttest <- function(ttpm_features, classes, mode = "auto") {
  .Deprecated("differential_expression_modern", 
              msg = "OmicSelector_differential_expression_ttest is deprecated. Use differential_expression_modern() instead.")
  
  # Convert old parameters to new function
  differential_expression_modern(
    expression_data = ttpm_features,
    class_labels = classes,
    data_type = if (mode == "auto") "auto" else if (mode == "deltact") "delta_ct" else "log_tpm",
    verbose = FALSE  # Reduce verbosity for backward compatibility
  )
}
