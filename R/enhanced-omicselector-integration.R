#' Enhanced OmicSelector Integration
#'
#' @description
#' Integration module that enhances the existing modernized OmicSelector functions
#' with comprehensive feature selection capabilities while maintaining full
#' backward compatibility and minimal dependencies.
#'
#' @author OmicSelector Team
#' @export

# =============================================================================
# ENHANCED OMICSELECTOR INTEGRATION
# =============================================================================

#' Enhanced OmicSelector with Comprehensive Feature Selection
#'
#' @description
#' Integrates the comprehensive feature selection module with the modernized
#' OmicSelector pipeline, providing access to 70+ feature selection methods
#' while maintaining optimal performance and minimal dependencies.
#'
#' @param wd Working directory path
#' @param m Methods vector (1-70 for comprehensive selection)
#' @param max_iterations Maximum iterations for iterative methods
#' @param prefer_no_features Maximum features to select per method
#' @param use_comprehensive Whether to use comprehensive FS (default: TRUE)
#' @param enable_advanced_methods Enable methods 31-70 (default: TRUE)
#' @param timeout_sec Timeout per method in seconds
#' @param parallel Use parallel processing where available
#' @param ... Additional parameters passed to original function
#'
#' @return Enhanced OmicSelector results with comprehensive feature selection
#' @export
#'
#' @examples
#' \dontrun{
#' # Use comprehensive feature selection (all 70 methods)
#' results <- enhanced_omicselector(
#'   wd = "path/to/data",
#'   m = 1:20,  # Use first 20 methods
#'   prefer_no_features = 15,
#'   use_comprehensive = TRUE
#' )
#'
#' # Access selected features by method
#' print(names(results$selected_features))
#' 
#' # Get best performing methods
#' best_methods <- get_best_fs_methods(results, top_n = 5)
#' }
enhanced_omicselector <- function(wd = getwd(),
                                 m = 1:20,
                                 max_iterations = 10,
                                 prefer_no_features = 20,
                                 use_comprehensive = TRUE,
                                 enable_advanced_methods = TRUE,
                                 timeout_sec = 1800,
                                 parallel = TRUE,
                                 ...) {
  
  # Log function start
  if (exists("log_info")) {
    log_info("Starting Enhanced OmicSelector with comprehensive feature selection")
  } else {
    cat("[INFO] Starting Enhanced OmicSelector with comprehensive feature selection\n")
  }
  
  # Set working directory
  original_wd <- getwd()
  setwd(wd)
  on.exit(setwd(original_wd))
  
  # Validate input parameters
  if (use_comprehensive && max(m) > 70) {
    warning("Method numbers > 70 not supported in comprehensive mode. Limiting to 1-70.")
    m <- m[m <= 70]
  }
  
  if (!enable_advanced_methods) {
    # Limit to basic methods (1-30) for faster execution
    m <- m[m <= 30]
  }
  
  # Load training data for feature selection
  train_data <- load_training_data_safe()
  if (is.null(train_data)) {
    stop("Failed to load training data from working directory")
  }
  
  # Initialize enhanced results object
  enhanced_results <- list(
    formulas = list(),
    selected_features = list(),
    method_performance = list(),
    execution_times = list(),
    errors = list(),
    summary = list(),
    metadata = list(
      enhanced_version = TRUE,
      comprehensive_fs = use_comprehensive,
      methods_requested = m,
      advanced_methods_enabled = enable_advanced_methods,
      total_methods = length(m),
      execution_start = Sys.time()
    )
  )
  
  # Execute feature selection based on mode
  if (use_comprehensive) {
    
    cat("Running comprehensive feature selection with", length(m), "methods...\n")
    
    # Use comprehensive feature selection module
    fs_results <- run_comprehensive_feature_selection(
      data = train_data,
      methods = m,
      config = list(
        prefer_no_features = prefer_no_features,
        max_iterations = max_iterations,
        timeout_sec = timeout_sec,
        parallel = parallel
      ),
      use_smote = TRUE,
      prefer_no_features = prefer_no_features,
      timeout_sec = timeout_sec,
      parallel = parallel
    )
    
    # Integrate comprehensive results
    enhanced_results$formulas <- fs_results$formulas
    enhanced_results$selected_features <- fs_results$selected_features
    enhanced_results$method_performance <- fs_results$method_performance
    enhanced_results$execution_times <- fs_results$execution_times
    enhanced_results$errors <- fs_results$errors
    
  } else {
    
    cat("Running standard feature selection...\n")
    
    # Use original/improved OmicSelector function
    if (exists("OmicSelector_OmicSelector")) {
      standard_results <- OmicSelector_OmicSelector(
        wd = wd,
        m = m,
        max_iterations = max_iterations,
        prefer_no_features = prefer_no_features,
        timeout_sec = timeout_sec,
        ...
      )
      
      # Convert standard results to enhanced format
      enhanced_results <- convert_standard_to_enhanced(standard_results, enhanced_results)
    } else {
      stop("Standard OmicSelector function not available")
    }
  }
  
  # Compute summary statistics
  enhanced_results$summary <- compute_enhanced_summary(enhanced_results)
  enhanced_results$metadata$execution_end <- Sys.time()
  enhanced_results$metadata$total_execution_time <- 
    as.numeric(enhanced_results$metadata$execution_end - enhanced_results$metadata$execution_start, units = "secs")
  
  # Set class for S3 methods
  class(enhanced_results) <- c("enhanced_omicselector", "omicselector_results", "list")
  
  # Print summary
  print_enhanced_summary(enhanced_results)
  
  return(enhanced_results)
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Load Training Data Safely
#'
#' @return Training data frame or NULL if failed
#' @keywords internal
load_training_data_safe <- function() {
  
  # Try to load from various file names
  possible_files <- c(
    "mixed_train.csv",
    "train.csv", 
    "training.csv",
    "mixed_training.csv"
  )
  
  for (file_name in possible_files) {
    if (file.exists(file_name)) {
      tryCatch({
        data <- utils::read.csv(file_name, stringsAsFactors = FALSE)
        if (ncol(data) > 1 && nrow(data) > 0) {
          cat("Loaded training data from:", file_name, "\n")
          return(data)
        }
      }, error = function(e) {
        # Continue to next file
      })
    }
  }
  
  return(NULL)
}

#' Convert Standard Results to Enhanced Format
#'
#' @param standard_results Results from standard OmicSelector
#' @param enhanced_results Enhanced results template
#' @return Enhanced results object
#' @keywords internal
convert_standard_to_enhanced <- function(standard_results, enhanced_results) {
  
  if (is.list(standard_results)) {
    
    # Extract formulas if available
    if ("formulas" %in% names(standard_results)) {
      enhanced_results$formulas <- standard_results$formulas
    }
    
    # Extract selected features from formulas
    if (length(enhanced_results$formulas) > 0) {
      enhanced_results$selected_features <- extract_features_from_formulas(enhanced_results$formulas)
    }
    
    # Copy other available results
    for (field in c("method_performance", "execution_times", "errors", "summary")) {
      if (field %in% names(standard_results)) {
        enhanced_results[[field]] <- standard_results[[field]]
      }
    }
  }
  
  return(enhanced_results)
}

#' Extract Features from Formula List
#'
#' @param formulas List of formulas
#' @return List of feature vectors
#' @keywords internal
extract_features_from_formulas <- function(formulas) {
  
  selected_features <- list()
  
  for (method_name in names(formulas)) {
    formula_obj <- formulas[[method_name]]
    
    if (inherits(formula_obj, "formula")) {
      # Extract feature names from formula
      formula_str <- as.character(formula_obj)[3]  # RHS of formula
      
      if (formula_str != "1") {  # Not empty formula
        # Split by + and clean up
        features <- strsplit(formula_str, " \\+ ")[[1]]
        features <- gsub("`", "", features)  # Remove backticks
        features <- trimws(features)  # Remove whitespace
        
        selected_features[[method_name]] <- features
      }
    }
  }
  
  return(selected_features)
}

#' Compute Enhanced Summary Statistics
#'
#' @param enhanced_results Enhanced results object
#' @return Summary statistics list
#' @keywords internal
compute_enhanced_summary <- function(enhanced_results) {
  
  summary_stats <- list(
    total_methods_attempted = enhanced_results$metadata$total_methods,
    successful_methods = length(enhanced_results$selected_features),
    failed_methods = length(enhanced_results$errors),
    total_unique_features = length(unique(unlist(enhanced_results$selected_features))),
    average_features_per_method = 0,
    execution_time_seconds = enhanced_results$metadata$total_execution_time %||% 0,
    feature_frequency = list(),
    top_methods_by_features = list(),
    performance_summary = list()
  )
  
  # Calculate average features per method
  if (length(enhanced_results$selected_features) > 0) {
    feature_counts <- sapply(enhanced_results$selected_features, length)
    summary_stats$average_features_per_method <- mean(feature_counts, na.rm = TRUE)
    
    # Feature frequency analysis
    all_features <- unlist(enhanced_results$selected_features)
    if (length(all_features) > 0) {
      feature_freq <- table(all_features)
      summary_stats$feature_frequency <- as.list(sort(feature_freq, decreasing = TRUE))
    }
    
    # Top methods by number of features
    summary_stats$top_methods_by_features <- names(sort(feature_counts, decreasing = TRUE))
  }
  
  return(summary_stats)
}

#' Print Enhanced Summary
#'
#' @param enhanced_results Enhanced results object
#' @keywords internal
print_enhanced_summary <- function(enhanced_results) {
  
  cat("\n")
  cat("=== Enhanced OmicSelector Results Summary ===\n")
  cat("Comprehensive Feature Selection:", enhanced_results$metadata$comprehensive_fs, "\n")
  cat("Methods Attempted:", enhanced_results$summary$total_methods_attempted, "\n")
  cat("Successful Methods:", enhanced_results$summary$successful_methods, "\n")
  cat("Failed Methods:", enhanced_results$summary$failed_methods, "\n")
  cat("Total Unique Features Found:", enhanced_results$summary$total_unique_features, "\n")
  cat("Average Features per Method:", round(enhanced_results$summary$average_features_per_method, 2), "\n")
  cat("Total Execution Time:", round(enhanced_results$summary$execution_time_seconds, 2), "seconds\n")
  
  # Show top features if available
  if (length(enhanced_results$summary$feature_frequency) > 0) {
    cat("\nTop 10 Most Frequently Selected Features:\n")
    top_features <- head(enhanced_results$summary$feature_frequency, 10)
    for (i in seq_along(top_features)) {
      cat(sprintf("  %s: %d methods\n", names(top_features)[i], top_features[i]))
    }
  }
  
  # Show top methods if available
  if (length(enhanced_results$summary$top_methods_by_features) > 0) {
    cat("\nTop 5 Methods by Number of Features:\n")
    top_methods <- head(enhanced_results$summary$top_methods_by_features, 5)
    for (method in top_methods) {
      n_features <- length(enhanced_results$selected_features[[method]])
      cat(sprintf("  %s: %d features\n", method, n_features))
    }
  }
  
  cat("\n")
}

# =============================================================================
# S3 METHODS FOR ENHANCED RESULTS
# =============================================================================

#' Print Method for Enhanced OmicSelector Results
#'
#' @param x Enhanced OmicSelector results
#' @param ... Additional arguments
#' @method print enhanced_omicselector
#' @export
print.enhanced_omicselector <- function(x, ...) {
  print_enhanced_summary(x)
}

#' Summary Method for Enhanced OmicSelector Results
#'
#' @param object Enhanced OmicSelector results
#' @param ... Additional arguments
#' @method summary enhanced_omicselector
#' @export
summary.enhanced_omicselector <- function(object, ...) {
  
  cat("Enhanced OmicSelector Results Summary\n")
  cat("=====================================\n\n")
  
  # Basic statistics
  cat("Execution Details:\n")
  cat("- Comprehensive Feature Selection:", object$metadata$comprehensive_fs, "\n")
  cat("- Advanced Methods Enabled:", object$metadata$advanced_methods_enabled, "\n")
  cat("- Methods Requested:", paste(head(object$metadata$methods_requested, 10), collapse = ", "))
  if (length(object$metadata$methods_requested) > 10) cat("...")
  cat("\n\n")
  
  # Performance statistics
  cat("Performance Statistics:\n")
  cat("- Methods Attempted:", object$summary$total_methods_attempted, "\n")
  cat("- Successful Methods:", object$summary$successful_methods, "\n")
  cat("- Success Rate:", round(100 * object$summary$successful_methods / object$summary$total_methods_attempted, 1), "%\n")
  cat("- Execution Time:", round(object$summary$execution_time_seconds, 2), "seconds\n\n")
  
  # Feature statistics
  cat("Feature Selection Statistics:\n")
  cat("- Total Unique Features:", object$summary$total_unique_features, "\n")
  cat("- Average Features per Method:", round(object$summary$average_features_per_method, 2), "\n")
  
  return(invisible(object$summary))
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

#' Get Best Feature Selection Methods
#'
#' @description
#' Analyze enhanced OmicSelector results to identify the best performing
#' feature selection methods based on various criteria.
#'
#' @param enhanced_results Enhanced OmicSelector results
#' @param top_n Number of top methods to return
#' @param criteria Ranking criteria ("features", "stability", "performance")
#' @return Data frame with top methods and their metrics
#' @export
get_best_fs_methods <- function(enhanced_results, top_n = 5, criteria = "features") {
  
  if (!inherits(enhanced_results, "enhanced_omicselector")) {
    stop("Input must be enhanced_omicselector results object")
  }
  
  if (length(enhanced_results$selected_features) == 0) {
    warning("No successful feature selection methods found")
    return(data.frame())
  }
  
  # Create analysis data frame
  methods_df <- data.frame(
    method = names(enhanced_results$selected_features),
    n_features = sapply(enhanced_results$selected_features, length),
    execution_time = sapply(names(enhanced_results$selected_features), function(m) {
      enhanced_results$execution_times[[m]] %||% NA
    }),
    stringsAsFactors = FALSE
  )
  
  # Add stability score (based on feature frequency)
  all_features <- unlist(enhanced_results$selected_features)
  feature_freq <- table(all_features)
  
  methods_df$stability_score <- sapply(enhanced_results$selected_features, function(features) {
    if (length(features) == 0) return(0)
    mean(feature_freq[features], na.rm = TRUE)
  })
  
  # Rank by criteria
  if (criteria == "features") {
    methods_df <- methods_df[order(methods_df$n_features, decreasing = TRUE), ]
  } else if (criteria == "stability") {
    methods_df <- methods_df[order(methods_df$stability_score, decreasing = TRUE), ]
  } else if (criteria == "performance") {
    # Combined score: features + stability - execution_time
    methods_df$combined_score <- scale(methods_df$n_features)[,1] + 
                                 scale(methods_df$stability_score)[,1] - 
                                 scale(methods_df$execution_time, center = TRUE, scale = TRUE)[,1]
    methods_df <- methods_df[order(methods_df$combined_score, decreasing = TRUE), ]
  }
  
  # Return top N
  result <- head(methods_df, top_n)
  rownames(result) <- NULL
  
  return(result)
}

#' Get Feature Overlap Analysis
#'
#' @description
#' Analyze overlap between features selected by different methods.
#'
#' @param enhanced_results Enhanced OmicSelector results
#' @param min_methods Minimum methods that must select a feature
#' @return Data frame with feature overlap analysis
#' @export
get_feature_overlap_analysis <- function(enhanced_results, min_methods = 2) {
  
  if (!inherits(enhanced_results, "enhanced_omicselector")) {
    stop("Input must be enhanced_omicselector results object")
  }
  
  if (length(enhanced_results$selected_features) == 0) {
    return(data.frame())
  }
  
  # Get feature frequency
  all_features <- unlist(enhanced_results$selected_features)
  feature_freq <- table(all_features)
  
  # Filter by minimum methods
  stable_features <- feature_freq[feature_freq >= min_methods]
  
  if (length(stable_features) == 0) {
    return(data.frame())
  }
  
  # Create analysis data frame
  overlap_df <- data.frame(
    feature = names(stable_features),
    frequency = as.numeric(stable_features),
    percentage = round(100 * as.numeric(stable_features) / length(enhanced_results$selected_features), 2),
    stringsAsFactors = FALSE
  )
  
  # Add methods that selected each feature
  overlap_df$selected_by <- sapply(overlap_df$feature, function(feat) {
    methods <- names(enhanced_results$selected_features)[
      sapply(enhanced_results$selected_features, function(x) feat %in% x)
    ]
    paste(methods, collapse = ", ")
  })
  
  # Order by frequency
  overlap_df <- overlap_df[order(overlap_df$frequency, decreasing = TRUE), ]
  rownames(overlap_df) <- NULL
  
  return(overlap_df)
}

# =============================================================================
# BACKWARD COMPATIBILITY ALIASES
# =============================================================================

#' @rdname enhanced_omicselector
#' @export
comprehensive_omicselector <- enhanced_omicselector

#' @rdname enhanced_omicselector  
#' @export
omicselector_enhanced <- enhanced_omicselector
