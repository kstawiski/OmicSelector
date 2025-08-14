#' Modular OmicSelector Integration
#'
#' @description
#' High-level API that integrates the modular feature selection framework
#' with existing OmicSelector functionality. Provides seamless transition
#' from the comprehensive system to the new modular architecture.
#'
#' @author OmicSelector Team

# =============================================================================
# HIGH-LEVEL API
# =============================================================================

#' Modular Feature Selection
#'
#' @description
#' Execute feature selection using the modular framework. Supports both
#' individual methods and multiple methods with automatic parallelization.
#'
#' @param data Data frame with class column as first column
#' @param methods Vector of method IDs or method names to execute
#' @param max_features Maximum features per method
#' @param config Named list of method-specific configurations
#' @param use_smote Apply SMOTE balancing
#' @param parallel Use parallel execution for multiple methods
#' @param timeout_sec Timeout per method
#' @param initialize_framework Initialize framework if not already done
#' @param progress Show progress information
#'
#' @return List of method results or single result if one method
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(iris)
#' iris_data <- data.frame(class = iris$Species, iris[1:4])
#' 
#' # Single method
#' result <- modular_feature_selection(iris_data, methods = 1, max_features = 3)
#' print(result$features)
#' 
#' # Multiple methods
#' results <- modular_feature_selection(
#'   iris_data, 
#'   methods = c(1, 2, 3), 
#'   max_features = 3,
#'   parallel = TRUE
#' )
#' 
#' # Using method names
#' results <- modular_feature_selection(
#'   iris_data,
#'   methods = c("Student's t-test", "Correlation", "Random Forest"),
#'   max_features = 3
#' )
#' }
modular_feature_selection <- function(data,
                                     methods = c(1, 2, 3, 4, 5),
                                     max_features = 20,
                                     config = list(),
                                     use_smote = FALSE,
                                     parallel = TRUE,
                                     timeout_sec = NULL,
                                     initialize_framework = TRUE,
                                     progress = TRUE) {
  
  # Initialize framework if requested
  if (initialize_framework) {
    framework_status <- initialize_fs_framework(verbose = progress)
    if (!framework_status$success) {
      stop("Failed to initialize feature selection framework")
    }
  }
  
  # Convert method names to IDs if necessary
  method_ids <- resolve_method_identifiers(methods)
  
  if (length(method_ids) == 0) {
    stop("No valid methods specified")
  }
  
  # Validate data
  if (!is.data.frame(data) || ncol(data) < 2) {
    stop("Data must be a data frame with at least 2 columns (class + features)")
  }
  
  if (progress) {
    cat("Executing modular feature selection...\n")
    cat("Methods:", paste(method_ids, collapse = ", "), "\n")
    cat("Max features per method:", max_features, "\n")
    cat("Data shape:", nrow(data), "samples x", ncol(data) - 1, "features\n")
    if (use_smote) cat("SMOTE balancing: enabled\n")
    cat("\n")
  }
  
  # Execute methods
  if (length(method_ids) == 1) {
    # Single method execution
    result <- execute_fs_method(
      method_id = method_ids[1],
      data = data,
      config = config,
      max_features = max_features,
      use_smote = use_smote,
      timeout_sec = timeout_sec
    )
    
    if (progress && !is.null(result)) {
      cat("✓ Method", method_ids[1], "completed successfully\n")
      cat("Selected features:", paste(result$features, collapse = ", "), "\n")
    }
    
    return(result)
    
  } else {
    # Multiple methods execution
    results <- execute_multiple_fs_methods(
      method_ids = method_ids,
      data = data,
      config = config,
      max_features = max_features,
      use_smote = use_smote,
      parallel = parallel,
      timeout_sec = timeout_sec,
      progress = progress
    )
    
    return(results)
  }
}

#' Enhanced Feature Selection with Analysis
#'
#' @description
#' Run feature selection with additional analysis including feature overlap,
#' stability analysis, and performance summaries.
#'
#' @param data Training data
#' @param methods Methods to execute
#' @param max_features Max features per method
#' @param config Method configurations
#' @param use_smote SMOTE balancing
#' @param analysis_options List of analysis options
#' @param parallel Parallel execution
#' @param progress Show progress
#'
#' @return Enhanced results with analysis
#' @export
enhanced_modular_feature_selection <- function(data,
                                              methods = c(1, 2, 3, 4, 5),
                                              max_features = 20,
                                              config = list(),
                                              use_smote = FALSE,
                                              analysis_options = list(
                                                feature_overlap = TRUE,
                                                stability_analysis = TRUE,
                                                performance_summary = TRUE,
                                                consensus_features = TRUE
                                              ),
                                              parallel = TRUE,
                                              progress = TRUE) {
  
  # Run core feature selection
  results <- modular_feature_selection(
    data = data,
    methods = methods,
    max_features = max_features,
    config = config,
    use_smote = use_smote,
    parallel = parallel,
    timeout_sec = NULL,
    initialize_framework = TRUE,
    progress = progress
  )
  
  if (is.null(results) || length(results) == 0) {
    warning("No successful feature selection results")
    return(NULL)
  }
  
  # Filter successful results
  successful_results <- Filter(function(x) !is.null(x), results)
  
  if (length(successful_results) == 0) {
    warning("No methods completed successfully")
    return(list(results = results, analysis = NULL))
  }
  
  if (progress) {
    cat("\nPerforming feature selection analysis...\n")
  }
  
  # Initialize analysis results
  analysis <- list()
  
  # Feature overlap analysis
  if (analysis_options$feature_overlap) {
    analysis$feature_overlap <- analyze_feature_overlap(successful_results, progress = progress)
  }
  
  # Stability analysis
  if (analysis_options$stability_analysis) {
    analysis$stability <- analyze_feature_stability(successful_results, progress = progress)
  }
  
  # Performance summary
  if (analysis_options$performance_summary) {
    analysis$performance_summary <- generate_performance_summary(successful_results, progress = progress)
  }
  
  # Consensus features
  if (analysis_options$consensus_features) {
    analysis$consensus_features <- identify_consensus_features(successful_results, progress = progress)
  }
  
  if (progress) {
    cat("Analysis complete!\n\n")
    
    # Print summary
    if (!is.null(analysis$consensus_features)) {
      cat("Top consensus features:\n")
      top_consensus <- head(analysis$consensus_features$features, 10)
      for (i in seq_along(top_consensus)) {
        feat <- top_consensus[i]
        count <- analysis$consensus_features$selection_counts[feat]
        cat(sprintf("  %d. %s (selected by %d/%d methods)\n", 
                   i, feat, count, length(successful_results)))
      }
    }
  }
  
  return(list(
    results = results,
    analysis = analysis,
    summary = list(
      n_methods_attempted = length(methods),
      n_methods_successful = length(successful_results),
      n_features_input = ncol(data) - 1,
      execution_timestamp = Sys.time()
    )
  ))
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Resolve Method Identifiers
#'
#' @description
#' Convert method names to IDs or validate existing IDs.
#'
#' @param methods Vector of method IDs or names
#' @return Vector of method IDs
resolve_method_identifiers <- function(methods) {
  
  if (length(methods) == 0) {
    return(integer(0))
  }
  
  # Get available methods
  available_methods <- list_fs_methods()
  
  if (nrow(available_methods) == 0) {
    warning("No methods available in registry")
    return(integer(0))
  }
  
  resolved_ids <- integer(0)
  
  for (method in methods) {
    if (is.numeric(method)) {
      # Already an ID, validate it exists
      if (method %in% available_methods$id) {
        resolved_ids <- c(resolved_ids, method)
      } else {
        warning("Method ID ", method, " not found in registry")
      }
    } else if (is.character(method)) {
      # Method name, convert to ID
      matching_row <- available_methods[available_methods$name == method, ]
      if (nrow(matching_row) > 0) {
        resolved_ids <- c(resolved_ids, matching_row$id[1])
      } else {
        warning("Method name '", method, "' not found in registry")
      }
    } else {
      warning("Invalid method identifier: ", method)
    }
  }
  
  return(unique(resolved_ids))
}

#' Analyze Feature Overlap
#'
#' @description
#' Analyze overlap between features selected by different methods.
#'
#' @param results List of method results
#' @param progress Show progress
#' @return Feature overlap analysis
analyze_feature_overlap <- function(results, progress = FALSE) {
  
  if (progress) cat("  - Computing feature overlap matrix...\n")
  
  # Extract feature lists
  feature_lists <- lapply(results, function(r) {
    if (!is.null(r) && !is.null(r$features)) r$features else character(0)
  })
  
  # Remove empty results
  feature_lists <- feature_lists[sapply(feature_lists, length) > 0]
  
  if (length(feature_lists) < 2) {
    return(list(
      overlap_matrix = NULL,
      jaccard_matrix = NULL,
      message = "Insufficient results for overlap analysis"
    ))
  }
  
  # Create overlap matrix
  n_methods <- length(feature_lists)
  method_names <- names(feature_lists)
  if (is.null(method_names)) {
    method_names <- paste("Method", seq_len(n_methods))
  }
  
  overlap_matrix <- matrix(0, nrow = n_methods, ncol = n_methods,
                          dimnames = list(method_names, method_names))
  jaccard_matrix <- matrix(0, nrow = n_methods, ncol = n_methods,
                          dimnames = list(method_names, method_names))
  
  for (i in seq_len(n_methods)) {
    for (j in seq_len(n_methods)) {
      features_i <- feature_lists[[i]]
      features_j <- feature_lists[[j]]
      
      intersection <- length(intersect(features_i, features_j))
      union_size <- length(union(features_i, features_j))
      
      overlap_matrix[i, j] <- intersection
      jaccard_matrix[i, j] <- if (union_size > 0) intersection / union_size else 0
    }
  }
  
  return(list(
    overlap_matrix = overlap_matrix,
    jaccard_matrix = jaccard_matrix,
    feature_lists = feature_lists
  ))
}

#' Analyze Feature Stability
#'
#' @description
#' Analyze stability of feature selection across methods.
#'
#' @param results List of method results
#' @param progress Show progress
#' @return Stability analysis
analyze_feature_stability <- function(results, progress = FALSE) {
  
  if (progress) cat("  - Analyzing feature stability...\n")
  
  # Extract all unique features
  all_features <- unique(unlist(lapply(results, function(r) {
    if (!is.null(r) && !is.null(r$features)) r$features else character(0)
  })))
  
  if (length(all_features) == 0) {
    return(list(message = "No features found in results"))
  }
  
  # Count selections per feature
  selection_counts <- sapply(all_features, function(feature) {
    sum(sapply(results, function(r) {
      if (!is.null(r) && !is.null(r$features)) {
        feature %in% r$features
      } else {
        FALSE
      }
    }))
  })
  
  # Calculate stability metrics
  n_methods <- length(Filter(function(x) !is.null(x), results))
  selection_frequency <- selection_counts / n_methods
  
  # Identify stable features (selected by majority of methods)
  stability_threshold <- 0.5
  stable_features <- all_features[selection_frequency >= stability_threshold]
  
  return(list(
    all_features = all_features,
    selection_counts = selection_counts,
    selection_frequency = selection_frequency,
    stable_features = stable_features,
    stability_threshold = stability_threshold,
    n_methods = n_methods
  ))
}

#' Generate Performance Summary
#'
#' @description
#' Generate summary of method performance and execution statistics.
#'
#' @param results List of method results
#' @param progress Show progress
#' @return Performance summary
generate_performance_summary <- function(results, progress = FALSE) {
  
  if (progress) cat("  - Generating performance summary...\n")
  
  # Extract method information
  method_info <- lapply(results, function(r) {
    if (is.null(r) || is.null(r$metadata)) {
      return(list(
        method_name = "Unknown",
        execution_time = NA,
        n_features_selected = 0,
        success = FALSE
      ))
    }
    
    list(
      method_name = r$metadata$method_name %||% "Unknown",
      execution_time = r$metadata$execution_time %||% NA,
      n_features_selected = length(r$features %||% character(0)),
      success = TRUE
    )
  })
  
  # Create summary data frame
  summary_df <- data.frame(
    method = sapply(method_info, function(x) x$method_name),
    execution_time_sec = sapply(method_info, function(x) x$execution_time),
    features_selected = sapply(method_info, function(x) x$n_features_selected),
    success = sapply(method_info, function(x) x$success),
    stringsAsFactors = FALSE
  )
  
  # Calculate summary statistics
  successful_methods <- summary_df[summary_df$success, ]
  
  summary_stats <- list(
    total_methods = nrow(summary_df),
    successful_methods = nrow(successful_methods),
    success_rate = nrow(successful_methods) / nrow(summary_df),
    avg_execution_time = mean(successful_methods$execution_time_sec, na.rm = TRUE),
    avg_features_selected = mean(successful_methods$features_selected, na.rm = TRUE),
    total_execution_time = sum(successful_methods$execution_time_sec, na.rm = TRUE)
  )
  
  return(list(
    method_details = summary_df,
    summary_statistics = summary_stats
  ))
}

#' Identify Consensus Features
#'
#' @description
#' Identify features that are consistently selected across methods.
#'
#' @param results List of method results
#' @param progress Show progress
#' @return Consensus features analysis
identify_consensus_features <- function(results, progress = FALSE) {
  
  if (progress) cat("  - Identifying consensus features...\n")
  
  # Use stability analysis as base
  stability_result <- analyze_feature_stability(results, progress = FALSE)
  
  if (is.null(stability_result$selection_counts)) {
    return(list(message = "No features found for consensus analysis"))
  }
  
  # Sort features by selection frequency
  sorted_features <- names(sort(stability_result$selection_counts, decreasing = TRUE))
  
  # Create consensus ranking
  consensus_ranking <- data.frame(
    feature = sorted_features,
    selection_count = stability_result$selection_counts[sorted_features],
    selection_frequency = stability_result$selection_frequency[sorted_features],
    rank = seq_along(sorted_features),
    stringsAsFactors = FALSE
  )
  
  # Identify different consensus levels
  high_consensus <- consensus_ranking[consensus_ranking$selection_frequency >= 0.75, ]
  medium_consensus <- consensus_ranking[consensus_ranking$selection_frequency >= 0.5 & 
                                       consensus_ranking$selection_frequency < 0.75, ]
  low_consensus <- consensus_ranking[consensus_ranking$selection_frequency >= 0.25 & 
                                    consensus_ranking$selection_frequency < 0.5, ]
  
  return(list(
    features = sorted_features,
    selection_counts = stability_result$selection_counts[sorted_features],
    consensus_ranking = consensus_ranking,
    high_consensus = high_consensus$feature,
    medium_consensus = medium_consensus$feature,
    low_consensus = low_consensus$feature,
    consensus_levels = list(
      high = nrow(high_consensus),
      medium = nrow(medium_consensus),
      low = nrow(low_consensus)
    )
  ))
}
