#' Improved OmicSelector_OmicSelector Function
#'
#' @description
#' Completely rewritten and optimized version of the original OmicSelector_OmicSelector
#' function. This version is specifically optimized for miRNA Case/Control analysis
#' with improved performance, error handling, and modern R practices.
#'
#' Key improvements:
#' - Optimized for "hsa_" prefixed miRNA features
#' - Case/Control binary classification focus  
#' - Modern error handling and validation
#' - Parallel processing optimization
#' - Memory-efficient operations
#' - Comprehensive logging and progress tracking
#' - Timeout protection for long-running methods
#' - Intermediate result saving
#' - Modern S3 class system
#'
#' @param wd Working directory with data files (mixed_train.csv, mixed_test.csv, mixed_validation.csv)
#' @param m Methods of feature selection to be performed (integer vector 1-70)
#' @param max_iterations Maximum number of iterations in selected methods
#' @param code_path Path to external scripts (legacy parameter, maintained for compatibility)
#' @param register_parallel Whether to use parallel processing
#' @param clx Computing cluster (legacy parameter, auto-managed now)
#' @param stamp Timestamp for marking output files
#' @param prefer_no_features Maximum number of features that can be selected
#' @param conda_path Path to conda binary (legacy parameter)
#' @param debug Enable debug mode with additional output
#' @param timeout_sec Timeout in seconds for each method
#' @param type Mode parameter (legacy parameter, auto-detected now)
#' @param class_column Name of the class column (default: "Class")
#' @param case_label Label for positive cases (default: "Case")
#' @param control_label Label for negative controls (default: "Control")
#' @param feature_prefix Prefix for feature columns (default: "hsa_")
#' @param p_threshold P-value threshold for significance (default: 0.05)
#' @param fc_threshold Fold change threshold (default: 1.0)
#' @param use_smote Whether to apply SMOTE for class balancing (default: TRUE)
#' @param save_intermediate Save intermediate results (default: TRUE)
#' @param verbose Enable verbose logging (default: TRUE)
#'
#' @return Modern OmicSelector results object (S3 class) with selected features,
#'   performance metrics, differential expression results, and execution metadata.
#'   Backward compatible with legacy code expecting formula lists.
#'
#' @examples
#' \dontrun{
#' # Basic usage (backward compatible)
#' results <- OmicSelector_OmicSelector(
#'   wd = "path/to/data",
#'   m = c(1, 2, 4, 11, 17)
#' )
#'
#' # Access selected features (new way)
#' selected_features <- results$selected_features
#' 
#' # Access formulas (legacy compatible)
#' formulas <- results$formulas
#'
#' # With custom parameters for miRNA analysis
#' results <- OmicSelector_OmicSelector(
#'   wd = "path/to/mirna_data",
#'   m = c(1, 2, 4, 11, 13, 17),
#'   max_iterations = 20,
#'   prefer_no_features = 15,
#'   p_threshold = 0.01,
#'   fc_threshold = 1.5,
#'   case_label = "Case",
#'   control_label = "Control",
#'   feature_prefix = "hsa_"
#' )
#' }
#'
#' @export
OmicSelector_OmicSelector <- function(wd = getwd(), 
                                     m = c(1:20),
                                     max_iterations = 10, 
                                     code_path = system.file("extdata", "", package = "OmicSelector"),
                                     register_parallel = TRUE, 
                                     clx = NULL, 
                                     stamp = as.numeric(Sys.time()),
                                     prefer_no_features = 20, 
                                     conda_path = "/home/konrad/anaconda3/bin/conda", 
                                     debug = FALSE,
                                     timeout_sec = 1800, 
                                     type = "auto",
                                     # New optimized parameters
                                     class_column = "Class",
                                     case_label = "Case",
                                     control_label = "Control", 
                                     feature_prefix = "hsa_",
                                     p_threshold = 0.05,
                                     fc_threshold = 1.0,
                                     use_smote = TRUE,
                                     save_intermediate = TRUE,
                                     verbose = TRUE) {
  
  # Record original working directory
  original_wd <- getwd()
  
  # Setup cleanup on exit
  on.exit({
    if (getwd() != original_wd) {
      setwd(original_wd)
    }
  })
  
  # Change to working directory
  if (wd != getwd()) {
    setwd(wd)
  }
  
  # Create output directory for results
  output_dir <- file.path(wd, "omicselector_results")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Setup parallel processing parameters
  parallel_cores <- if (register_parallel) {
    if (!is.null(clx)) {
      # Use provided cluster
      NULL  # Let optimized function auto-detect
    } else {
      # Auto-detect optimal cores
      min(parallel::detectCores() - 1, 6)
    }
  } else {
    1  # Sequential processing
  }
  
  # Convert timeout from seconds to minutes
  timeout_minutes <- ceiling(timeout_sec / 60)
  
  if (verbose) {
    cat("=============================================================\n")
    cat("OmicSelector - Optimized miRNA Feature Selection Analysis\n")
    cat("=============================================================\n")
    cat("Working directory:", wd, "\n")
    cat("Methods to run:", length(m), "(", paste(head(m, 10), collapse = ", "), 
        if(length(m) > 10) "..." else "", ")\n")
    cat("Max features per method:", prefer_no_features, "\n")
    cat("Max iterations:", max_iterations, "\n")
    cat("Timeout per method:", timeout_minutes, "minutes\n")
    cat("Class column:", class_column, "\n")
    cat("Case label:", case_label, "\n")
    cat("Control label:", control_label, "\n") 
    cat("Feature prefix:", feature_prefix, "\n")
    cat("P-value threshold:", p_threshold, "\n")
    cat("Fold change threshold:", fc_threshold, "\n")
    cat("Use SMOTE:", use_smote, "\n")
    cat("Parallel processing:", register_parallel, "\n")
    cat("=============================================================\n\n")
  }
  
  # Call the optimized implementation
  tryCatch({
    results <- omicselector_optimized(
      data_path = wd,
      methods = m,
      class_column = class_column,
      case_label = case_label,
      control_label = control_label,
      feature_prefix = feature_prefix,
      max_features = prefer_no_features,
      max_iterations = max_iterations,
      p_threshold = p_threshold,
      fc_threshold = fc_threshold,
      parallel_cores = parallel_cores,
      use_smote = use_smote,
      timeout_minutes = timeout_minutes,
      save_intermediate = save_intermediate,
      output_dir = output_dir,
      verbose = verbose,
      seed = if (is.numeric(stamp)) stamp else 42
    )
    
    # Safety check: ensure results is a list
    if (!is.list(results)) {
      stop("Results object is not a list (it's a ", class(results)[1], ")")
    }
    
    # Create backward-compatible formula list
    formulas <- convert_to_legacy_formulas(results$selected_features, class_column)
    
    # Add legacy compatibility
    results$formulas <- formulas
    results$stamp <- stamp
    results$working_directory <- wd
    
    # Add legacy class for compatibility
    class(results) <- c("omicselector_results", "omicselector_legacy", "list")
    
    if (verbose) {
      cat("\n=============================================================\n")
      cat("Analysis completed successfully!\n")
      # Safe access with defaults
      runtime <- if (is.list(results$metadata) && !is.null(results$metadata$total_runtime)) {
        round(results$metadata$total_runtime, 2)
      } else 0
      
      successful <- if (is.list(results$summary) && !is.null(results$summary$successful_methods)) {
        results$summary$successful_methods
      } else 0
      
      total <- if (is.list(results$summary) && !is.null(results$summary$total_methods)) {
        results$summary$total_methods
      } else 0
      
      cat("Runtime:", runtime, "minutes\n")
      cat("Successful methods:", successful, "/", total, "\n")
      cat("Results saved to:", output_dir, "\n")
      cat("=============================================================\n")
    }
    
    return(results)
    
  }, error = function(e) {
    if (verbose) {
      cat("\nERROR: Analysis failed with message:\n")
      cat(e$message, "\n")
      cat("\nFalling back to legacy implementation...\n")
    }
    
    # Fallback to legacy implementation if optimized version fails
    return(run_legacy_omicselector(
      wd = wd, m = m, max_iterations = max_iterations,
      stamp = stamp, prefer_no_features = prefer_no_features,
      debug = debug, timeout_sec = timeout_sec,
      register_parallel = register_parallel, clx = clx,
      verbose = verbose
    ))
  })
}

#' Convert Selected Features to Legacy Formula Format
#'
#' @param selected_features List of feature vectors by method
#' @param class_column Character name of class column
#' @return List of formulas compatible with legacy code
#' @keywords internal
convert_to_legacy_formulas <- function(selected_features, class_column = "Class") {
  
  # Safety check
  if (!is.list(selected_features)) {
    warning("selected_features is not a list, creating empty formulas")
    return(list())
  }
  
  formulas <- list()
  
  for (method_name in names(selected_features)) {
    features <- selected_features[[method_name]]
    
    if (length(features) > 0) {
      # Create formula string
      formula_str <- paste(class_column, "~", paste(paste0("`", features, "`"), collapse = " + "))
      formulas[[method_name]] <- as.formula(formula_str)
    } else {
      # Empty feature set
      formulas[[method_name]] <- as.formula(paste(class_column, "~ 1"))
    }
  }
  
  return(formulas)
}

#' Legacy Fallback Implementation
#'
#' @description
#' Simplified fallback version that maintains basic functionality
#' if the optimized version fails for any reason.
#'
#' @param ... Parameters from main function
#' @return Legacy-style results
#' @keywords internal
run_legacy_omicselector <- function(wd, m, max_iterations, stamp, prefer_no_features, 
                                   debug, timeout_sec, register_parallel, clx, verbose) {
  
  if (verbose) cat("Running legacy fallback implementation...\n")
  
  # Basic data loading
  tryCatch({
    train_data <- read.csv(file.path(wd, "mixed_train.csv"), stringsAsFactors = FALSE)
    
    # Ensure Class column is factor
    if ("Class" %in% colnames(train_data)) {
      train_data$Class <- factor(train_data$Class, levels = c("Control", "Case"))
    }
    
    # Get hsa_ features
    hsa_features <- colnames(train_data)[startsWith(colnames(train_data), "hsa_")]
    
    if (length(hsa_features) == 0) {
      stop("No hsa_ features found in training data")
    }
    
    # Simple feature selection methods
    formulas <- list()
    
    # Method 1: All features
    if (1 %in% m) {
      formulas[["all"]] <- as.formula(paste("Class ~", paste(paste0("`", hsa_features, "`"), collapse = " + ")))
    }
    
    # Method 2: Simple t-test based selection
    if (2 %in% m) {
      if ("Class" %in% colnames(train_data)) {
        # Simple p-value based selection
        pvalues <- numeric(length(hsa_features))
        names(pvalues) <- hsa_features
        
        for (i in seq_along(hsa_features)) {
          feature <- hsa_features[i]
          case_vals <- train_data[train_data$Class == "Case", feature]
          control_vals <- train_data[train_data$Class == "Control", feature]
          
          if (length(case_vals) > 1 && length(control_vals) > 1) {
            tryCatch({
              t_result <- t.test(case_vals, control_vals)
              pvalues[i] <- t_result$p.value
            }, error = function(e) {
              pvalues[i] <- 1.0
            })
          } else {
            pvalues[i] <- 1.0
          }
        }
        
        # Select top features by p-value
        n_top <- min(prefer_no_features, length(pvalues))
        if (n_top > 0) {
          top_features <- names(sort(pvalues))[seq_len(n_top)]
          formulas[["sig"]] <- as.formula(paste("Class ~", paste(paste0("`", top_features, "`"), collapse = " + ")))
        }
      }
    }
    
    # Create legacy-compatible results
    results <- list(
      formulas = formulas,
      stamp = stamp,
      working_directory = wd,
      metadata = list(
        start_time = Sys.time(),
        methods_requested = m,
        fallback_mode = TRUE
      )
    )
    
    class(results) <- c("omicselector_legacy", "list")
    
    if (verbose) cat("Legacy fallback completed with", length(formulas), "methods\n")
    
    return(results)
    
  }, error = function(e) {
    if (verbose) cat("Legacy fallback also failed:", e$message, "\n")
    
    # Ultimate fallback - return empty structure
    return(list(
      formulas = list(),
      stamp = stamp,
      error = e$message,
      metadata = list(
        start_time = Sys.time(),
        methods_requested = m,
        failed = TRUE
      )
    ))
  })
}

#' Print Method for Legacy OmicSelector Results
#'
#' @param x OmicSelector legacy results object
#' @param ... Additional arguments
#' @export
print.omicselector_legacy <- function(x, ...) {
  cat("OmicSelector Results (Legacy Compatible)\n")
  cat("======================================\n\n")
  
  if (!is.null(x$metadata$failed) && x$metadata$failed) {
    cat("Analysis failed:", x$error %||% "Unknown error", "\n")
    return(invisible(x))
  }
  
  if (!is.null(x$metadata$fallback_mode) && x$metadata$fallback_mode) {
    cat("Note: Results from fallback mode\n\n")
  }
  
  cat("Working directory:", x$working_directory %||% "Unknown", "\n")
  cat("Methods requested:", length(x$metadata$methods_requested), "\n")
  cat("Formulas created:", length(x$formulas), "\n")
  
  if (length(x$formulas) > 0) {
    cat("\nAvailable formulas:\n")
    for (name in names(x$formulas)) {
      formula_str <- deparse(x$formulas[[name]])
      # Truncate long formulas
      if (nchar(formula_str) > 60) {
        formula_str <- paste0(substr(formula_str, 1, 57), "...")
      }
      cat("  ", name, ":", formula_str, "\n")
    }
  }
  
  invisible(x)
}

#' Summary Method for Legacy OmicSelector Results
#'
#' @param object OmicSelector legacy results object
#' @param ... Additional arguments
#' @export
summary.omicselector_legacy <- function(object, ...) {
  print(object)
  invisible(object)
}
