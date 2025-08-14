#' Modular Feature Selection Framework
#'
#' @description
#' A plugin-based, modular framework for feature selection methods that allows
#' easy addition of new methods without modifying core code. Each method is
#' self-contained and auto-registers with the system.
#'
#' @author OmicSelector Team
#' @export

# =============================================================================
# FEATURE SELECTION METHOD REGISTRY
# =============================================================================

# Global registry for feature selection methods
.fs_method_registry <- new.env(parent = emptyenv())

#' Register a Feature Selection Method
#'
#' @description
#' Register a new feature selection method with the framework. This function
#' is typically called automatically when method modules are loaded.
#'
#' @param method_id Unique method identifier (integer)
#' @param method_name Human-readable method name
#' @param method_function Function that implements the method
#' @param category Method category (statistical, wrapper, embedded, etc.)
#' @param dependencies Required packages (optional)
#' @param description Method description
#' @param parameters Default parameters for the method
#' @param complexity Computational complexity (low, medium, high)
#' @param supports_smote Whether method supports SMOTE data
#' @param timeout_default Default timeout in seconds
#'
#' @return Invisible TRUE on success
#' @export
register_fs_method <- function(method_id,
                              method_name,
                              method_function,
                              category = "statistical",
                              dependencies = character(0),
                              description = "",
                              parameters = list(),
                              complexity = "medium",
                              supports_smote = TRUE,
                              timeout_default = 300) {
  
  # Validate inputs
  if (!is.numeric(method_id) || method_id <= 0) {
    stop("method_id must be a positive integer")
  }
  
  if (!is.function(method_function)) {
    stop("method_function must be a function")
  }
  
  # Create method definition
  method_def <- list(
    id = method_id,
    name = method_name,
    function_obj = method_function,
    category = category,
    dependencies = dependencies,
    description = description,
    parameters = parameters,
    complexity = complexity,
    supports_smote = supports_smote,
    timeout_default = timeout_default,
    registered_at = Sys.time()
  )
  
  # Store in registry
  .fs_method_registry[[as.character(method_id)]] <- method_def
  
  invisible(TRUE)
}

#' Get Registered Feature Selection Method
#'
#' @param method_id Method ID to retrieve
#' @return Method definition list or NULL
#' @export
get_fs_method <- function(method_id) {
  .fs_method_registry[[as.character(method_id)]]
}

#' List All Registered Methods
#'
#' @param category Filter by category (optional)
#' @param complexity Filter by complexity (optional)
#' @return Data frame of registered methods
#' @export
list_fs_methods <- function(category = NULL, complexity = NULL) {
  
  if (length(.fs_method_registry) == 0) {
    return(data.frame())
  }
  
  # Extract method information
  methods_df <- data.frame(
    id = sapply(ls(.fs_method_registry), function(x) .fs_method_registry[[x]]$id),
    name = sapply(ls(.fs_method_registry), function(x) .fs_method_registry[[x]]$name),
    category = sapply(ls(.fs_method_registry), function(x) .fs_method_registry[[x]]$category),
    complexity = sapply(ls(.fs_method_registry), function(x) .fs_method_registry[[x]]$complexity),
    supports_smote = sapply(ls(.fs_method_registry), function(x) .fs_method_registry[[x]]$supports_smote),
    dependencies = sapply(ls(.fs_method_registry), function(x) {
      deps <- .fs_method_registry[[x]]$dependencies
      if (length(deps) == 0) "" else paste(deps, collapse = ", ")
    }),
    description = sapply(ls(.fs_method_registry), function(x) .fs_method_registry[[x]]$description),
    stringsAsFactors = FALSE
  )
  
  # Apply filters
  if (!is.null(category)) {
    methods_df <- methods_df[methods_df$category == category, ]
  }
  
  if (!is.null(complexity)) {
    methods_df <- methods_df[methods_df$complexity == complexity, ]
  }
  
  # Sort by ID
  methods_df <- methods_df[order(methods_df$id), ]
  rownames(methods_df) <- NULL
  
  return(methods_df)
}

#' Check Method Dependencies
#'
#' @param method_id Method ID to check
#' @return List with availability status and missing packages
#' @export
check_method_dependencies <- function(method_id) {
  
  method_def <- get_fs_method(method_id)
  if (is.null(method_def)) {
    return(list(available = FALSE, missing = character(0), error = "Method not found"))
  }
  
  dependencies <- method_def$dependencies
  if (length(dependencies) == 0) {
    return(list(available = TRUE, missing = character(0), error = NULL))
  }
  
  # Check each dependency
  missing <- character(0)
  for (pkg in dependencies) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }
  
  return(list(
    available = length(missing) == 0,
    missing = missing,
    error = if (length(missing) > 0) paste("Missing packages:", paste(missing, collapse = ", ")) else NULL
  ))
}

# =============================================================================
# STANDARDIZED METHOD INTERFACE
# =============================================================================

#' Standard Feature Selection Method Interface
#'
#' @description
#' All feature selection methods must implement this standardized interface.
#' This ensures consistency and allows the framework to work with any method.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list with parameters
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with components:
#'   - features: Character vector of selected feature names
#'   - scores: Optional numeric vector of feature scores/importance
#'   - metadata: List with method execution information
#'   - performance: Optional performance metrics
#'
#' @export
fs_method_interface <- function(data, config = list(), max_features = 20, use_smote = FALSE, timeout_sec = 300) {
  stop("This is an interface definition. Implement in specific methods.")
}

# =============================================================================
# MODULAR EXECUTION ENGINE
# =============================================================================

#' Execute Feature Selection Method
#'
#' @description
#' Execute a registered feature selection method with proper error handling,
#' timeout protection, and dependency checking.
#'
#' @param method_id Method ID to execute
#' @param data Training data
#' @param config Configuration parameters
#' @param max_features Maximum features to select
#' @param use_smote Apply SMOTE balancing
#' @param timeout_sec Timeout in seconds (uses method default if NULL)
#'
#' @return Method result list or NULL on failure
#' @export
execute_fs_method <- function(method_id, 
                             data, 
                             config = list(), 
                             max_features = 20, 
                             use_smote = FALSE, 
                             timeout_sec = NULL) {
  
  # Get method definition
  method_def <- get_fs_method(method_id)
  if (is.null(method_def)) {
    warning("Method ", method_id, " not found in registry")
    return(NULL)
  }
  
  # Check dependencies
  dep_check <- check_method_dependencies(method_id)
  if (!dep_check$available) {
    warning("Method ", method_id, " (", method_def$name, ") dependencies not available: ", dep_check$error)
    return(NULL)
  }
  
  # Use method default timeout if not specified
  if (is.null(timeout_sec)) {
    timeout_sec <- method_def$timeout_default
  }
  
  # Check SMOTE support
  if (use_smote && !method_def$supports_smote) {
    warning("Method ", method_id, " (", method_def$name, ") does not support SMOTE. Using original data.")
    use_smote <- FALSE
  }
  
  # Execute with error handling and timeout
  result <- tryCatch({
    
    # Add execution metadata
    execution_metadata <- list(
      method_id = method_id,
      method_name = method_def$name,
      execution_start = Sys.time(),
      config_used = config,
      max_features = max_features,
      use_smote = use_smote,
      timeout_sec = timeout_sec
    )
    
    # Execute method with timeout protection
    if (requireNamespace("future", quietly = TRUE)) {
      # Use future for timeout
      f <- future::future({
        method_def$function_obj(
          data = data,
          config = config,
          max_features = max_features,
          use_smote = use_smote,
          timeout_sec = timeout_sec
        )
      })
      
      method_result <- future::value(f, timeout = timeout_sec)
    } else {
      # Direct execution without timeout
      method_result <- method_def$function_obj(
        data = data,
        config = config,
        max_features = max_features,
        use_smote = use_smote,
        timeout_sec = timeout_sec
      )
    }
    
    # Add execution metadata to result
    execution_metadata$execution_end <- Sys.time()
    execution_metadata$execution_time <- as.numeric(
      execution_metadata$execution_end - execution_metadata$execution_start, 
      units = "secs"
    )
    
    if (is.list(method_result)) {
      method_result$metadata <- execution_metadata
    } else {
      method_result <- list(
        features = method_result,
        metadata = execution_metadata
      )
    }
    
    return(method_result)
    
  }, error = function(e) {
    warning("Method ", method_id, " (", method_def$name, ") failed: ", e$message)
    return(NULL)
  })
  
  return(result)
}

#' Execute Multiple Feature Selection Methods
#'
#' @description
#' Execute multiple feature selection methods in parallel or sequentially.
#'
#' @param method_ids Vector of method IDs to execute
#' @param data Training data
#' @param config Configuration parameters
#' @param max_features Maximum features per method
#' @param use_smote Apply SMOTE balancing
#' @param parallel Use parallel execution
#' @param timeout_sec Timeout per method
#' @param progress Show progress bar
#'
#' @return List of method results
#' @export
execute_multiple_fs_methods <- function(method_ids,
                                       data,
                                       config = list(),
                                       max_features = 20,
                                       use_smote = FALSE,
                                       parallel = TRUE,
                                       timeout_sec = NULL,
                                       progress = TRUE) {
  
  if (progress) {
    cat("Executing", length(method_ids), "feature selection methods...\n")
  }
  
  # Create progress bar
  if (progress && requireNamespace("progress", quietly = TRUE)) {
    pb <- progress::progress_bar$new(
      format = "Methods [:bar] :percent :elapsed",
      total = length(method_ids),
      width = 60
    )
  } else {
    pb <- NULL
  }
  
  # Execute methods
  results <- list()
  
  if (parallel && requireNamespace("future", quietly = TRUE) && requireNamespace("furrr", quietly = TRUE)) {
    # Parallel execution
    results <- furrr::future_map(method_ids, function(method_id) {
      if (!is.null(pb)) pb$tick()
      execute_fs_method(
        method_id = method_id,
        data = data,
        config = config,
        max_features = max_features,
        use_smote = use_smote,
        timeout_sec = timeout_sec
      )
    }, .names = as.character(method_ids))
  } else {
    # Sequential execution
    for (method_id in method_ids) {
      if (!is.null(pb)) pb$tick()
      
      result <- execute_fs_method(
        method_id = method_id,
        data = data,
        config = config,
        max_features = max_features,
        use_smote = use_smote,
        timeout_sec = timeout_sec
      )
      
      results[[as.character(method_id)]] <- result
    }
  }
  
  if (progress) {
    successful <- sum(sapply(results, function(x) !is.null(x)))
    cat("\nCompleted:", successful, "/", length(method_ids), "methods successful\n")
  }
  
  return(results)
}

# =============================================================================
# METHOD DISCOVERY AND LOADING
# =============================================================================

#' Load Feature Selection Method Modules
#'
#' @description
#' Automatically discover and load feature selection method modules from
#' the fs_methods directory or package.
#'
#' @param methods_dir Directory containing method modules
#' @param pattern File pattern for method modules
#' @param verbose Show loading messages
#'
#' @return Number of methods loaded
#' @export
load_fs_method_modules <- function(methods_dir = "R/fs_methods", 
                                  pattern = "^fs_method_.*\\.R$", 
                                  verbose = TRUE) {
  
  initial_count <- length(ls(.fs_method_registry))
  
  # Check if methods directory exists
  if (!dir.exists(methods_dir)) {
    if (verbose) {
      cat("Methods directory", methods_dir, "not found. Creating it...\n")
    }
    dir.create(methods_dir, recursive = TRUE)
    return(0)
  }
  
  # Find method files
  method_files <- list.files(methods_dir, pattern = pattern, full.names = TRUE)
  
  if (length(method_files) == 0) {
    if (verbose) {
      cat("No method modules found in", methods_dir, "\n")
    }
    return(0)
  }
  
  if (verbose) {
    cat("Loading", length(method_files), "method modules...\n")
  }
  
  # Load each method file
  for (file_path in method_files) {
    tryCatch({
      source(file_path)
      if (verbose) {
        cat("  ✓ Loaded:", basename(file_path), "\n")
      }
    }, error = function(e) {
      if (verbose) {
        cat("  ✗ Failed to load:", basename(file_path), "-", e$message, "\n")
      }
    })
  }
  
  final_count <- length(ls(.fs_method_registry))
  methods_loaded <- final_count - initial_count
  
  if (verbose) {
    cat("Loaded", methods_loaded, "new methods. Total methods:", final_count, "\n")
  }
  
  return(methods_loaded)
}

#' Initialize Feature Selection Framework
#'
#' @description
#' Initialize the modular feature selection framework by loading all
#' available method modules and core methods.
#'
#' @param methods_dir Directory containing method modules
#' @param load_core_methods Load built-in core methods
#' @param verbose Show initialization messages
#'
#' @return List with initialization status
#' @export
initialize_fs_framework <- function(methods_dir = "R/fs_methods",
                                   load_core_methods = TRUE,
                                   verbose = TRUE) {
  
  if (verbose) {
    cat("Initializing Modular Feature Selection Framework...\n")
    cat("================================================\n")
  }
  
  # Clear registry
  rm(list = ls(.fs_method_registry), envir = .fs_method_registry)
  
  # Load core methods if requested
  core_methods_loaded <- 0
  if (load_core_methods) {
    if (verbose) cat("Loading core methods...\n")
    core_methods_loaded <- load_core_fs_methods(verbose = verbose)
  }
  
  # Load method modules
  plugin_methods_loaded <- load_fs_method_modules(methods_dir = methods_dir, verbose = verbose)
  
  # Summary
  total_methods <- length(ls(.fs_method_registry))
  
  if (verbose) {
    cat("\n")
    cat("Framework Initialization Complete!\n")
    cat("==================================\n")
    cat("Core methods loaded:", core_methods_loaded, "\n")
    cat("Plugin methods loaded:", plugin_methods_loaded, "\n")
    cat("Total methods available:", total_methods, "\n")
    
    if (total_methods > 0) {
      cat("\nAvailable method categories:\n")
      methods_df <- list_fs_methods()
      categories <- table(methods_df$category)
      for (cat_name in names(categories)) {
        cat("  -", cat_name, ":", categories[cat_name], "methods\n")
      }
    }
  }
  
  return(list(
    success = TRUE,
    core_methods_loaded = core_methods_loaded,
    plugin_methods_loaded = plugin_methods_loaded,
    total_methods = total_methods
  ))
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Get Method Statistics
#'
#' @return List with registry statistics
#' @export
get_fs_framework_stats <- function() {
  
  methods_df <- list_fs_methods()
  
  if (nrow(methods_df) == 0) {
    return(list(
      total_methods = 0,
      categories = character(0),
      complexity_distribution = character(0),
      dependency_analysis = list()
    ))
  }
  
  # Category analysis
  categories <- table(methods_df$category)
  
  # Complexity analysis
  complexity <- table(methods_df$complexity)
  
  # Dependency analysis
  all_deps <- unlist(strsplit(methods_df$dependencies[methods_df$dependencies != ""], ", "))
  dep_freq <- if (length(all_deps) > 0) table(all_deps) else table(character(0))
  
  return(list(
    total_methods = nrow(methods_df),
    categories = categories,
    complexity_distribution = complexity,
    dependency_analysis = list(
      unique_packages = length(dep_freq),
      most_common = if (length(dep_freq) > 0) names(sort(dep_freq, decreasing = TRUE))[seq_len(min(5, length(dep_freq)))] else character(0)
    )
  ))
}
