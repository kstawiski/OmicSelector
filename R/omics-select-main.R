#' Modern OmicSelector Main Function
#'
#' @description
#' Refactored and modernized version of the main OmicSelector feature selection function.
#' This function orchestrates the entire feature selection workflow with improved
#' error handling, logging, and modular design.
#'
#' @param wd Character string. Working directory containing required data files
#'   (mixed_train.csv, mixed_test.csv, mixed_validation.csv)
#' @param methods Integer vector. Feature selection methods to apply (1-70)
#' @param config_name Character string. Configuration profile to use 
#'   ("default", "development", "hpc")
#' @param config_override List. Additional configuration parameters to override
#' @param stamp Character string. Unique identifier for this run
#' @param debug Logical. Enable debug mode with additional output
#'
#' @return OmicSelector results object (S3 class)
#'
#' @examples
#' \dontrun{
#' # Basic usage with default settings
#' results <- omics_select(
#'   wd = "path/to/data",
#'   methods = c(1, 2, 3, 11, 17)
#' )
#'
#' # Use development configuration for faster testing
#' results <- omics_select(
#'   wd = "path/to/data", 
#'   methods = c(1, 2, 3),
#'   config_name = "development"
#' )
#'
#' # Custom configuration override
#' results <- omics_select(
#'   wd = "path/to/data",
#'   methods = c(1:5),
#'   config_override = list(max_iterations = 5, cores = 2)
#' )
#' }
#'
#' @export
omics_select <- function(wd = getwd(),
                        methods = NULL,
                        config_name = "default",
                        config_override = list(),
                        stamp = NULL,
                        debug = FALSE) {
  
  # Initialize timing
  start_time <- Sys.time()
  
  # Generate unique stamp if not provided
  if (is.null(stamp)) {
    stamp <- format(start_time, "omics_%Y%m%d_%H%M%S")
  }
  
  # Load and validate configuration
  config <- get_config(config_name, config_file = NULL, !!!config_override)
  config <- validate_config(config)
  
  if (debug) {
    config$log_level <- "DEBUG"
  }
  
  # Setup logging
  log_dir <- file.path(wd, "temp")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  
  log_file <- if (config$log_to_file) file.path(log_dir, paste0(stamp, ".log")) else NULL
  setup_logging(config$log_level, log_file)
  
  log_info("Starting OmicSelector feature selection")
  log_info("Working directory: {wd}")
  log_info("Configuration: {config_name}")
  log_info("Stamp: {stamp}")
  
  tryCatch({
    # Validate inputs
    validate_working_directory(wd)
    
    if (is.null(methods)) {
      methods <- config$default_methods
      log_info("Using default methods: {paste(methods, collapse = ', ')}")
    }
    methods <- validate_methods(methods)
    
    # Setup parallel processing
    parallel_config <- validate_parallel_config(config$cores, config$register_parallel)
    setup_parallel_processing(parallel_config)
    
    # Load and validate data
    log_info("Loading data files")
    data_list <- load_omics_data(wd)
    
    # Initialize results object
    results <- initialize_results_object(stamp, config, methods, data_list)
    
    # Run feature selection methods
    log_info("Running {length(methods)} feature selection methods")
    results <- run_feature_selection_pipeline(results, methods, config)
    
    # Finalize results
    end_time <- Sys.time()
    results$runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))
    results$completed_at <- end_time
    
    log_info("Feature selection completed in {round(results$runtime, 2)} minutes")
    
    # Save results
    save_results(results, wd)
    
    return(results)
    
  }, error = function(e) {
    log_error("Feature selection failed: {e$message}")
    stop("OmicSelector failed: ", e$message, call. = FALSE)
  }, finally = {
    cleanup_parallel_processing()
  })
}

#' Load Omics Data Files
#'
#' @param wd Working directory path
#' @return List containing loaded data frames
#' @keywords internal
load_omics_data <- function(wd) {
  
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(wd)
  
  data_files <- c(
    "train" = "mixed_train.csv",
    "test" = "mixed_test.csv",
    "validation" = "mixed_validation.csv"
  )
  
  data_list <- list()
  
  for (name in names(data_files)) {
    file_path <- data_files[name]
    
    if (file.exists(file_path)) {
      log_debug("Loading {name} data from {file_path}")
      data_list[[name]] <- safe_read_csv(file_path)
      
      if (is.null(data_list[[name]])) {
        stop("Failed to load ", file_path, call. = FALSE)
      }
      
      validate_data_frame(data_list[[name]], min_rows = 5, min_cols = 2)
      log_info("Loaded {name} data: {nrow(data_list[[name]])} rows, {ncol(data_list[[name]])} columns")
      
    } else if (name %in% c("train", "test")) {
      stop("Required file missing: ", file_path, call. = FALSE)
    } else {
      log_warn("Optional file missing: {file_path}")
    }
  }
  
  return(data_list)
}

#' Initialize Results Object
#'
#' @param stamp Character string with run identifier
#' @param config Configuration list
#' @param methods Integer vector with methods
#' @param data_list List with loaded data
#' @return OmicSelector results object
#' @keywords internal
initialize_results_object <- function(stamp, config, methods, data_list) {
  
  results <- list(
    stamp = stamp,
    config = config,
    methods = methods,
    data_info = list(
      train_samples = nrow(data_list$train),
      train_features = ncol(data_list$train) - 1,  # excluding class column
      test_samples = if (!is.null(data_list$test)) nrow(data_list$test) else 0,
      validation_samples = if (!is.null(data_list$validation)) nrow(data_list$validation) else 0
    ),
    formulas = list(),
    method_results = list(),
    errors = list(),
    warnings = list(),
    started_at = Sys.time(),
    completed_at = NULL,
    runtime = NULL
  )
  
  class(results) <- c("OmicSelector", "list")
  
  return(results)
}

#' Setup Parallel Processing
#'
#' @param parallel_config Parallel configuration list
#' @keywords internal
setup_parallel_processing <- function(parallel_config) {
  
  if (parallel_config$register_parallel && parallel_config$cores > 1) {
    
    log_info("Setting up parallel processing with {parallel_config$cores} cores")
    
    if (requireNamespace("future", quietly = TRUE) && requireNamespace("furrr", quietly = TRUE)) {
      # Use modern future/furrr approach
      future::plan(future::multisession, workers = parallel_config$cores)
      log_debug("Using future/furrr parallel backend")
      
    } else if (requireNamespace("doParallel", quietly = TRUE)) {
      # Fallback to doParallel
      .omics_cluster <- parallel::makeCluster(parallel_config$cores)
      doParallel::registerDoParallel(.omics_cluster)
      # Store cluster in global environment for cleanup
      assign(".omics_cluster", .omics_cluster, envir = .GlobalEnv)
      log_debug("Using doParallel backend")
      
    } else {
      log_warn("No parallel backend available, running sequentially")
    }
  } else {
    log_info("Running in sequential mode")
  }
}

#' Cleanup Parallel Processing
#'
#' @keywords internal
cleanup_parallel_processing <- function() {
  
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
  }
  
  # Clean up doParallel cluster if it exists
  if (exists(".omics_cluster", envir = .GlobalEnv)) {
    cluster <- get(".omics_cluster", envir = .GlobalEnv)
    try(parallel::stopCluster(cluster), silent = TRUE)
    rm(".omics_cluster", envir = .GlobalEnv)
  }
  
  log_debug("Parallel processing cleaned up")
}

#' Save Results to Files
#'
#' @param results OmicSelector results object
#' @param wd Working directory
#' @keywords internal
save_results <- function(results, wd) {
  
  output_dir <- file.path(wd, "temp")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Save complete results object
  results_file <- file.path(output_dir, paste0(results$stamp, "_results.rds"))
  safe_save_rds(results_file, results)
  
  # Save formulas for backward compatibility
  if (length(results$formulas) > 0) {
    formulas_file <- file.path(wd, "featureselection_formulas_final.RDS")
    safe_save_rds(formulas_file, results$formulas)
  }
  
  # Save summary report
  create_summary_report(results, wd)
  
  log_info("Results saved to {output_dir}")
}

#' Create Summary Report
#'
#' @param results OmicSelector results object
#' @param wd Working directory
#' @keywords internal
create_summary_report <- function(results, wd) {
  
  report_file <- file.path(wd, "temp", paste0(results$stamp, "_summary.txt"))
  
  tryCatch({
    sink(report_file)
    
    cat("OmicSelector Feature Selection Summary\n")
    cat("=====================================\n\n")
    cat("Run ID:", results$stamp, "\n")
    cat("Started:", format(results$started_at), "\n")
    cat("Completed:", format(results$completed_at), "\n")
    cat("Runtime:", round(results$runtime, 2), "minutes\n\n")
    
    cat("Data Information:\n")
    cat("- Training samples:", results$data_info$train_samples, "\n")
    cat("- Training features:", results$data_info$train_features, "\n")
    cat("- Test samples:", results$data_info$test_samples, "\n")
    cat("- Validation samples:", results$data_info$validation_samples, "\n\n")
    
    cat("Methods Requested:", length(results$methods), "\n")
    cat("Methods Completed:", length(results$formulas), "\n")
    cat("Methods Failed:", length(results$errors), "\n\n")
    
    if (length(results$formulas) > 0) {
      cat("Successfully Generated Formulas:\n")
      for (i in seq_along(results$formulas)) {
        cat("  ", i, ":", results$formulas[[i]], "\n")
      }
    }
    
    if (length(results$errors) > 0) {
      cat("\nErrors Encountered:\n")
      for (i in seq_along(results$errors)) {
        cat("  Method", names(results$errors)[i], ":", results$errors[[i]], "\n")
      }
    }
    
    sink()
    log_info("Summary report saved to {report_file}")
    
  }, error = function(e) {
    if (sink.number() > 0) sink()
    log_warn("Failed to create summary report: {e$message}")
  })
}
