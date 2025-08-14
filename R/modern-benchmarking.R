#' Modern Benchmarking System for OmicSelector
#'
#' @description
#' Modernized benchmarking function with improved parallel processing,
#' better error handling, and integration with tidymodels ecosystem.
#'
#' @param wd Working directory containing data files
#' @param formulas List of formulas to benchmark (from feature selection)
#' @param algorithms Character vector of machine learning algorithms to test
#' @param validation_strategy Character: "holdout", "cv", or "bootstrap"
#' @param cv_folds Integer: number of CV folds (default 10)
#' @param cv_repeats Integer: number of CV repeats (default 5)
#' @param search_iterations Integer: hyperparameter search iterations
#' @param parallel_processing Logical: use parallel processing
#' @param output_file Character: output CSV file name
#' @param config Configuration object (optional)
#'
#' @return Benchmarking results object
#'
#' @examples
#' \dontrun{
#' # Basic benchmarking
#' results <- omics_benchmark(
#'   wd = "path/to/data",
#'   formulas = list("class ~ feature1 + feature2"),
#'   algorithms = c("rf", "svm")
#' )
#' 
#' # Advanced benchmarking with tidymodels
#' results <- omics_benchmark(
#'   wd = "path/to/data", 
#'   formulas = feature_selection_results$formulas,
#'   algorithms = c("random_forest", "svm_radial", "logistic_reg"),
#'   validation_strategy = "cv",
#'   cv_folds = 10
#' )
#' }
#'
#' @export
omics_benchmark <- function(wd = getwd(),
                           formulas = NULL,
                           algorithms = c("random_forest", "svm_radial", "logistic_reg"),
                           validation_strategy = "holdout",
                           cv_folds = 10,
                           cv_repeats = 5,
                           search_iterations = 100,
                           parallel_processing = TRUE,
                           output_file = "benchmark_results.csv",
                           config = NULL) {
  
  start_time <- Sys.time()
  
  log_info("Starting OmicSelector benchmarking")
  log_info("Working directory: {wd}")
  log_info("Algorithms: {paste(algorithms, collapse = ', ')}")
  log_info("Validation strategy: {validation_strategy}")
  
  # Load configuration
  if (is.null(config)) {
    config <- get_config()
  }
  
  # Setup logging
  stamp <- format(start_time, "benchmark_%Y%m%d_%H%M%S")
  log_dir <- file.path(wd, "temp")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  
  # Load formulas if not provided
  if (is.null(formulas)) {
    formulas_file <- file.path(wd, "featureselection_formulas_final.RDS")
    if (file.exists(formulas_file)) {
      formulas <- readRDS(formulas_file)
      log_info("Loaded {length(formulas)} formulas from {formulas_file}")
    } else {
      stop("No formulas provided and no formula file found", call. = FALSE)
    }
  }
  
  # Validate inputs
  validate_working_directory(wd, required_files = c("mixed_train.csv", "mixed_test.csv"))
  
  tryCatch({
    # Load and prepare data
    data_list <- load_benchmark_data(wd)
    
    # Setup parallel processing
    if (parallel_processing) {
      parallel_config <- validate_parallel_config(config$cores, TRUE)
      setup_parallel_processing(parallel_config)
    }
    
    # Initialize results object
    benchmark_results <- initialize_benchmark_results(
      stamp = stamp,
      formulas = formulas,
      algorithms = algorithms,
      config = config,
      data_info = get_data_info(data_list)
    )
    
    # Run benchmarking
    log_info("Running benchmark for {length(formulas)} formulas with {length(algorithms)} algorithms")
    benchmark_results <- run_benchmark_pipeline(
      benchmark_results, 
      data_list, 
      validation_strategy,
      cv_folds,
      cv_repeats,
      search_iterations
    )
    
    # Finalize and save results
    end_time <- Sys.time()
    benchmark_results$runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))
    benchmark_results$completed_at <- end_time
    
    save_benchmark_results(benchmark_results, wd, output_file)
    
    log_info("Benchmarking completed in {round(benchmark_results$runtime, 2)} minutes")
    
    return(benchmark_results)
    
  }, error = function(e) {
    log_error("Benchmarking failed: {e$message}")
    stop("Benchmarking failed: ", e$message, call. = FALSE)
    
  }, finally = {
    if (parallel_processing) {
      cleanup_parallel_processing()
    }
  })
}

#' Load Benchmark Data
#'
#' @param wd Working directory
#' @return List with train, test, and optional validation data
#' @keywords internal
load_benchmark_data <- function(wd) {
  
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(wd)
  
  # Load training data
  train_data <- safe_read_csv("mixed_train.csv")
  if (is.null(train_data)) {
    stop("Failed to load training data", call. = FALSE)
  }
  
  # Load test data
  test_data <- safe_read_csv("mixed_test.csv")
  if (is.null(test_data)) {
    stop("Failed to load test data", call. = FALSE)
  }
  
  # Load validation data (optional)
  validation_data <- NULL
  if (file.exists("mixed_validation.csv")) {
    validation_data <- safe_read_csv("mixed_validation.csv")
    if (!is.null(validation_data)) {
      log_info("Loaded validation set: {nrow(validation_data)} samples")
    }
  }
  
  log_info("Loaded benchmark data - Train: {nrow(train_data)}, Test: {nrow(test_data)}")
  
  return(list(
    train = train_data,
    test = test_data,
    validation = validation_data
  ))
}

#' Initialize Benchmark Results Object
#'
#' @param stamp Run identifier
#' @param formulas List of formulas
#' @param algorithms Algorithm names
#' @param config Configuration
#' @param data_info Data information
#' @return Benchmark results object
#' @keywords internal
initialize_benchmark_results <- function(stamp, formulas, algorithms, config, data_info) {
  
  results <- list(
    stamp = stamp,
    formulas = formulas,
    algorithms = algorithms,
    config = config,
    data_info = data_info,
    results_matrix = NULL,
    detailed_results = list(),
    model_performances = list(),
    best_models = list(),
    errors = list(),
    started_at = Sys.time(),
    completed_at = NULL,
    runtime = NULL
  )
  
  class(results) <- c("OmicBenchmark", "list")
  
  return(results)
}

#' Get Data Information
#'
#' @param data_list List with data frames
#' @return Data information list
#' @keywords internal
get_data_info <- function(data_list) {
  
  info <- list(
    train_samples = nrow(data_list$train),
    train_features = ncol(data_list$train) - 1,
    test_samples = nrow(data_list$test),
    validation_samples = if (!is.null(data_list$validation)) nrow(data_list$validation) else 0,
    class_distribution_train = table(data_list$train[, 1]),
    class_distribution_test = table(data_list$test[, 1])
  )
  
  return(info)
}

#' Run Benchmark Pipeline
#'
#' @param benchmark_results Benchmark results object
#' @param data_list Data list
#' @param validation_strategy Validation strategy
#' @param cv_folds CV folds
#' @param cv_repeats CV repeats
#' @param search_iterations Search iterations
#' @return Updated benchmark results
#' @keywords internal
run_benchmark_pipeline <- function(benchmark_results, 
                                 data_list,
                                 validation_strategy,
                                 cv_folds,
                                 cv_repeats,
                                 search_iterations) {
  
  formulas <- benchmark_results$formulas
  algorithms <- benchmark_results$algorithms
  
  # Initialize results matrix
  n_combinations <- length(formulas) * length(algorithms)
  results_matrix <- data.frame(
    formula_id = integer(n_combinations),
    algorithm = character(n_combinations),
    accuracy = numeric(n_combinations),
    auc = numeric(n_combinations),
    sensitivity = numeric(n_combinations),
    specificity = numeric(n_combinations),
    precision = numeric(n_combinations),
    f1_score = numeric(n_combinations),
    balanced_accuracy = numeric(n_combinations),
    stringsAsFactors = FALSE
  )
  
  # Progress tracking
  pb <- create_progress_bar(n_combinations, "Running benchmark combinations")
  combination_index <- 0
  
  # Main benchmarking loop
  for (formula_idx in seq_along(formulas)) {
    formula <- formulas[[formula_idx]]
    
    for (algo_idx in seq_along(algorithms)) {
      algorithm <- algorithms[algo_idx]
      combination_index <- combination_index + 1
      
      log_debug("Testing formula {formula_idx} with algorithm {algorithm}")
      
      tryCatch({
        # Run single benchmark combination
        perf_metrics <- run_single_benchmark(
          formula = formula,
          algorithm = algorithm,
          data_list = data_list,
          validation_strategy = validation_strategy,
          cv_folds = cv_folds,
          cv_repeats = cv_repeats,
          search_iterations = search_iterations
        )
        
        # Store results
        results_matrix[combination_index, ] <- list(
          formula_id = formula_idx,
          algorithm = algorithm,
          accuracy = perf_metrics$accuracy,
          auc = perf_metrics$auc,
          sensitivity = perf_metrics$sensitivity,
          specificity = perf_metrics$specificity,
          precision = perf_metrics$precision,
          f1_score = perf_metrics$f1_score,
          balanced_accuracy = perf_metrics$balanced_accuracy
        )
        
        # Store detailed results
        key <- paste(formula_idx, algorithm, sep = "_")
        benchmark_results$detailed_results[[key]] <- perf_metrics
        
      }, error = function(e) {
        error_msg <- paste("Formula", formula_idx, "with", algorithm, "failed:", e$message)
        benchmark_results$errors[[paste(formula_idx, algorithm, sep = "_")]] <- error_msg
        log_error(error_msg)
        
        # Fill with NA values
        results_matrix[combination_index, ] <- list(
          formula_id = formula_idx,
          algorithm = algorithm,
          accuracy = NA,
          auc = NA,
          sensitivity = NA,
          specificity = NA,
          precision = NA,
          f1_score = NA,
          balanced_accuracy = NA
        )
      })
      
      update_progress_bar(pb, combination_index)
    }
  }
  
  finish_progress_bar(pb)
  
  # Store results matrix
  benchmark_results$results_matrix <- results_matrix
  
  # Find best models
  benchmark_results$best_models <- find_best_models(results_matrix)
  
  log_info("Benchmark pipeline completed")
  return(benchmark_results)
}

#' Run Single Benchmark Combination
#'
#' @param formula Model formula
#' @param algorithm Algorithm name
#' @param data_list Data list
#' @param validation_strategy Validation strategy
#' @param cv_folds CV folds
#' @param cv_repeats CV repeats
#' @param search_iterations Search iterations
#' @return Performance metrics
#' @keywords internal
run_single_benchmark <- function(formula,
                                algorithm,
                                data_list,
                                validation_strategy,
                                cv_folds,
                                cv_repeats,
                                search_iterations) {
  
  # Prepare data based on formula
  train_data <- prepare_model_data(data_list$train, formula)
  test_data <- prepare_model_data(data_list$test, formula)
  
  # Choose modeling approach
  if (requireNamespace("tidymodels", quietly = TRUE) && 
      algorithm %in% c("random_forest", "svm_radial", "logistic_reg")) {
    
    # Use tidymodels approach
    return(run_tidymodels_benchmark(
      train_data, test_data, algorithm, 
      validation_strategy, cv_folds, cv_repeats, search_iterations
    ))
    
  } else {
    # Fallback to caret approach
    return(run_caret_benchmark(
      train_data, test_data, algorithm,
      validation_strategy, cv_folds, cv_repeats
    ))
  }
}

#' Save Benchmark Results
#'
#' @param benchmark_results Benchmark results object
#' @param wd Working directory
#' @param output_file Output file name
#' @keywords internal
save_benchmark_results <- function(benchmark_results, wd, output_file) {
  
  # Save complete results object
  results_file <- file.path(wd, "temp", paste0(benchmark_results$stamp, "_benchmark_results.rds"))
  safe_save_rds(results_file, benchmark_results)
  
  # Save results matrix as CSV
  csv_file <- file.path(wd, output_file)
  write.csv(benchmark_results$results_matrix, csv_file, row.names = FALSE)
  
  # Create summary report
  create_benchmark_summary(benchmark_results, wd)
  
  log_info("Benchmark results saved to {wd}")
}

#' Find Best Models
#'
#' @param results_matrix Results matrix
#' @return List of best models by different metrics
#' @keywords internal
find_best_models <- function(results_matrix) {
  
  # Remove rows with all NA values
  complete_results <- results_matrix[complete.cases(results_matrix[, c("accuracy", "auc")]), ]
  
  if (nrow(complete_results) == 0) {
    return(list())
  }
  
  best_models <- list(
    best_accuracy = complete_results[which.max(complete_results$accuracy), ],
    best_auc = complete_results[which.max(complete_results$auc), ],
    best_f1 = complete_results[which.max(complete_results$f1_score), ],
    best_balanced_accuracy = complete_results[which.max(complete_results$balanced_accuracy), ]
  )
  
  return(best_models)
}

# Additional helper functions for tidymodels and caret implementations would go here
# These are simplified versions - full implementation would include complete 
# hyperparameter tuning and cross-validation logic
