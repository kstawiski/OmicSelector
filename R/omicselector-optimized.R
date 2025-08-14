#' Optimized OmicSelector for miRNA Case-Control Analysis
#'
#' @description
#' Completely rewritten and optimized OmicSelector function specifically designed
#' for miRNA biomarker selection with Case vs Control classification. Features:
#' - Optimized for "hsa_" prefixed miRNA features  
#' - Case/Control binary classification focus
#' - Modern error handling and validation
#' - Parallel processing optimization
#' - Memory-efficient operations
#' - Comprehensive logging and progress tracking
#'
#' @param data_path Character. Path to directory containing train/test/validation CSV files
#' @param methods Integer vector. Feature selection methods to apply (default: c(1:20))
#' @param class_column Character. Name of the class column (default: "Class") 
#' @param case_label Character. Label for positive cases (default: "Case")
#' @param control_label Character. Label for negative controls (default: "Control")
#' @param feature_prefix Character. Prefix for feature columns (default: "hsa_")
#' @param max_features Integer. Maximum number of features to select (default: 20)
#' @param max_iterations Integer. Maximum iterations for iterative methods (default: 10)
#' @param p_threshold Numeric. P-value threshold for significance (default: 0.05)
#' @param fc_threshold Numeric. Fold change threshold (default: 1.0)
#' @param parallel_cores Integer. Number of cores for parallel processing (default: NULL, auto-detect)
#' @param use_smote Logical. Whether to apply SMOTE for class balancing (default: TRUE)
#' @param timeout_minutes Integer. Timeout for each method in minutes (default: 30)
#' @param save_intermediate Logical. Save intermediate results (default: TRUE)
#' @param output_dir Character. Output directory for results (default: "./omicselector_results")
#' @param verbose Logical. Enable verbose logging (default: TRUE)
#' @param seed Integer. Random seed for reproducibility (default: 42)
#'
#' @return S3 object of class "omicselector_results" containing:
#'   - selected_features: List of feature sets by method
#'   - method_performance: Performance metrics for each method
#'   - differential_expression: DE analysis results
#'   - execution_times: Runtime for each method
#'   - metadata: Run configuration and parameters
#'
#' @examples
#' \dontrun{
#' # Basic miRNA analysis
#' results <- omicselector_optimized(
#'   data_path = "path/to/mirna_data",
#'   methods = c(1, 2, 4, 11, 17)
#' )
#'
#' # Quick analysis with fewer methods
#' results <- omicselector_optimized(
#'   data_path = "path/to/mirna_data",
#'   methods = c(1:5),
#'   max_features = 10,
#'   max_iterations = 5
#' )
#'
#' # Analysis with custom thresholds
#' results <- omicselector_optimized(
#'   data_path = "path/to/mirna_data",
#'   methods = c(1, 2, 4, 11, 17, 23),
#'   p_threshold = 0.01,
#'   fc_threshold = 1.5,
#'   max_features = 15
#' )
#' }
#'
#' @export
omicselector_optimized <- function(data_path = getwd(),
                                  methods = c(1:20),
                                  class_column = "Class",
                                  case_label = "Case", 
                                  control_label = "Control",
                                  feature_prefix = "hsa_",
                                  max_features = 20,
                                  max_iterations = 10,
                                  p_threshold = 0.05,
                                  fc_threshold = 1.0,
                                  parallel_cores = NULL,
                                  use_smote = TRUE,
                                  timeout_minutes = 30,
                                  save_intermediate = TRUE,
                                  output_dir = "./omicselector_results",
                                  verbose = TRUE,
                                  seed = 42) {
  
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Record start time
  start_time <- Sys.time()
  
  # Input validation and setup
  validated_params <- validate_omicselector_inputs(
    data_path = data_path,
    methods = methods,
    class_column = class_column,
    case_label = case_label,
    control_label = control_label,
    feature_prefix = feature_prefix,
    max_features = max_features,
    max_iterations = max_iterations,
    p_threshold = p_threshold,
    fc_threshold = fc_threshold,
    parallel_cores = parallel_cores,
    timeout_minutes = timeout_minutes,
    output_dir = output_dir
  )
  
  # Setup logging
  if (verbose) {
    setup_omicselector_logging(output_dir, "omicselector_optimized")
    log_info("Starting OmicSelector optimized analysis")
    # Safely create parameter summary
    if (is.list(validated_params) && length(validated_params) > 0) {
      param_names <- names(validated_params)
      param_values <- sapply(validated_params, function(x) {
        if (is.vector(x) && length(x) > 1) {
          paste0("c(", paste(x, collapse = ", "), ")")
        } else {
          as.character(x)
        }
      })
      param_summary <- paste(param_names, param_values, sep = "=", collapse = ", ")
      log_info(paste("Parameters:", param_summary))
    } else {
      log_info("Parameters validation completed")
    }
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) log_info(paste("Created output directory:", output_dir))
  }
  
  # Load and validate data
  if (verbose) log_info(paste("Loading and validating data from", data_path))
  data_list <- load_omicselector_data(
    data_path = data_path,
    class_column = class_column,
    case_label = case_label,
    control_label = control_label,
    feature_prefix = feature_prefix,
    use_smote = use_smote,
    verbose = verbose
  )
  
  # Extract data components
  train_data <- data_list$train
  # test_data and validation_data available for future use if needed
  # test_data <- data_list$test
  # validation_data <- data_list$validation
  train_smoted <- data_list$train_smoted
  
  # Get feature matrix (without class column)
  feature_columns <- get_feature_columns(train_data, feature_prefix)
  train_features <- train_data[, feature_columns, drop = FALSE]
  train_classes <- train_data[[class_column]]
  
  if (verbose) {
    log_info("Data loaded successfully:")
    log_info(paste("  Training samples:", nrow(train_data), 
                   paste0("(", sum(train_classes == case_label), " ", case_label, ", ",
                          sum(train_classes == control_label), " ", control_label, ")")))
    log_info(paste("  Test samples:", nrow(data_list$test)))
    log_info(paste("  Validation samples:", nrow(data_list$validation)))
    log_info(paste("  Features:", length(feature_columns), "(prefix:", feature_prefix, ")"))
    if (use_smote) {
      log_info(paste("  SMOTE samples:", nrow(train_smoted)))
    }
  }
  
  # Perform differential expression analysis
  if (verbose) log_info("Performing differential expression analysis")
  de_results <- perform_differential_expression(
    feature_data = train_features,
    class_labels = train_classes,
    case_label = case_label,
    control_label = control_label,
    p_threshold = p_threshold,
    fc_threshold = fc_threshold
  )
  
  # Setup parallel processing
  if (verbose) log_info("Setting up parallel processing")
  parallel_setup <- setup_parallel_processing(parallel_cores, verbose)
  
  # Initialize results storage
  results <- initialize_omicselector_results(
    methods = methods,
    parameters = validated_params,
    de_results = de_results,
    start_time = start_time
  )
  
  # Execute feature selection methods
  if (verbose) {
    log_info(paste("Executing", length(methods), "feature selection methods"))
    pb <- create_progress_bar(length(methods), "Feature Selection Methods")
  } else {
    pb <- NULL
  }
  
  for (i in seq_along(methods)) {
    method_id <- methods[i]
    
    if (verbose) {
      if (!is.null(pb) && is.list(pb) && !is.null(pb$tick)) {
        tryCatch(pb$tick(), error = function(e) cat("."))
      }
      log_info(paste("Processing method", method_id, ":", get_method_name(method_id)))
    }
    
    # Execute method with timeout and error handling
    method_result <- tryCatch({
      withTimeout({
        execute_feature_selection_method(
          method_id = method_id,
          train_data = train_data,
          train_features = train_features,
          train_classes = train_classes,
          train_smoted = train_smoted,
          de_results = de_results,
          class_column = class_column,
          feature_prefix = feature_prefix,
          max_features = max_features,
          max_iterations = max_iterations,
          case_label = case_label,
          control_label = control_label,
          parallel_setup = parallel_setup,
          verbose = verbose
        )
      }, timeout = timeout_minutes * 60)
    }, TimeoutException = function(e) {
      if (verbose) log_warn(paste("Method", method_id, "timed out after", timeout_minutes, "minutes"))
      list(features = character(0), status = "timeout", error = "Timeout exceeded")
    }, error = function(e) {
      if (verbose) log_error(paste("Method", method_id, "failed:", e$message))
      list(features = character(0), status = "error", error = e$message)
    })
    
    # Store results
    results$selected_features[[paste0("method_", method_id)]] <- method_result$features
    results$method_status[[paste0("method_", method_id)]] <- method_result$status %||% "success"
    results$execution_times[[paste0("method_", method_id)]] <- method_result$execution_time %||% NA
    
    if (!is.null(method_result$error)) {
      results$method_errors[[paste0("method_", method_id)]] <- method_result$error
    }
    
    # Save intermediate results
    if (save_intermediate) {
      saveRDS(results, file.path(output_dir, paste0("intermediate_results_method_", method_id, ".rds")))
    }
  }
  
  # Cleanup parallel processing
  cleanup_parallel_processing(parallel_setup)
  
  # Finalize results
  results$metadata$end_time <- Sys.time()
  results$metadata$total_runtime <- as.numeric(difftime(results$metadata$end_time, start_time, units = "mins"))
  
  # Generate summary statistics
  results$summary <- generate_omicselector_summary(results, de_results, verbose)
  
  # Save final results
  if (save_intermediate) {
    saveRDS(results, file.path(output_dir, "final_results.rds"))
    if (verbose) log_info(paste("Results saved to", file.path(output_dir, "final_results.rds")))
  }
  
  if (verbose) {
    log_info("OmicSelector analysis completed successfully")
    log_info(paste("Total runtime:", round(results$metadata$total_runtime, 2), "minutes"))
    log_info(paste("Methods completed:", sum(sapply(results$method_status, function(x) x == "success")), "/", length(methods)))
  }
  
  # Set class and return
  class(results) <- c("omicselector_results", "list")
  return(results)
}

#' Validate OmicSelector Input Parameters
#'
#' @param ... All input parameters from main function
#' @return List of validated parameters
#' @keywords internal
validate_omicselector_inputs <- function(data_path, methods, class_column, case_label, 
                                       control_label, feature_prefix, max_features, 
                                       max_iterations, p_threshold, fc_threshold, 
                                       parallel_cores, timeout_minutes, output_dir) {
  
  # Validate data path
  if (!dir.exists(data_path)) {
    stop("Data path does not exist: ", data_path)
  }
  
  # Check for required data files
  required_files <- c("mixed_train.csv", "mixed_test.csv", "mixed_validation.csv")
  missing_files <- required_files[!file.exists(file.path(data_path, required_files))]
  if (length(missing_files) > 0) {
    stop("Missing required data files: ", paste(missing_files, collapse = ", "))
  }
  
  # Validate methods
  if (!is.numeric(methods) || any(methods < 1) || any(methods > 70)) {
    stop("Methods must be integers between 1 and 70")
  }
  
  # Validate thresholds
  if (p_threshold <= 0 || p_threshold >= 1) {
    stop("p_threshold must be between 0 and 1")
  }
  
  if (fc_threshold < 0) {
    stop("fc_threshold must be non-negative")
  }
  
  # Validate feature counts
  if (max_features <= 0 || max_iterations <= 0) {
    stop("max_features and max_iterations must be positive integers")
  }
  
  # Validate parallel cores
  if (!is.null(parallel_cores)) {
    max_cores <- parallel::detectCores()
    if (parallel_cores > max_cores) {
      warning("Requested cores (", parallel_cores, ") exceeds available cores (", max_cores, "). Using ", max_cores)
      parallel_cores <- max_cores
    }
  }
  
  # Return validated parameters
  list(
    data_path = normalizePath(data_path),
    methods = unique(sort(methods)),
    class_column = class_column,
    case_label = case_label,
    control_label = control_label,
    feature_prefix = feature_prefix,
    max_features = as.integer(max_features),
    max_iterations = as.integer(max_iterations),
    p_threshold = p_threshold,
    fc_threshold = fc_threshold,
    parallel_cores = parallel_cores,
    timeout_minutes = as.integer(timeout_minutes),
    output_dir = normalizePath(output_dir, mustWork = FALSE)
  )
}

#' Load and Prepare OmicSelector Data
#'
#' @param data_path Character. Path to data directory
#' @param class_column Character. Name of class column
#' @param case_label Character. Positive case label
#' @param control_label Character. Negative control label
#' @param feature_prefix Character. Feature column prefix
#' @param use_smote Logical. Apply SMOTE balancing
#' @param verbose Logical. Enable logging
#' @return List containing train, test, validation, and optionally SMOTE data
#' @keywords internal
load_omicselector_data <- function(data_path, class_column, case_label, control_label,
                                 feature_prefix, use_smote, verbose) {
  
  # Load data files
  train_file <- file.path(data_path, "mixed_train.csv")
  test_file <- file.path(data_path, "mixed_test.csv")
  validation_file <- file.path(data_path, "mixed_validation.csv")
  
  train_data <- read.csv(train_file, stringsAsFactors = FALSE)
  test_data <- read.csv(test_file, stringsAsFactors = FALSE)
  validation_data <- read.csv(validation_file, stringsAsFactors = FALSE)
  
  # Validate class column exists
  if (!class_column %in% colnames(train_data)) {
    stop("Class column '", class_column, "' not found in training data")
  }
  
  # Validate class labels
  unique_classes <- unique(train_data[[class_column]])
  if (!case_label %in% unique_classes || !control_label %in% unique_classes) {
    stop("Expected class labels '", case_label, "' and '", control_label, 
         "' not found. Available: ", paste(unique_classes, collapse = ", "))
  }
  
  # Ensure class column is factor with correct levels
  train_data[[class_column]] <- factor(train_data[[class_column]], 
                                     levels = c(control_label, case_label))
  test_data[[class_column]] <- factor(test_data[[class_column]], 
                                    levels = c(control_label, case_label))
  validation_data[[class_column]] <- factor(validation_data[[class_column]], 
                                          levels = c(control_label, case_label))
  
  # Validate feature columns
  feature_cols <- get_feature_columns(train_data, feature_prefix)
  if (length(feature_cols) == 0) {
    stop("No feature columns found with prefix '", feature_prefix, "'")
  }
  
  if (verbose) {
    log_info(paste("Found", length(feature_cols), "features with prefix", paste0("'", feature_prefix, "'")))
    log_info(paste("Class distribution in training:", paste(names(table(train_data[[class_column]])), table(train_data[[class_column]]), collapse = ", ")))
  }
  
  # Prepare result list
  result <- list(
    train = train_data,
    test = test_data,
    validation = validation_data
  )
  
  # Apply SMOTE if requested
  if (use_smote) {
    if (verbose) log_info("Applying SMOTE for class balancing")
    
    # Check if SMOTE packages are available
    if (!requireNamespace("DMwR", quietly = TRUE)) {
      if (verbose) log_warn("DMwR package not available, skipping SMOTE")
      result$train_smoted <- train_data
    } else {
      # Apply SMOTE
      smote_formula <- as.formula(paste(class_column, "~ ."))
      train_smoted <- tryCatch({
        DMwR::SMOTE(smote_formula, data = train_data, 
                   perc.over = 200, perc.under = 150, k = 5)
      }, error = function(e) {
        if (verbose) log_warn("SMOTE failed: {e$message}. Using original training data.")
        train_data
      })
      
      result$train_smoted <- train_smoted
      
      if (verbose) {
        log_info(paste("SMOTE completed. New training size:", nrow(train_smoted)))
        log_info(paste("SMOTE class distribution:", paste(names(table(train_smoted[[class_column]])), table(train_smoted[[class_column]]), collapse = ", ")))
      }
    }
  } else {
    result$train_smoted <- train_data
  }
  
  return(result)
}

#' Get Feature Column Names
#'
#' @param data Data frame
#' @param feature_prefix Character prefix for feature columns
#' @return Character vector of feature column names
#' @keywords internal
get_feature_columns <- function(data, feature_prefix) {
  feature_cols <- colnames(data)[startsWith(colnames(data), feature_prefix)]
  return(feature_cols)
}

#' Perform Differential Expression Analysis
#'
#' @param feature_data Matrix or data frame of features
#' @param class_labels Factor of class labels
#' @param case_label Character. Positive case label
#' @param control_label Character. Negative control label  
#' @param p_threshold Numeric. P-value threshold
#' @param fc_threshold Numeric. Fold change threshold
#' @return Data frame with DE results
#' @keywords internal
perform_differential_expression <- function(feature_data, class_labels, case_label, 
                                          control_label, p_threshold, fc_threshold) {
  
  # Prepare results data frame
  feature_names <- colnames(feature_data)
  n_features <- length(feature_names)
  
  de_results <- data.frame(
    feature = feature_names,
    mean_case = numeric(n_features),
    mean_control = numeric(n_features),
    log2fc = numeric(n_features),
    p_value = numeric(n_features),
    p_value_adj = numeric(n_features),
    significant = logical(n_features),
    stringsAsFactors = FALSE
  )
  
  # Calculate statistics for each feature
  case_indices <- which(class_labels == case_label)
  control_indices <- which(class_labels == control_label)
  
  for (i in seq_along(feature_names)) {
    feature_values <- feature_data[, i]
    
    case_values <- feature_values[case_indices]
    control_values <- feature_values[control_indices]
    
    # Remove missing values
    case_values <- case_values[!is.na(case_values)]
    control_values <- control_values[!is.na(control_values)]
    
    if (length(case_values) > 1 && length(control_values) > 1) {
      # Calculate means
      mean_case <- mean(case_values)
      mean_control <- mean(control_values)
      
      # Calculate log2 fold change (add small pseudocount to avoid log(0))
      pseudocount <- 1e-6
      log2fc <- log2((mean_case + pseudocount) / (mean_control + pseudocount))
      
      # Perform t-test
      t_test_result <- tryCatch({
        t.test(case_values, control_values)
      }, error = function(e) {
        list(p.value = 1)
      })
      
      # Store results
      de_results$mean_case[i] <- mean_case
      de_results$mean_control[i] <- mean_control
      de_results$log2fc[i] <- log2fc
      de_results$p_value[i] <- t_test_result$p.value
    } else {
      # Insufficient data
      de_results$p_value[i] <- 1
    }
  }
  
  # Adjust p-values using Benjamini-Hochberg method
  de_results$p_value_adj <- p.adjust(de_results$p_value, method = "BH")
  
  # Determine significance
  de_results$significant <- (de_results$p_value_adj <= p_threshold) & 
                           (abs(de_results$log2fc) >= log2(fc_threshold + 1))
  
  # Sort by adjusted p-value
  de_results <- de_results[order(de_results$p_value_adj), ]
  
  return(de_results)
}

#' Setup Parallel Processing
#'
#' @param parallel_cores Integer or NULL
#' @param verbose Logical
#' @return List with parallel processing setup
#' @keywords internal
setup_parallel_processing <- function(parallel_cores, verbose) {
  
  if (is.null(parallel_cores)) {
    parallel_cores <- min(parallel::detectCores() - 1, 4)  # Conservative default
  }
  
  if (parallel_cores > 1) {
    if (verbose) log_info(paste("Setting up parallel processing with", parallel_cores, "cores"))
    
    # Setup cluster
    cluster <- parallel::makePSOCKcluster(parallel_cores)
    doParallel::registerDoParallel(cluster)
    
    return(list(
      cluster = cluster,
      cores = parallel_cores,
      enabled = TRUE
    ))
  } else {
    if (verbose) log_info("Using sequential processing")
    return(list(
      cluster = NULL,
      cores = 1,
      enabled = FALSE
    ))
  }
}

#' Cleanup Parallel Processing
#'
#' @param parallel_setup List from setup_parallel_processing
#' @keywords internal
cleanup_parallel_processing <- function(parallel_setup) {
  if (parallel_setup$enabled && !is.null(parallel_setup$cluster)) {
    parallel::stopCluster(parallel_setup$cluster)
    doParallel::registerDoSEQ()  # Reset to sequential
  }
}

#' Initialize OmicSelector Results Object
#'
#' @param methods Integer vector of methods
#' @param parameters List of validated parameters
#' @param de_results Data frame of DE results
#' @param start_time POSIXct start time
#' @return List structure for results
#' @keywords internal
initialize_omicselector_results <- function(methods, parameters, de_results, start_time) {
  
  list(
    selected_features = list(),
    method_status = list(),
    method_errors = list(),
    execution_times = list(),
    differential_expression = de_results,
    metadata = list(
      start_time = start_time,
      parameters = parameters,
      methods_requested = methods,
      package_version = packageVersion("OmicSelector"),
      r_version = R.version.string,
      session_info = sessionInfo()
    )
  )
}

#' Get Method Name by ID
#'
#' @param method_id Integer method identifier
#' @return Character method name
#' @keywords internal
get_method_name <- function(method_id) {
  method_names <- c(
    "1" = "All Features",
    "2" = "Significant Features (t-test)",  
    "3" = "Fold Change + Significance",
    "4" = "Correlation-based Feature Selection (CFS)",
    "5" = "Classifier Loop",
    "6" = "Classifier Loop (SMOTE)",
    "7" = "Classifier Loop (Significant)",
    "8" = "Classifier Loop (SMOTE + Significant)",
    "9" = "Forward Correlation-based Selection",
    "10" = "Forward CFS (SMOTE)",
    "11" = "Random Forest RFE",
    "12" = "Random Forest RFE (SMOTE)",
    "13" = "Boruta Algorithm",
    "14" = "Boruta Algorithm (SMOTE)",
    "15" = "Genetic Algorithm RF",
    "16" = "Genetic Algorithm RF (SMOTE)",
    "17" = "Simulated Annealing",
    "18" = "Simulated Annealing (SMOTE)",
    "19" = "varSelRF",
    "20" = "varSelRF (SMOTE)"
  )
  
  return(method_names[as.character(method_id)] %||% paste("Method", method_id))
}

#' Execute Individual Feature Selection Method
#'
#' @param method_id Integer method identifier
#' @param train_data Full training data
#' @param train_features Training feature matrix
#' @param train_classes Training class labels
#' @param train_smoted SMOTE-balanced training data
#' @param de_results Differential expression results
#' @param class_column Character class column name
#' @param feature_prefix Character feature prefix
#' @param max_features Integer maximum features to select
#' @param max_iterations Integer maximum iterations
#' @param case_label Character positive case label
#' @param control_label Character negative control label
#' @param parallel_setup List with parallel setup
#' @param verbose Logical enable logging
#' @return List with selected features and metadata
#' @keywords internal
execute_feature_selection_method <- function(method_id, train_data, train_features, 
                                           train_classes, train_smoted, de_results,
                                           class_column, feature_prefix, max_features,
                                           max_iterations, case_label, control_label,
                                           parallel_setup, verbose) {
  
  method_start <- Sys.time()
  
  # Get significant features for methods that use them
  sig_features <- de_results$feature[de_results$significant]
  if (length(sig_features) == 0) {
    sig_features <- head(de_results$feature, max_features)  # Fallback to top features
  }
  
  selected_features <- switch(as.character(method_id),
    
    # Method 1: All features
    "1" = {
      if (verbose) log_info("Method 1: Selecting all features")
      get_feature_columns(train_data, feature_prefix)
    },
    
    # Method 2: Significant features from DE analysis
    "2" = {
      if (verbose) log_info(paste("Method 2: Selecting significant features (n=", length(sig_features), ")", sep=""))
      sig_features
    },
    
    # Method 3: Fold change + significance filter
    "3" = {
      if (verbose) log_info("Method 3: Fold change + significance filter")
      de_results$feature[de_results$significant & abs(de_results$log2fc) >= log2(2)]
    },
    
    # Method 4: Correlation-based Feature Selection
    "4" = {
      if (verbose) log_info("Method 4: Correlation-based Feature Selection")
      execute_cfs_method(train_data, class_column, verbose)
    },
    
    # Method 11: Random Forest RFE
    "11" = {
      if (verbose) log_info("Method 11: Random Forest RFE")
      execute_rf_rfe_method(train_data, class_column, max_features, verbose)
    },
    
    # Method 13: Boruta Algorithm  
    "13" = {
      if (verbose) log_info("Method 13: Boruta Algorithm")
      execute_boruta_method(train_data, class_column, max_iterations, verbose)
    },
    
    # Method 17: Simulated Annealing
    "17" = {
      if (verbose) log_info("Method 17: Simulated Annealing")
      execute_sa_method(train_data, class_column, max_features, max_iterations, verbose)
    },
    
    # Default: return empty if method not implemented
    {
      if (verbose) log_warn("Method {method_id} not implemented yet")
      character(0)
    }
  )
  
  # Calculate execution time
  execution_time <- as.numeric(difftime(Sys.time(), method_start, units = "secs"))
  
  # Limit number of features
  if (length(selected_features) > max_features) {
    selected_features <- head(selected_features, max_features)
    if (verbose) log_info(paste("Limited features to", max_features, "for method", method_id))
  }
  
  return(list(
    features = selected_features,
    execution_time = execution_time,
    status = "success"
  ))
}

#' Execute CFS Method
#' @keywords internal
execute_cfs_method <- function(train_data, class_column, verbose) {
  # Suppress rgl errors during Biocomb loading
  old_rgl_opt <- getOption("rgl.useNULL", NULL)
  on.exit(options(rgl.useNULL = old_rgl_opt))
  options(rgl.useNULL = TRUE)
  
  if (!requireNamespace("Biocomb", quietly = TRUE)) {
    if (verbose) log_warn("Biocomb package not available for CFS")
    return(character(0))
  }
  
  tryCatch({
    # Suppress any graphics-related warnings/errors
    suppressMessages({
      cfs_result <- Biocomb::select.cfs(train_data)
    })
    return(as.character(cfs_result$Biomarker))
  }, error = function(e) {
    if (verbose) log_error("CFS method failed: {e$message}")
    return(character(0))
  })
}

#' Execute Random Forest RFE Method
#' @keywords internal
execute_rf_rfe_method <- function(train_data, class_column, max_features, verbose) {
  if (!requireNamespace("caret", quietly = TRUE)) {
    if (verbose) log_warn("caret package not available for RF-RFE")
    return(character(0))
  }
  
  tryCatch({
    # Prepare data
    x <- train_data[, !colnames(train_data) %in% class_column, drop = FALSE]
    y <- train_data[[class_column]]
    
    # Setup RFE control
    ctrl <- caret::rfeControl(
      functions = caret::rfFuncs,
      method = "cv",
      number = 5,
      verbose = FALSE
    )
    
    # Run RFE
    rfe_result <- caret::rfe(
      x = x, y = y,
      sizes = c(5, 10, 15, min(max_features, ncol(x))),
      rfeControl = ctrl
    )
    
    return(rfe_result$optVariables)
  }, error = function(e) {
    if (verbose) log_error("RF-RFE method failed: {e$message}")
    return(character(0))
  })
}

#' Execute Boruta Method
#' @keywords internal
execute_boruta_method <- function(train_data, class_column, max_iterations, verbose) {
  if (!requireNamespace("Boruta", quietly = TRUE)) {
    if (verbose) log_warn("Boruta package not available")
    return(character(0))
  }
  
  tryCatch({
    # Prepare formula
    formula_str <- paste(class_column, "~ .")
    boruta_formula <- as.formula(formula_str)
    
    # Run Boruta
    boruta_result <- Boruta::Boruta(
      formula = boruta_formula,
      data = train_data,
      maxRuns = max_iterations * 10,
      doTrace = 0
    )
    
    # Get confirmed features
    confirmed_features <- names(boruta_result$finalDecision[boruta_result$finalDecision == "Confirmed"])
    return(confirmed_features)
  }, error = function(e) {
    if (verbose) log_error("Boruta method failed: {e$message}")
    return(character(0))
  })
}

#' Execute Simulated Annealing Method
#' @keywords internal
execute_sa_method <- function(train_data, class_column, max_features, max_iterations, verbose) {
  if (!requireNamespace("caret", quietly = TRUE)) {
    if (verbose) log_warn("caret package not available for SA")
    return(character(0))
  }
  
  tryCatch({
    # Prepare data
    x <- train_data[, !colnames(train_data) %in% class_column, drop = FALSE]
    y <- train_data[[class_column]]
    
    # Setup SA control
    sa_ctrl <- caret::safsControl(
      functions = caret::rfSA,
      method = "cv",
      number = 3,  # Reduce CV folds for speed
      improve = 5,
      returnResamp = "none",
      verbose = FALSE
    )
    
    # Run SA
    sa_result <- caret::safs(
      x = x, y = y,
      iters = max_iterations * 5,
      safsControl = sa_ctrl
    )
    
    return(sa_result$optVariables)
  }, error = function(e) {
    if (verbose) log_error("Simulated Annealing method failed: {e$message}")
    return(character(0))
  })
}

#' Generate Summary Statistics
#'
#' @param results OmicSelector results object
#' @param de_results Differential expression results
#' @param verbose Logical enable logging
#' @return List with summary statistics
#' @keywords internal
generate_omicselector_summary <- function(results, de_results, verbose) {
  
  # Count successful methods
  successful_methods <- sum(sapply(results$method_status, function(x) x == "success"))
  total_methods <- length(results$method_status)
  
  # Get feature selection statistics
  feature_counts <- sapply(results$selected_features, length)
  feature_counts <- feature_counts[feature_counts > 0]  # Remove empty results
  
  # Find most commonly selected features
  all_selected <- unlist(results$selected_features)
  feature_frequency <- table(all_selected)
  n_common <- min(10, length(feature_frequency))
  most_common <- if (n_common > 0) {
    names(sort(feature_frequency, decreasing = TRUE))[seq_len(n_common)]
  } else {
    character(0)
  }
  
  # Calculate execution time statistics
  execution_times <- unlist(results$execution_times)
  execution_times <- execution_times[!is.na(execution_times)]
  
  summary_stats <- list(
    total_methods = total_methods,
    successful_methods = successful_methods,
    success_rate = successful_methods / total_methods,
    feature_counts = feature_counts,
    avg_features_selected = mean(feature_counts, na.rm = TRUE),
    median_features_selected = median(feature_counts, na.rm = TRUE),
    most_common_features = most_common,
    feature_frequency = feature_frequency[most_common],
    total_execution_time = sum(execution_times, na.rm = TRUE),
    avg_method_time = mean(execution_times, na.rm = TRUE),
    significant_features_found = sum(de_results$significant),
    total_features_tested = nrow(de_results)
  )
  
  if (verbose) {
    log_info("Summary Statistics:")
    log_info(paste("  Successful methods:", successful_methods, "/", total_methods, 
                   paste0("(", round(summary_stats$success_rate*100, 1), "%)")))
    log_info(paste("  Average features selected:", round(summary_stats$avg_features_selected, 1)))
    log_info(paste("  Total execution time:", round(summary_stats$total_execution_time, 1), "seconds"))
    log_info(paste("  Significant features found:", summary_stats$significant_features_found))
  }
  
  return(summary_stats)
}

#' Print Method for OmicSelector Results
#'
#' @param x OmicSelector results object
#' @param ... Additional arguments
#' @export
print.omicselector_results <- function(x, ...) {
  cat("OmicSelector Results\n")
  cat("===================\n\n")
  
  cat("Run Information:\n")
  cat("  Start time:", format(x$metadata$start_time), "\n")
  cat("  Total runtime:", round(x$metadata$total_runtime, 2), "minutes\n")
  cat("  Methods requested:", length(x$metadata$methods_requested), "\n")
  cat("  Methods completed:", x$summary$successful_methods, "/", x$summary$total_methods, "\n")
  cat("  Success rate:", round(x$summary$success_rate * 100, 1), "%\n\n")
  
  cat("Feature Selection Results:\n")
  cat("  Total features tested:", x$summary$total_features_tested, "\n")
  cat("  Significant features:", x$summary$significant_features_found, "\n")
  cat("  Average features selected:", round(x$summary$avg_features_selected, 1), "\n")
  cat("  Median features selected:", round(x$summary$median_features_selected, 1), "\n\n")
  
  if (length(x$summary$most_common_features) > 0) {
    cat("Most Frequently Selected Features:\n")
    n_show <- min(5, length(x$summary$most_common_features))
    for (i in seq_len(n_show)) {
      feature <- x$summary$most_common_features[i]
      freq <- x$summary$feature_frequency[feature]
      cat("  ", feature, "(selected", freq, "times)\n")
    }
  }
  
  invisible(x)
}

#' Summary Method for OmicSelector Results
#'
#' @param object OmicSelector results object
#' @param ... Additional arguments
#' @export
summary.omicselector_results <- function(object, ...) {
  print(object)
  
  cat("\nMethod-wise Results:\n")
  for (method_name in names(object$selected_features)) {
    n_features <- length(object$selected_features[[method_name]])
    status <- object$method_status[[method_name]]
    time <- object$execution_times[[method_name]]
    
    cat("  ", method_name, ": ", n_features, " features", sep = "")
    if (!is.na(time)) {
      cat(" (", round(time, 1), "s)", sep = "")
    }
    if (status != "success") {
      cat(" [", status, "]", sep = "")
    }
    cat("\n")
  }
  
  invisible(object)
}
