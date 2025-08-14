#' Input Validation and Error Handling
#'
#' @description
#' Comprehensive input validation functions for OmicSelector
#'
#' @name validation
NULL

#' Validate Working Directory
#'
#' @param wd Character string with working directory path
#' @param required_files Character vector with required file names
#'
#' @return Logical, TRUE if valid
#' @keywords internal
validate_working_directory <- function(wd, required_files = c("mixed_train.csv", "mixed_test.csv")) {
  
  if (!dir.exists(wd)) {
    stop("Working directory does not exist: ", wd, call. = FALSE)
  }
  
  # Check for required files
  missing_files <- character(0)
  for (file in required_files) {
    if (!file.exists(file.path(wd, file))) {
      missing_files <- c(missing_files, file)
    }
  }
  
  if (length(missing_files) > 0) {
    stop("Required files missing from working directory: ", 
         paste(missing_files, collapse = ", "), call. = FALSE)
  }
  
  TRUE
}

#' Validate Feature Selection Methods
#'
#' @param methods Integer vector with method numbers
#' @param max_method Integer, maximum valid method number
#'
#' @return Validated integer vector
#' @keywords internal
validate_methods <- function(methods, max_method = 70) {
  
  if (!is.numeric(methods)) {
    stop("Methods must be numeric vector", call. = FALSE)
  }
  
  methods <- as.integer(methods)
  
  if (any(methods < 1) || any(methods > max_method)) {
    stop("Methods must be between 1 and ", max_method, call. = FALSE)
  }
  
  if (any(duplicated(methods))) {
    warning("Duplicate methods detected, removing duplicates")
    methods <- unique(methods)
  }
  
  methods
}

#' Validate Data Frame
#'
#' @param data Data frame to validate
#' @param min_rows Integer, minimum number of rows required
#' @param min_cols Integer, minimum number of columns required
#' @param required_cols Character vector, names of required columns
#'
#' @return Logical, TRUE if valid
#' @keywords internal
validate_data_frame <- function(data, min_rows = 10, min_cols = 2, required_cols = NULL) {
  
  if (!is.data.frame(data)) {
    stop("Data must be a data frame", call. = FALSE)
  }
  
  if (nrow(data) < min_rows) {
    stop("Data must have at least ", min_rows, " rows", call. = FALSE)
  }
  
  if (ncol(data) < min_cols) {
    stop("Data must have at least ", min_cols, " columns", call. = FALSE)
  }
  
  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
      stop("Required columns missing: ", paste(missing_cols, collapse = ", "), 
           call. = FALSE)
    }
  }
  
  TRUE
}

#' Check Optional Dependencies
#'
#' @param packages Character vector with package names
#' @param action Character string: "stop", "warn", or "skip"
#'
#' @return Logical vector indicating which packages are available
#' @keywords internal
check_optional_dependencies <- function(packages, action = "warn") {
  
  available <- sapply(packages, function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
  })
  
  missing <- packages[!available]
  
  if (length(missing) > 0) {
    msg <- paste("Optional packages not available:", paste(missing, collapse = ", "))
    
    switch(action,
      "stop" = stop(msg, call. = FALSE),
      "warn" = warning(msg, call. = FALSE),
      "skip" = message(msg)
    )
  }
  
  available
}

#' Validate Parallel Configuration
#'
#' @param cores Integer, number of cores to use
#' @param register_parallel Logical, whether to register parallel backend
#'
#' @return List with validated parallel configuration
#' @keywords internal
validate_parallel_config <- function(cores, register_parallel = TRUE) {
  
  max_cores <- parallel::detectCores()
  if (is.na(max_cores)) max_cores <- 1
  
  if (cores > max_cores) {
    warning("Requested cores (", cores, ") exceeds available (", max_cores, 
            "), using ", max_cores)
    cores <- max_cores
  }
  
  if (cores < 1) {
    warning("Cores must be at least 1, setting to 1")
    cores <- 1
  }
  
  list(
    cores = cores,
    register_parallel = register_parallel,
    max_cores = max_cores
  )
}

#' Safe File Operations
#'
#' @description
#' Wrapper functions for safe file operations with proper error handling
#'
#' @param file Character string with file path
#' @param data Data to save
#' @param ... Additional arguments passed to save function
#'
#' @return Logical, TRUE if successful
#' @keywords internal
safe_save_rds <- function(file, data, ...) {
  tryCatch({
    saveRDS(data, file, ...)
    log_debug("Saved RDS file: ", file)
    TRUE
  }, error = function(e) {
    log_error("Failed to save RDS file ", file, ": ", e$message)
    FALSE
  })
}

#' @rdname safe_save_rds
#' @keywords internal
safe_read_rds <- function(file) {
  tryCatch({
    data <- readRDS(file)
    log_debug("Loaded RDS file: ", file)
    data
  }, error = function(e) {
    log_error("Failed to read RDS file ", file, ": ", e$message)
    NULL
  })
}

#' @rdname safe_save_rds
#' @keywords internal
safe_read_csv <- function(file, ...) {
  tryCatch({
    if (requireNamespace("readr", quietly = TRUE)) {
      data <- readr::read_csv(file, ...)
    } else {
      data <- read.csv(file, ...)
    }
    log_debug("Loaded CSV file: ", file)
    data
  }, error = function(e) {
    log_error("Failed to read CSV file ", file, ": ", e$message)
    NULL
  })
}
