#' Modern Logging System for OmicSelector
#'
#' @description
#' Provides centralized logging functionality with different levels and outputs
#'
#' @param message Character string with the log message
#' @param level Character string with log level (TRACE, DEBUG, INFO, WARN, ERROR, FATAL)
#' @param logger_name Character string with logger name, defaults to "OmicSelector"
#' @param ... Additional parameters passed to logging function
#'
#' @examples
#' \dontrun{
#' log_info("Starting feature selection")
#' log_warn("Some features were removed due to low variance")
#' log_error("Failed to load data file")
#' }
#'
#' @export
log_info <- function(message, ...) {
  omics_log("INFO", message, ...)
}

#' @rdname log_info
#' @export
log_warn <- function(message, ...) {
  omics_log("WARN", message, ...)
}

#' @rdname log_info
#' @export
log_error <- function(message, ...) {
  omics_log("ERROR", message, ...)
}

#' @rdname log_info
#' @export
log_debug <- function(message, ...) {
  omics_log("DEBUG", message, ...)
}

#' Internal logging function
#'
#' @param level Log level
#' @param message Log message
#' @param ... Additional parameters
#' @keywords internal
omics_log <- function(level, message, ...) {
  
  # Check if logger package is available
  if (requireNamespace("logger", quietly = TRUE)) {
    switch(level,
      "TRACE" = logger::log_trace(message, ...),
      "DEBUG" = logger::log_debug(message, ...),
      "INFO" = logger::log_info(message, ...),
      "WARN" = logger::log_warn(message, ...),
      "ERROR" = logger::log_error(message, ...),
      "FATAL" = logger::log_fatal(message, ...)
    )
  } else {
    # Fallback to base R messaging
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    formatted_msg <- paste0("[", timestamp, "] [", level, "] ", message)
    
    switch(level,
      "ERROR" = stop(formatted_msg, call. = FALSE),
      "WARN" = warning(formatted_msg, call. = FALSE),
      message(formatted_msg)
    )
  }
}

#' Setup Logging Configuration
#'
#' @description
#' Initializes logging system with specified configuration
#'
#' @param log_level Character string specifying minimum log level
#' @param log_file Character string with path to log file, NULL for console only
#' @param logger_name Character string with logger name
#'
#' @examples
#' \dontrun{
#' setup_logging("DEBUG", "omicselector.log")
#' setup_logging("INFO")  # Console only
#' }
#'
#' @export
setup_logging <- function(log_level = "INFO", log_file = NULL, logger_name = "OmicSelector") {
  
  if (requireNamespace("logger", quietly = TRUE)) {
    # Set log level
    logger::log_threshold(log_level, namespace = logger_name)
    
    # Setup file logging if specified
    if (!is.null(log_file)) {
      logger::log_appender(logger::appender_tee(log_file), namespace = logger_name)
    }
    
    # Setup console formatting
    logger::log_formatter(logger::formatter_glue_or_sprintf, namespace = logger_name)
    
    log_info("Logging system initialized", level = log_level, file = log_file %||% "console")
  } else {
    message("Logger package not available, using basic logging")
  }
}

#' Progress Bar Utilities
#'
#' @description
#' Modern progress reporting using cli package when available
#'
#' @param total Integer, total number of steps
#' @param message Character string with progress message
#' @param clear Logical, whether to clear progress bar when complete
#'
#' @examples
#' \dontrun{
#' pb <- create_progress_bar(100, "Processing features")
#' for (i in 1:100) {
#'   update_progress_bar(pb, i)
#'   Sys.sleep(0.1)
#' }
#' finish_progress_bar(pb)
#' }
#'
#' @export
create_progress_bar <- function(total, message = "Processing", clear = TRUE) {
  
  if (requireNamespace("cli", quietly = TRUE)) {
    cli::cli_progress_bar(
      total = total,
      format = paste("{message} {cli::pb_bar} {cli::pb_percent} [{cli::pb_elapsed}]"),
      clear = clear
    )
  } else {
    # Fallback progress bar
    list(
      total = total,
      current = 0,
      message = message,
      type = "base"
    )
  }
}

#' @rdname create_progress_bar
#' @export
update_progress_bar <- function(pb, value = NULL, increment = 1) {
  
  if (is.list(pb) && pb$type == "base") {
    # Base R progress bar
    pb$current <- pb$current + increment
    if (pb$current %% 10 == 0 || pb$current == pb$total) {
      cat("\r", pb$message, ": ", pb$current, "/", pb$total, 
          " (", round(100 * pb$current / pb$total, 1), "%)", sep = "")
      if (pb$current == pb$total) cat("\n")
    }
  } else if (requireNamespace("cli", quietly = TRUE)) {
    if (is.null(value)) {
      cli::cli_progress_update(id = pb, inc = increment)
    } else {
      cli::cli_progress_update(id = pb, set = value)
    }
  }
}

#' @rdname create_progress_bar
#' @export
finish_progress_bar <- function(pb) {
  
  if (requireNamespace("cli", quietly = TRUE) && !is.list(pb)) {
    cli::cli_progress_done(id = pb)
  } else if (is.list(pb) && pb$type == "base") {
    cat("\n", pb$message, " completed!\n")
  }
}

# Utility function for null-coalescing
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
