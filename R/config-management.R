#' Configuration Management for OmicSelector
#'
#' @description
#' This function provides centralized configuration management for OmicSelector.
#' It loads default settings and allows for environment-specific overrides.
#'
#' @param config_name Character string specifying the configuration environment
#'   (default, development, hpc). Default is "default".
#' @param config_file Path to custom configuration file. If NULL, uses built-in
#'   configuration.
#' @param ... Additional configuration parameters to override
#'
#' @return A list containing configuration parameters
#'
#' @examples
#' \dontrun{
#' # Load default configuration
#' config <- get_config()
#' 
#' # Load development configuration
#' config_dev <- get_config("development")
#' 
#' # Override specific parameters
#' config_custom <- get_config(max_iterations = 5, cores = 2)
#' }
#'
#' @export
get_config <- function(config_name = "default", config_file = NULL, ...) {
  
  # Try to load config package, fallback to defaults if not available
  if (requireNamespace("config", quietly = TRUE)) {
    if (is.null(config_file)) {
      config_file <- system.file("config", "default.yml", package = "OmicSelector")
    }
    
    if (file.exists(config_file)) {
      cfg <- config::get(config = config_name, file = config_file)
    } else {
      cfg <- get_default_config()
    }
  } else {
    cfg <- get_default_config()
  }
  
  # Override with any additional parameters
  override_params <- list(...)
  if (length(override_params) > 0) {
    cfg <- modifyList(cfg, override_params)
  }
  
  return(cfg)
}

#' Get Default Configuration
#'
#' @description
#' Returns hardcoded default configuration when config package is not available
#'
#' @return List of default configuration parameters
#' @keywords internal
get_default_config <- function() {
  list(
    max_iterations = 10,
    timeout_sec = 172800,
    prefer_no_features = 11,
    cores = max(1, parallel::detectCores() - 1),
    log_level = "INFO",
    log_to_file = TRUE,
    default_methods = c(1, 2, 3, 11, 17, 20),
    search_iters = 2000,
    keras_epochs = 5000,
    holdout_validation = TRUE,
    register_parallel = TRUE,
    use_future = TRUE,
    plot_theme = "minimal",
    save_plots = TRUE,
    plot_format = "png",
    plot_dpi = 300
  )
}

#' Validate Configuration
#'
#' @description
#' Validates configuration parameters and ensures they are within acceptable ranges
#'
#' @param config Configuration list to validate
#'
#' @return Validated configuration list
#' @keywords internal
validate_config <- function(config) {
  
  # Validate numeric parameters
  if (!is.numeric(config$max_iterations) || config$max_iterations < 1) {
    stop("max_iterations must be a positive integer")
  }
  
  if (!is.numeric(config$timeout_sec) || config$timeout_sec < 1) {
    stop("timeout_sec must be a positive number")
  }
  
  if (!is.numeric(config$cores) || config$cores < 1) {
    config$cores <- 1
    warning("cores must be a positive integer, setting to 1")
  }
  
  # Ensure cores doesn't exceed available
  max_cores <- parallel::detectCores()
  if (config$cores > max_cores) {
    config$cores <- max_cores
    warning(paste("cores reduced to maximum available:", max_cores))
  }
  
  return(config)
}
