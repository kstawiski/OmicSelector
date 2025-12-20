#' @title Parallelization Support for OmicSelector 2.0
#'
#' @description
#' Functions for configuring and using future-based parallelization.
#' Uses the `future` package for backend-agnostic parallel execution.
#'
#' @details
#' OmicSelector uses `future` for parallelization because:
#' - Works with local cores, HPC clusters, and cloud backends
#' - Integrates natively with mlr3 ecosystem
#' - Supports lazy evaluation and proper error handling
#'
#' @name parallel
NULL


#' @title Configure Parallelization for OmicSelector
#'
#' @description
#' Sets up parallel processing using the future package.
#' Call this before running benchmarks or stability computations.
#'
#' @param workers Number of workers (cores) to use. Default is available cores - 1.
#' @param plan Parallelization plan: "multisession" (default, works everywhere),
#'   "multicore" (Unix only, more memory efficient), "sequential" (no parallelism)
#' @param memory_limit Per-worker memory limit in bytes (future.globals.maxSize)
#'
#' @return Invisible previous plan (for restoration)
#'
#' @examples
#' \dontrun{
#' # Use 4 cores
#' setup_parallel(workers = 4)
#'
#' # Run your analysis
#' result <- service$run(parallel = TRUE)
#'
#' # Reset to sequential
#' reset_parallel()
#' }
#'
#' @export
setup_parallel <- function(workers = NULL,
                           plan = c("multisession", "multicore", "sequential"),
                           memory_limit = 500 * 1024^2) {

  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' is required for parallelization. Install with: install.packages('future')",
         call. = FALSE)
  }

  plan <- match.arg(plan)

  # Default workers: available cores minus 1
  if (is.null(workers)) {
    workers <- max(1, parallel::detectCores() - 1)
  }

  # Set memory limit for globals
  options(future.globals.maxSize = memory_limit)

  # Store previous plan for restoration
  old_plan <- future::plan()

  if (plan == "sequential") {
    future::plan(future::sequential)
    message("Parallelization disabled (sequential mode)")
  } else if (plan == "multicore") {
    if (.Platform$OS.type == "windows") {
      warning("'multicore' not available on Windows. Using 'multisession' instead.")
      future::plan(future::multisession, workers = workers)
    } else {
      future::plan(future::multicore, workers = workers)
    }
    message(sprintf("Parallelization enabled: %s with %d workers", plan, workers))
  } else {
    future::plan(future::multisession, workers = workers)
    message(sprintf("Parallelization enabled: multisession with %d workers", workers))
  }

  invisible(old_plan)
}


#' @title Reset Parallelization to Sequential
#'
#' @description
#' Resets the future plan to sequential processing.
#' Useful for cleanup after parallel runs.
#'
#' @return Invisible previous plan
#'
#' @export
reset_parallel <- function() {
  if (!requireNamespace("future", quietly = TRUE)) {
    return(invisible(NULL))
  }

  old_plan <- future::plan()
  future::plan(future::sequential)
  message("Parallelization reset to sequential mode")
  invisible(old_plan)
}


#' @title Get Current Parallelization Status
#'
#' @description
#' Returns information about the current parallelization configuration.
#'
#' @return A list with parallelization details
#'
#' @export
get_parallel_status <- function() {
  if (!requireNamespace("future", quietly = TRUE)) {
    return(list(
      available = FALSE,
      plan = "sequential",
      workers = 1,
      message = "Package 'future' not installed"
    ))
  }

  plan_info <- future::plan()
  plan_name <- class(plan_info)[1]

  # Try to get number of workers
  workers <- tryCatch({
    future::nbrOfWorkers()
  }, error = function(e) 1)

  list(
    available = TRUE,
    plan = plan_name,
    workers = workers,
    is_parallel = !inherits(plan_info, "sequential"),
    memory_limit = getOption("future.globals.maxSize", NA)
  )
}


#' @title Run Function in Parallel with Split-Aware Keys
#'
#' @description
#' Internal helper for parallel execution with proper handling of
#' split indices to prevent data leakage in cached computations.
#'
#' @param items List of items to process
#' @param fun Function to apply to each item
#' @param split_id Current CV split identifier (for cache keys)
#' @param ... Additional arguments passed to fun
#'
#' @return List of results
#'
#' @keywords internal
parallel_lapply <- function(items, fun, split_id = NULL, ...) {

  if (!requireNamespace("future.apply", quietly = TRUE)) {
    # Fallback to sequential
    return(lapply(items, fun, ...))
  }

  # Add split_id to environment for cache key generation
  if (!is.null(split_id)) {
    # This ensures any caching inside `fun` can access the split context
    environment(fun)$.omic_split_id <- split_id
  }

  future.apply::future_lapply(items, fun, ..., future.seed = TRUE)
}


#' @title With Parallel Scope
#'
#' @description
#' Execute code within a parallelization scope, automatically cleaning up.
#'
#' @param expr Expression to evaluate
#' @param workers Number of workers
#' @param plan Parallelization plan
#'
#' @return Result of expression
#'
#' @export
with_parallel <- function(expr, workers = NULL, plan = "multisession") {
  old_plan <- setup_parallel(workers = workers, plan = plan)
  on.exit(future::plan(old_plan), add = TRUE)
  force(expr)
}
