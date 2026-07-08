#' @title Split-Aware Caching for OmicSelector
#'
#' @description
#' Caching utilities that prevent data leakage by including split identifiers
#' in cache keys. Uses memoise with cachem for efficient LRU caching.
#'
#' @details
#' CRITICAL for zero-leakage: Cache keys MUST include:
#' - Split/fold identifier (prevents cross-fold contamination)
#' - Preprocessing parameters (ensures correct pipeline state)
#' - Feature subset hash (for filter-specific caching)
#'
#' This module ensures that cached filter scores or model artifacts
#' from one CV fold are never reused in another fold.
#'
#' @name cache
NULL


#' @title Create Split-Aware Cache
#'
#' @description
#' Creates a cache instance configured for OmicSelector's zero-leakage requirements.
#' Uses an LRU (Least Recently Used) eviction policy.
#'
#' @param max_size Maximum cache size in bytes. Default: 200MB.
#' @param max_age Maximum age of cached items in seconds. Default: 3600 (1 hour).
#' @param dir Optional directory for disk-based caching. NULL for memory-only.
#'
#' @return A cachem cache object
#'
#' @export
create_omic_cache <- function(max_size = 200 * 1024^2,
                               max_age = 3600,
                               dir = NULL) {

  if (!requireNamespace("cachem", quietly = TRUE)) {
    warning("cachem not installed. Caching disabled.")
    return(NULL)
  }

  if (is.null(dir)) {
    # Memory-only cache
    cache <- cachem::cache_mem(
      max_size = max_size,
      max_age = max_age
    )
  } else {
    # Disk-backed cache for persistence
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    cache <- cachem::cache_disk(
      dir = dir,
      max_size = max_size,
      max_age = max_age
    )
  }

  # Add metadata
  attr(cache, "omic_cache") <- TRUE
  attr(cache, "created") <- Sys.time()

  return(cache)
}


#' @title Generate Split-Aware Cache Key
#'
#' @description
#' Creates a cache key that includes split information to prevent leakage.
#' This is the CRITICAL function for safe caching in nested CV.
#'
#' @param operation Name of the operation being cached (e.g., "filter_anova")
#' @param split_id Identifier for the current CV split (REQUIRED)
#' @param features Character vector of features involved
#' @param params Named list of parameters affecting the result
#'
#' @return A character string cache key
#'
#' @details
#' The key structure: "operation::split_id::features_hash::params_hash"
#'
#' Example: "filter_anova::outer3_inner2::a1b2c3d4::x5y6z7"
#'
#' @export
generate_cache_key <- function(operation,
                                split_id,
                                features = NULL,
                                params = NULL) {

  if (is.null(split_id) || split_id == "") {
    stop("split_id is REQUIRED for safe caching. Cannot cache without CV fold context.",
         call. = FALSE)
  }

  key_parts <- c(operation, split_id)

  # Hash features if provided

if (!is.null(features) && length(features) > 0) {
    features_sorted <- sort(features)
    features_hash <- substr(
      digest::digest(features_sorted, algo = "xxhash32"),
      1, 8
    )
    key_parts <- c(key_parts, features_hash)
  }

  # Hash params if provided
  if (!is.null(params) && length(params) > 0) {
    params_hash <- substr(
      digest::digest(params, algo = "xxhash32"),
      1, 8
    )
    key_parts <- c(key_parts, params_hash)
  }

  paste(key_parts, collapse = "::")
}


#' @title Cached Filter Computation
#'
#' @description
#' Wrapper for filter score computation with split-aware caching.
#' Automatically handles cache key generation and invalidation.
#'
#' @param filter_fun Function that computes filter scores
#' @param task mlr3 Task object
#' @param filter_name Name of the filter (e.g., "anova", "mrmr")
#' @param split_id Current CV split identifier (REQUIRED)
#' @param cache Cache object from create_omic_cache()
#'
#' @return Filter scores (from cache or freshly computed)
#'
#' @export
cached_filter <- function(filter_fun,
                           task,
                           filter_name,
                           split_id,
                           cache = NULL) {

  if (is.null(cache)) {
    # No cache, compute directly
    return(filter_fun(task))
  }

  # Generate cache key
  key <- generate_cache_key(
    operation = paste0("filter_", filter_name),
    split_id = split_id,
    features = task$feature_names
  )

  # Check cache
  cached_value <- cache$get(key)
  if (!is.null(cached_value)) {
    return(cached_value)
  }

  # Compute and cache
  result <- filter_fun(task)
  cache$set(key, result)

  return(result)
}


#' @title Memoize Function with Split Context
#'
#' @description
#' Creates a memoized version of a function that includes split_id in the cache key.
#' This is the recommended way to add caching to custom functions.
#'
#' @param f Function to memoize
#' @param cache Cache object (optional, creates new if NULL)
#'
#' @return A memoized function that requires split_id as first argument
#'
#' @details
#' The returned function has signature: function(split_id, ...)
#' The split_id is used in the cache key but passed through to f.
#'
#' @examples
#' \dontrun{
#' # Original function
#' expensive_computation <- function(data, param) {
#'   # ... expensive operation ...
#' }
#'
#' # Memoized version
#' cached_computation <- memoize_with_split(expensive_computation)
#'
#' # Usage in CV loop
#' for (i in 1:5) {
#'   split_id <- paste0("fold_", i)
#'   result <- cached_computation(split_id, data = fold_data, param = 0.5)
#' }
#' }
#'
#' @export
memoize_with_split <- function(f, cache = NULL) {

  if (!requireNamespace("memoise", quietly = TRUE)) {
    warning("memoise not installed. Returning unmemoized function.")
    return(function(split_id, ...) f(...))
  }

  # Create cache if not provided
  if (is.null(cache)) {
    cache <- create_omic_cache()
  }

  # Create wrapper that incorporates split_id into memoization
  function(split_id, ...) {
    # The split_id becomes part of the cache key through memoise
    if (is.null(cache)) {
      return(f(...))
    }

    key <- generate_cache_key(
      operation = deparse(substitute(f)),
      split_id = split_id,
      params = list(...)
    )

    cached <- cache$get(key)
    if (!is.null(cached)) {
      return(cached)
    }

    result <- f(...)
    cache$set(key, result)
    return(result)
  }
}


#' @title Clear OmicSelector Cache
#'
#' @description
#' Clears all cached values. Use between experiments or when data changes.
#'
#' @param cache Cache object to clear
#'
#' @export
clear_cache <- function(cache) {
  if (is.null(cache)) {
    return(invisible(NULL))
  }

  if (inherits(cache, "cache_mem") || inherits(cache, "cache_disk")) {
    cache$reset()
    message("Cache cleared")
  }

  invisible(NULL)
}


#' @title Get Cache Statistics
#'
#' @description
#' Returns information about cache usage.
#'
#' @param cache Cache object
#'
#' @return A list with cache statistics
#'
#' @export
cache_stats <- function(cache) {
  if (is.null(cache)) {
    return(list(enabled = FALSE))
  }

  list(
    enabled = TRUE,
    type = class(cache)[1],
    created = attr(cache, "created"),
    info = tryCatch(
      cache$info(),
      error = function(e) list(error = e$message)
    )
  )
}
