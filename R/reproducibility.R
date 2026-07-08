#' @title Internal Reproducibility Helpers
#' @description
#' Internal helpers for consistent seed and thread handling.
#' @keywords internal
set_reproducible <- function(seed = NULL, threads = 1) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.null(threads) &&
      requireNamespace("mlr3", quietly = TRUE)) {
    mlr3::set_threads(threads)
  }

  invisible(TRUE)
}

#' @title Build a Stable Cache Key
#' @description
#' Internal helper for cache key hashing.
#' @keywords internal
hash_cache_key <- function(key) {
  tmp <- tempfile(pattern = "omicselector_key_", fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(as.character(key), con = tmp)
  tools::md5sum(tmp)
}
