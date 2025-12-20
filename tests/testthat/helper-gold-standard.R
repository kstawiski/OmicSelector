# Helper functions for testthat tests
# This file is automatically sourced by testthat before running tests

# Set SKIP_DATA_GENERATION to prevent auto-generation when sourcing
SKIP_DATA_GENERATION <- TRUE

# Try multiple paths to find the Gold Standard dataset generator
possible_paths <- c(
  # When running from package root
  "data-raw/make_gold_standard.R",
  # When running from tests/testthat/
  "../../data-raw/make_gold_standard.R",
  # Absolute path fallback
  file.path(getwd(), "data-raw", "make_gold_standard.R")
)

# Also try to find package root via rprojroot if available
if (requireNamespace("rprojroot", quietly = TRUE)) {
  tryCatch({
    pkg_root <- rprojroot::find_package_root_file()
    possible_paths <- c(
      file.path(pkg_root, "data-raw", "make_gold_standard.R"),
      possible_paths
    )
  }, error = function(e) NULL)
}

gold_standard_loaded <- FALSE
for (path in possible_paths) {
  if (file.exists(path)) {
    tryCatch({
      source(path, local = TRUE)
      gold_standard_loaded <- TRUE
      message("Gold Standard helper loaded from: ", path)
      break
    }, error = function(e) {
      warning("Failed to source ", path, ": ", e$message)
    })
  }
}

# Only warn, don't stop - some tests don't need gold standard
if (!gold_standard_loaded) {
  warning("Gold Standard generator not found. Related tests will be skipped.")
}

# Helper to check if gold standard functions are available
gold_standard_available <- function() {
  exists("generate_gold_standard") && exists("evaluate_feature_selection")
}
