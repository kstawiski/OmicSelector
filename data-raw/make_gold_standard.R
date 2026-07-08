#' Generate Gold Standard Synthetic Datasets
#'
#' This script generates the OmicSelector Gold Standard synthetic datasets
#' for use in package data and testing. The datasets are designed to detect
#' data leakage in feature selection and cross-validation pipelines.
#'
#' @note This script assumes the package is loaded or sourced.
#' Run with: Rscript data-raw/make_gold_standard.R
#' Or interactively after devtools::load_all()

# Load the package (or source files during development)
if (!requireNamespace("OmicSelector", quietly = TRUE)) {
  # During development, use devtools
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".")
  } else {
    stop("Please install OmicSelector or run devtools::load_all() first.")
  }
} else {
  library(OmicSelector)
}

# =========================================================================
# GENERATE AND SAVE DATASETS
# =========================================================================

if (interactive() || !exists("SKIP_DATA_GENERATION")) {
  cat("Generating Gold Standard datasets...\n")

  # Small dataset for quick CI tests
  gold_standard_small <- generate_gold_standard(
    n_patients = 50,
    samples_per_patient = 2,
    n_features = 200,
    n_causal = 10,
    n_trap = 15,
    seed = 42
  )

  # Full dataset for thorough testing
  gold_standard_full <- generate_gold_standard(
    n_patients = 100,
    samples_per_patient = 3,
    n_features = 1000,
    n_causal = 20,
    n_trap = 30,
    seed = 42
  )

  cat("Saving datasets to data/...\n")

  save(gold_standard_small, file = "data/gold_standard_small.rda", compress = "xz")
  save(gold_standard_full, file = "data/gold_standard_full.rda", compress = "xz")

  # Note: Functions are now exported from the package, no need to save them
  # The old gold_standard_functions.rda is deprecated

  cat("Done!\n")
  print(gold_standard_small)
}
