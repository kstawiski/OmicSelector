# Run testthat tests for OmicSelector package

if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(OmicSelector)

  test_check("OmicSelector")
}
