#!/usr/bin/env Rscript

# Script to regenerate NAMESPACE and documentation
# Run this before testing the package

cat("Regenerating documentation with roxygen2...\n")

# Load required package
if (!requireNamespace("roxygen2", quietly = TRUE)) {
  stop("roxygen2 is required. Install with: install.packages('roxygen2')")
}

# Set working directory to package root
setwd("/home/user/OmicSelector")

# Roxygenize the package
roxygen2::roxygenise()

cat("Documentation regenerated successfully!\n")
cat("NAMESPACE file updated.\n")
cat("man/ files updated.\n")
