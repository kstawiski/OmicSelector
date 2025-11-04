#!/usr/bin/env Rscript

# Install essential packages for testing OmicSelector 2.0
cat("Installing required packages for OmicSelector 2.0 testing...\n")
cat("This may take several minutes...\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Essential packages for testing
essential_packages <- c(
  "devtools",
  "roxygen2",
  "testthat",
  "dplyr",
  "ggplot2"
)

cat("Phase 1: Installing essential packages...\n")
for (pkg in essential_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing", pkg, "...\n")
    install.packages(pkg, quiet = TRUE)
  } else {
    cat("  ✓", pkg, "already installed\n")
  }
}

cat("\nPhase 2: Checking installation...\n")
success <- TRUE
for (pkg in essential_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("  ✓", pkg, "\n")
  } else {
    cat("  ✗", pkg, "FAILED\n")
    success <- FALSE
  }
}

if (success) {
  cat("\n✓ All essential packages installed successfully!\n")
} else {
  cat("\n✗ Some packages failed to install\n")
  quit(status = 1)
}
