#!/usr/bin/env Rscript

# Fast installation of minimal dependencies for testing
cat("Installing minimal test dependencies...\n\n")

options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Only install what we absolutely need for basic tests
minimal_deps <- c("dplyr")

for (pkg in minimal_deps) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, dependencies = FALSE, quiet = TRUE)
  }
}

cat("\n✓ Minimal dependencies installed\n")
