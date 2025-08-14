#!/usr/bin/env Rscript

# Simple test of modular framework
source('R/modular-fs-framework.R')
source('R/core-fs-methods.R')
source('R/modular-integration.R')

# Create test data
test_data <- data.frame(
  class = factor(rep(c("A", "B"), each = 30)),
  feature1 = rnorm(60, mean = rep(c(0, 1), each = 30)),
  feature2 = rnorm(60),
  feature3 = rnorm(60, mean = rep(c(1, 0), each = 30)),
  feature4 = rnorm(60),
  feature5 = rnorm(60, mean = rep(c(0.5, -0.5), each = 30))
)

cat("Testing modular framework API...\n")
results <- modular_feature_selection(
  data = test_data,
  methods = c(1, 2),
  max_features = 3,
  parallel = FALSE,
  progress = TRUE
)

cat("Success! Got results for", length(results), "methods\n")
for (i in seq_along(results)) {
  result <- results[[i]]
  if (!is.null(result)) {
    cat("Method", i, ":", length(result$features), "features:", 
        paste(result$features, collapse = ", "), "\n")
  } else {
    cat("Method", i, ": failed\n")
  }
}

cat("\nModular framework test completed successfully!\n")
