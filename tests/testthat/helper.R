# Test helper functions
library(testthat)

# Create sample test data
create_test_data <- function(n_samples = 100, n_features = 50, n_classes = 2) {
  set.seed(123)
  
  # Create feature matrix
  features <- matrix(
    rnorm(n_samples * n_features), 
    nrow = n_samples, 
    ncol = n_features
  )
  
  # Create meaningful feature names
  colnames(features) <- paste0("feature_", sprintf("%03d", 1:n_features))
  
  # Create class labels
  classes <- factor(sample(0:(n_classes-1), n_samples, replace = TRUE))
  
  # Combine into data frame
  data <- data.frame(
    class = classes,
    features,
    stringsAsFactors = FALSE
  )
  
  return(data)
}

# Create test configuration
create_test_config <- function() {
  list(
    max_iterations = 2,
    timeout_sec = 30,
    cores = 1,
    prefer_no_features = 5
  )
}

# Skip tests if optional dependencies are not available
skip_if_no_keras <- function() {
  testthat::skip_if_not_installed("keras")
  testthat::skip_if_not_installed("tensorflow")
  testthat::skip_if_not_installed("reticulate")
}

skip_if_no_bioc <- function() {
  testthat::skip_if_not_installed("edgeR")
  testthat::skip_if_not_installed("limma")
  testthat::skip_if_not_installed("Biobase")
}

skip_if_no_tidymodels <- function() {
  testthat::skip_if_not_installed("tidymodels")
  testthat::skip_if_not_installed("recipes")
  testthat::skip_if_not_installed("workflows")
}
