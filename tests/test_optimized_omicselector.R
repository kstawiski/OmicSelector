# Test for Optimized OmicSelector Functions
# ===========================================

# Load required libraries
suppressMessages({
  library(testthat)
  library(dplyr)
})

# Create test data for miRNA Case-Control analysis
create_test_mirna_data <- function(n_samples = 100, n_features = 50, seed = 42) {
  set.seed(seed)
  
  # Create miRNA feature names with hsa_ prefix
  feature_names <- paste0("hsa_miR_", sample(1:1000, n_features))
  
  # Create class labels (Case/Control)
  n_case <- floor(n_samples * 0.6)
  n_control <- n_samples - n_case
  class_labels <- c(rep("Case", n_case), rep("Control", n_control))
  
  # Create feature matrix with some differential expression
  feature_data <- matrix(rnorm(n_samples * n_features), nrow = n_samples, ncol = n_features)
  colnames(feature_data) <- feature_names
  
  # Make some features more expressed in cases
  case_indices <- which(class_labels == "Case")
  n_diff <- min(10, n_features)  # 10 differentially expressed features
  diff_features <- 1:n_diff
  
  for (i in diff_features) {
    feature_data[case_indices, i] <- feature_data[case_indices, i] + rnorm(length(case_indices), mean = 2, sd = 0.5)
  }
  
  # Combine into data frame
  data <- data.frame(
    Class = factor(class_labels, levels = c("Control", "Case")),
    feature_data,
    stringsAsFactors = FALSE
  )
  
  return(data)
}

# Create test data files
setup_test_environment <- function(test_dir = tempdir()) {
  # Create test directory
  if (!dir.exists(test_dir)) {
    dir.create(test_dir, recursive = TRUE)
  }
  
  # Generate test datasets
  set.seed(42)
  train_data <- create_test_mirna_data(n_samples = 120, n_features = 30)
  test_data <- create_test_mirna_data(n_samples = 40, n_features = 30)
  validation_data <- create_test_mirna_data(n_samples = 40, n_features = 30)
  
  # Save as CSV files
  write.csv(train_data, file.path(test_dir, "mixed_train.csv"), row.names = FALSE)
  write.csv(test_data, file.path(test_dir, "mixed_test.csv"), row.names = FALSE)
  write.csv(validation_data, file.path(test_dir, "mixed_validation.csv"), row.names = FALSE)
  
  return(test_dir)
}

# Test the optimized helper functions
test_helper_functions <- function() {
  cat("Testing helper functions...\n")
  
  # Test data creation
  test_data <- create_test_mirna_data(n_samples = 50, n_features = 20)
  
  # Test get_feature_columns
  hsa_features <- get_feature_columns(test_data, "hsa_")
  stopifnot(length(hsa_features) == 20)
  stopifnot(all(startsWith(hsa_features, "hsa_")))
  
  # Test validate_mirna_features
  valid_features <- validate_mirna_features(hsa_features, "hsa_")
  stopifnot(all(valid_features))  # All should be valid
  
  # Test invalid features
  invalid_features <- c("hsa_miR_123", "invalid_feature", "hsa_let_456")
  valid_check <- validate_mirna_features(invalid_features, "hsa_")
  stopifnot(sum(valid_check) == 2)  # Only 2 should be valid
  
  # Test fast_ttest
  feature_matrix <- as.matrix(test_data[, hsa_features])
  de_results <- fast_ttest(feature_matrix, test_data$Class, "Case", "Control")
  stopifnot(nrow(de_results) == length(hsa_features))
  stopifnot(all(c("feature", "mean_case", "mean_control", "log2fc", "p_value") %in% colnames(de_results)))
  
  # Test create_feature_formula
  selected_features <- hsa_features[1:5]
  formula_obj <- create_feature_formula(selected_features, "Class")
  stopifnot(class(formula_obj) == "formula")
  
  cat("✓ Helper functions tests passed\n")
}

# Test the optimized main function
test_optimized_function <- function() {
  cat("Testing optimized OmicSelector function...\n")
  
  # Setup test environment
  test_dir <- setup_test_environment()
  
  # Test with minimal methods for speed
  test_methods <- c(1, 2, 3)
  
  tryCatch({
    # Run optimized function
    results <- omicselector_optimized(
      data_path = test_dir,
      methods = test_methods,
      max_features = 10,
      max_iterations = 3,
      timeout_minutes = 5,
      parallel_cores = 1,  # Sequential for testing
      verbose = FALSE,
      seed = 42
    )
    
    # Validate results structure
    stopifnot(inherits(results, "omicselector_results"))
    stopifnot(!is.null(results$selected_features))
    stopifnot(!is.null(results$differential_expression))
    stopifnot(!is.null(results$metadata))
    stopifnot(!is.null(results$summary))
    
    # Check that methods were executed
    expected_methods <- paste0("method_", test_methods)
    stopifnot(all(expected_methods %in% names(results$selected_features)))
    
    # Check differential expression results
    de_results <- results$differential_expression
    stopifnot(nrow(de_results) > 0)
    stopifnot(all(c("feature", "p_value", "log2fc") %in% colnames(de_results)))
    
    # Check that some features were selected
    total_features <- sum(sapply(results$selected_features, length))
    stopifnot(total_features > 0)
    
    cat("✓ Optimized function tests passed\n")
    
  }, error = function(e) {
    cat("✗ Optimized function test failed:", e$message, "\n")
    stop(e)
  })
  
  # Cleanup
  unlink(test_dir, recursive = TRUE)
}

# Test the improved OmicSelector_OmicSelector wrapper
test_improved_wrapper <- function() {
  cat("Testing improved OmicSelector_OmicSelector wrapper...\n")
  
  # Setup test environment
  test_dir <- setup_test_environment()
  
  tryCatch({
    # Test the improved wrapper
    results <- OmicSelector_OmicSelector(
      wd = test_dir,
      m = c(1, 2),
      max_iterations = 3,
      prefer_no_features = 10,
      register_parallel = FALSE
    )
    
    # Check backward compatibility - should have formulas
    stopifnot(!is.null(results$formulas))
    stopifnot(length(results$formulas) > 0)
    
    # Check that formulas are actual formula objects
    for (formula_obj in results$formulas) {
      stopifnot(inherits(formula_obj, "formula"))
    }
    
    # Check new features are also available
    stopifnot(!is.null(results$selected_features))
    stopifnot(!is.null(results$metadata))
    
    cat("✓ Improved wrapper tests passed\n")
    
  }, error = function(e) {
    cat("✗ Improved wrapper test failed:", e$message, "\n")
    stop(e)
  })
  
  # Cleanup
  unlink(test_dir, recursive = TRUE)
}

# Test performance comparison
test_performance <- function() {
  cat("Testing performance improvements...\n")
  
  # Setup test environment
  test_dir <- setup_test_environment()
  
  # Test with minimal methods but measure time
  start_time <- Sys.time()
  
  perf_results <- omicselector_optimized(
    data_path = test_dir,
    methods = c(1, 2, 3, 4),  # A few methods
    max_features = 15,
    max_iterations = 5,
    timeout_minutes = 10,
    parallel_cores = 1,
    verbose = FALSE,
    seed = 42
  )
  
  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat("✓ Performance test completed in", round(runtime, 2), "seconds\n")
  
  # Check that it completed in reasonable time (should be fast for small test)
  stopifnot(runtime < 60)  # Should complete in under 1 minute for test data
  
  # Validate results were actually generated
  stopifnot(!is.null(perf_results$selected_features))
  
  # Cleanup
  unlink(test_dir, recursive = TRUE)
}

# Test error handling
test_error_handling <- function() {
  cat("Testing error handling...\n")
  
  # Test with invalid data path
  tryCatch({
    results <- omicselector_optimized(
      data_path = "/nonexistent/path",
      methods = c(1, 2),
      verbose = FALSE
    )
    stop("Should have failed with invalid path")
  }, error = function(e) {
    # Expected to fail
    cat("✓ Correctly caught invalid path error\n")
  })
  
  # Test with invalid methods
  test_dir <- setup_test_environment()
  
  tryCatch({
    error_results <- omicselector_optimized(
      data_path = test_dir,
      methods = c(999),  # Invalid method
      verbose = FALSE
    )
    # Use error_results to avoid unused variable warning
    if (!is.null(error_results)) {
      stop("Should have failed with invalid methods")
    }
  }, error = function(e) {
    # Expected to fail
    cat("✓ Correctly caught invalid methods error\n")
  })
  
  # Cleanup
  unlink(test_dir, recursive = TRUE)
}

# Run all tests
run_all_tests <- function() {
  cat("Running OmicSelector Optimization Tests\n")
  cat("======================================\n\n")
  
  start_time <- Sys.time()
  
  # Run individual test suites
  test_helper_functions()
  test_optimized_function()
  test_improved_wrapper()
  test_performance()
  test_error_handling()
  
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat("\n======================================\n")
  cat("All tests completed successfully!\n")
  cat("Total test time:", round(total_time, 2), "seconds\n")
  cat("======================================\n")
  
  return(TRUE)
}

# Example usage function
demo_optimized_omicselector <- function() {
  cat("OmicSelector Optimization Demo\n")
  cat("=============================\n\n")
  
  # Create demo data
  test_dir <- setup_test_environment()
  cat("Created demo miRNA dataset in:", test_dir, "\n")
  
  # Show data structure
  train_data <- read.csv(file.path(test_dir, "mixed_train.csv"))
  cat("Training data:", nrow(train_data), "samples,", ncol(train_data) - 1, "features\n")
  cat("Class distribution:", table(train_data$Class), "\n")
  cat("miRNA features:", sum(startsWith(colnames(train_data), "hsa_")), "\n\n")
  
  # Run analysis
  cat("Running optimized OmicSelector analysis...\n")
  results <- OmicSelector_OmicSelector(
    wd = test_dir,
    m = c(1, 2, 3, 4, 11),  # Selection of methods
    max_iterations = 5,
    prefer_no_features = 10,
    register_parallel = FALSE
  )
  
  # Show results
  cat("\nResults summary:\n")
  print(results)
  
  # Show selected features
  cat("\nSelected features by method:\n")
  for (method_name in names(results$selected_features)) {
    features <- results$selected_features[[method_name]]
    cat(method_name, ":", length(features), "features\n")
    if (length(features) > 0) {
      cat("  Top features:", paste(head(features, 3), collapse = ", "), "\n")
    }
  }
  
  # Cleanup
  unlink(test_dir, recursive = TRUE)
  
  cat("\nDemo completed successfully!\n")
  return(results)
}

# Run tests if script is executed directly
if (!interactive()) {
  run_all_tests()
}
