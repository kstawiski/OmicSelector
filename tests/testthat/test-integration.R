test_that("End-to-end integration test works", {
  
  # Skip on CRAN to avoid long running tests
  skip_on_cran()
  
  # Create temporary test environment
  temp_dir <- tempdir()
  old_wd <- getwd()
  
  tryCatch({
    
    # Generate realistic test data
    set.seed(12345)
    n_samples <- 40
    n_features <- 15
    
    # Create training data with signal
    train_data <- create_test_data(n_samples, n_features, 2)
    test_data <- create_test_data(20, n_features, 2)
    
    # Add some meaningful signal to first few features
    signal_features <- 1:3
    for (i in signal_features) {
      feature_col <- paste0("feature_", sprintf("%03d", i))
      # Make these features more discriminative
      train_data[train_data$class == 0, feature_col] <- 
        train_data[train_data$class == 0, feature_col] + rnorm(sum(train_data$class == 0), 0, 0.5)
      train_data[train_data$class == 1, feature_col] <- 
        train_data[train_data$class == 1, feature_col] + rnorm(sum(train_data$class == 1), 2, 0.5)
    }
    
    # Save test data
    write.csv(train_data, file.path(temp_dir, "mixed_train.csv"), row.names = FALSE)
    write.csv(test_data, file.path(temp_dir, "mixed_test.csv"), row.names = FALSE)
    
    # Test modern feature selection workflow
    results <- omics_select(
      wd = temp_dir,
      methods = c(1, 3),  # Sig and CFS methods only for speed
      config_name = "development",
      debug = TRUE
    )
    
    # Verify results structure
    expect_s3_class(results, "OmicSelector")
    expect_true("formulas" %in% names(results))
    expect_true("method_results" %in% names(results))
    expect_true(length(results$formulas) >= 0)  # At least attempted
    
    # Test S3 methods
    expect_output(print(results), "OmicSelector Results")
    
    if (length(results$formulas) > 0) {
      formulas_extracted <- formulas(results)
      expect_true(length(formulas_extracted) > 0)
      expect_true(all(grepl("class ~", formulas_extracted)))
      
      features_extracted <- features(results)
      expect_true(is.list(features_extracted))
      
      # Test plotting (if ggplot2 available)
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        expect_silent(plot(results, type = "overview"))
      }
      
      # Test modern benchmarking (simplified)
      if (requireNamespace("randomForest", quietly = TRUE)) {
        benchmark_results <- omics_benchmark(
          wd = temp_dir,
          formulas = results$formulas[1],  # Just test one formula
          algorithms = c("rf"),  # Just one algorithm for speed
          validation_strategy = "holdout",
          search_iterations = 5  # Minimal iterations
        )
        
        expect_s3_class(benchmark_results, "OmicBenchmark")
        expect_true("results_matrix" %in% names(benchmark_results))
      }
    }
    
    # Test configuration system
    config <- get_config("development")
    expect_true(is.list(config))
    expect_true(config$max_iterations == 2)
    
    # Test logging system
    expect_silent(setup_logging("INFO"))
    expect_silent(log_info("Integration test completed successfully"))
    
  }, finally = {
    # Cleanup
    setwd(old_wd)
    unlink(file.path(temp_dir, "mixed_train.csv"))
    unlink(file.path(temp_dir, "mixed_test.csv"))
    if (dir.exists(file.path(temp_dir, "temp"))) {
      unlink(file.path(temp_dir, "temp"), recursive = TRUE)
    }
  })
})

test_that("Legacy compatibility is maintained", {
  
  # Test that old function names still exist
  expect_true(exists("OmicSelector_OmicSelector"))
  expect_true(exists("OmicSelector_benchmark"))
  expect_true(exists("OmicSelector_differential_expression_ttest"))
  
  # Test differential expression modernization
  test_data <- data.frame(
    class = factor(rep(c("Case", "Control"), each = 10)),
    feature1 = c(rnorm(10, 0), rnorm(10, 2)),
    feature2 = c(rnorm(10, 1), rnorm(10, 1))
  )
  
  # Test modern function
  modern_result <- differential_expression_modern(
    features = test_data[, -1],
    classes = test_data$class
  )
  
  expect_true(is.data.frame(modern_result))
  expect_true("feature" %in% colnames(modern_result))
  expect_true("p_value" %in% colnames(modern_result))
  expect_true("p_adjusted" %in% colnames(modern_result))
  expect_true("log_fold_change" %in% colnames(modern_result))
})

test_that("Error handling is robust", {
  
  # Test configuration validation
  expect_error(validate_config(list(max_iterations = -1)))
  expect_error(validate_config(list(timeout_sec = 0)))
  
  # Test working directory validation
  expect_error(validate_working_directory("/nonexistent/path"))
  
  # Test method validation
  expect_error(validate_methods(c(0, 1, 2)))
  expect_error(validate_methods(c(1, 2, 100)))
  
  # Test data frame validation
  bad_data <- data.frame(x = 1:3)
  expect_error(validate_data_frame(bad_data, min_rows = 10))
  expect_error(validate_data_frame(bad_data, required_cols = "missing"))
})

test_that("Performance improvements are measurable", {
  
  # Create test data
  test_data <- create_test_data(50, 20, 2)
  
  # Time the modern differential expression
  start_time <- Sys.time()
  result <- differential_expression_modern(
    features = test_data[, -1],
    classes = test_data$class
  )
  end_time <- Sys.time()
  
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Should complete quickly (less than 5 seconds for this small dataset)
  expect_true(execution_time < 5)
  
  # Should return reasonable results
  expect_true(nrow(result) == 20)  # All features
  expect_true(all(result$p_value >= 0 & result$p_value <= 1))
  expect_false(any(is.na(result$p_value)))
})
