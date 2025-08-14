test_that("Data loading functions work", {
  
  # Create temporary test data
  temp_dir <- tempdir()
  old_wd <- getwd()
  setwd(temp_dir)
  
  # Create mock CSV files
  train_data <- data.frame(
    class = factor(c("A", "B", "A", "B", "A", "B")),
    feature1 = rnorm(6),
    feature2 = rnorm(6),
    feature3 = rnorm(6)
  )
  
  test_data <- data.frame(
    class = factor(c("A", "B", "A", "B")),
    feature1 = rnorm(4),
    feature2 = rnorm(4),
    feature3 = rnorm(4)
  )
  
  write.csv(train_data, "mixed_train.csv", row.names = FALSE)
  write.csv(test_data, "mixed_test.csv", row.names = FALSE)
  
  # Test data loading
  expect_silent(data_list <- load_omics_data(temp_dir))
  expect_type(data_list, "list")
  expect_true("train" %in% names(data_list))
  expect_true("test" %in% names(data_list))
  
  # Test data validation
  expect_equal(nrow(data_list$train), 6)
  expect_equal(ncol(data_list$train), 4)
  
  # Clean up
  setwd(old_wd)
  unlink(file.path(temp_dir, "mixed_train.csv"))
  unlink(file.path(temp_dir, "mixed_test.csv"))
})

test_that("Method name mapping works", {
  
  # Test known methods
  expect_equal(get_method_name(1), "Sig (t-test p<0.05)")
  expect_equal(get_method_name(17), "Boruta")
  
  # Test unknown method
  expect_equal(get_method_name(99), "Method 99")
})

test_that("Individual feature selection methods work", {
  
  # Create test data
  test_data <- data.frame(
    class = factor(rep(c("A", "B"), each = 10)),
    feature1 = c(rnorm(10, 0), rnorm(10, 2)),  # Should be significant
    feature2 = rnorm(20),                       # Should not be significant
    feature3 = c(rnorm(10, 0), rnorm(10, 1.5)) # Moderately significant
  )
  
  config <- get_config("development")
  
  # Test Sig method
  result_sig <- method_sig(test_data, config)
  expect_type(result_sig, "list")
  expect_true("formula" %in% names(result_sig))
  expect_true("features" %in% names(result_sig))
  
  if (!is.null(result_sig$features)) {
    expect_true(length(result_sig$features) > 0)
    expect_true(all(result_sig$features %in% colnames(test_data)[-1]))
  }
  
  # Test CFS method (simplified version)
  result_cfs <- method_cfs(test_data, config)
  expect_type(result_cfs, "list")
  
  # Test Random Forest RFE (if randomForest is available)
  skip_if_not_installed("randomForest")
  result_rf <- method_rf_rfe(test_data, config)
  expect_type(result_rf, "list")
})

test_that("Feature selection pipeline integration", {
  
  # Create temporary directory with test data
  temp_dir <- tempdir()
  old_wd <- getwd()
  
  # Create comprehensive test data
  set.seed(123)
  n_samples <- 50
  n_features <- 20
  
  # Create data with some meaningful signal
  train_data <- create_test_data(n_samples, n_features, 2)
  test_data <- create_test_data(20, n_features, 2)
  
  # Save to temporary directory
  write.csv(train_data, file.path(temp_dir, "mixed_train.csv"), row.names = FALSE)
  write.csv(test_data, file.path(temp_dir, "mixed_test.csv"), row.names = FALSE)
  
  # Test pipeline with development config
  config <- get_config("development")
  methods <- c(1, 3)  # Sig and CFS methods
  
  # Initialize results object
  data_list <- list(train = train_data, test = test_data)
  results <- initialize_results_object("test_pipeline", config, methods, data_list)
  
  # Test pipeline execution
  expect_silent(
    results_updated <- run_feature_selection_pipeline(results, methods, config)
  )
  
  expect_type(results_updated, "list")
  expect_s3_class(results_updated, "OmicSelector")
  
  # Clean up
  setwd(old_wd)
  unlink(file.path(temp_dir, "mixed_train.csv"))
  unlink(file.path(temp_dir, "mixed_test.csv"))
})
