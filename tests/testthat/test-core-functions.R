test_that("Configuration management works", {
  
  # Test default config loading
  config <- get_config()
  expect_type(config, "list")
  expect_true("max_iterations" %in% names(config))
  expect_true("cores" %in% names(config))
  
  # Test config override
  config_custom <- get_config(max_iterations = 5, cores = 2)
  expect_equal(config_custom$max_iterations, 5)
  expect_equal(config_custom$cores, 2)
  
  # Test config validation
  expect_error(validate_config(list(max_iterations = -1)))
  expect_error(validate_config(list(timeout_sec = 0)))
  
  validated <- validate_config(list(
    max_iterations = 10,
    timeout_sec = 3600,
    cores = 2
  ))
  expect_equal(validated$max_iterations, 10)
})

test_that("Logging system works", {
  
  # Test basic logging functions exist
  expect_true(exists("log_info"))
  expect_true(exists("log_warn"))
  expect_true(exists("log_error"))
  expect_true(exists("log_debug"))
  
  # Test logging setup
  expect_silent(setup_logging("INFO"))
  
  # Test progress bar creation
  pb <- create_progress_bar(10, "Test")
  expect_true(!is.null(pb))
  
  # Clean up
  if (requireNamespace("cli", quietly = TRUE) && !is.list(pb)) {
    finish_progress_bar(pb)
  }
})

test_that("Validation functions work", {
  
  # Test method validation
  expect_equal(validate_methods(c(1, 2, 3)), c(1, 2, 3))
  expect_error(validate_methods(c(0, 1, 2)))
  expect_error(validate_methods(c(1, 2, 100)))
  expect_warning(validate_methods(c(1, 1, 2)))
  
  # Test data frame validation
  test_data <- data.frame(
    class = factor(c("A", "B", "A", "B")),
    feature1 = c(1, 2, 3, 4),
    feature2 = c(5, 6, 7, 8)
  )
  
  expect_true(validate_data_frame(test_data))
  expect_error(validate_data_frame(test_data, min_rows = 10))
  expect_error(validate_data_frame(test_data, required_cols = "missing_col"))
  
  # Test parallel config validation
  parallel_config <- validate_parallel_config(2, TRUE)
  expect_type(parallel_config, "list")
  expect_true("cores" %in% names(parallel_config))
})

test_that("Dependency checking works", {
  
  # Test with packages that should be available
  available <- check_optional_dependencies("stats", action = "skip")
  expect_true(available["stats"])
  
  # Test with non-existent package
  expect_warning(
    check_optional_dependencies("nonexistent_package_12345", action = "warn")
  )
})
