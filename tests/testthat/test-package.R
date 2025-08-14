test_that("Package loads successfully", {
  expect_true(requireNamespace("OmicSelector", quietly = TRUE))
})

test_that("Main functions are available", {
  expect_true(exists("OmicSelector_OmicSelector"))
  expect_true(exists("OmicSelector_benchmark"))
  expect_true(exists("OmicSelector_differential_expression_ttest"))
})

test_that("Package has required dependencies", {
  # Core dependencies should be available
  required_packages <- c("dplyr", "ggplot2", "caret", "randomForest")
  
  for (pkg in required_packages) {
    expect_true(
      requireNamespace(pkg, quietly = TRUE),
      info = paste("Required package", pkg, "is not available")
    )
  }
})
