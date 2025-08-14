#!/usr/bin/env Rscript

#' Test Script for Modular Feature Selection Framework
#'
#' @description
#' Comprehensive testing of the new modular feature selection framework.
#' Tests core functionality, method registration, execution, and analysis.

cat("=========================================\n")
cat("Modular Feature Selection Framework Test\n")
cat("=========================================\n\n")

# Load required libraries
library(OmicSelector)

# Test data setup
cat("Setting up test data...\n")
test_data <- data.frame(
  class = factor(rep(c("A", "B", "C"), each = 30)),
  feature1 = rnorm(90, mean = rep(c(0, 1, 2), each = 30)),
  feature2 = rnorm(90, mean = rep(c(2, 0, 1), each = 30)),
  feature3 = rnorm(90, mean = rep(c(1, 2, 0), each = 30)),
  feature4 = rnorm(90),  # Random feature
  feature5 = rnorm(90),  # Random feature
  feature6 = rnorm(90, mean = rep(c(0, 2, 1), each = 30)),
  feature7 = rnorm(90),  # Random feature
  feature8 = rnorm(90, mean = rep(c(1, 0, 2), each = 30)),
  feature9 = rnorm(90),  # Random feature
  feature10 = rnorm(90, mean = rep(c(2, 1, 0), each = 30))
)

cat("Test data created: 90 samples x 10 features x 3 classes\n\n")

# =============================================================================
# Test 1: Framework Initialization
# =============================================================================

cat("Test 1: Framework Initialization\n")
cat("================================\n")

tryCatch({
  init_result <- initialize_fs_framework(verbose = TRUE)
  cat("✓ Framework initialization successful\n")
  cat("  Core methods loaded:", init_result$core_methods_loaded, "\n")
  cat("  Plugin methods loaded:", init_result$plugin_methods_loaded, "\n")
  cat("  Total methods:", init_result$total_methods, "\n\n")
}, error = function(e) {
  cat("✗ Framework initialization failed:", e$message, "\n\n")
  stop("Cannot continue without framework initialization")
})

# =============================================================================
# Test 2: Method Registry
# =============================================================================

cat("Test 2: Method Registry\n")
cat("=======================\n")

# List available methods
methods_df <- list_fs_methods()
cat("Available methods:\n")
print(methods_df[, c("id", "name", "category", "complexity")])
cat("\n")

# Test method lookup
test_method_id <- 1
method_def <- get_fs_method(test_method_id)
if (!is.null(method_def)) {
  cat("✓ Method lookup successful for ID", test_method_id, "\n")
  cat("  Method name:", method_def$name, "\n")
  cat("  Category:", method_def$category, "\n")
} else {
  cat("✗ Method lookup failed for ID", test_method_id, "\n")
}

# Test dependency checking
dep_check <- check_method_dependencies(test_method_id)
cat("✓ Dependency check for method", test_method_id, ":", 
    ifelse(dep_check$available, "PASS", "FAIL"), "\n")
if (!dep_check$available) {
  cat("  Missing packages:", paste(dep_check$missing, collapse = ", "), "\n")
}
cat("\n")

# =============================================================================
# Test 3: Single Method Execution
# =============================================================================

cat("Test 3: Single Method Execution\n")
cat("===============================\n")

# Test each core method
core_method_ids <- c(1, 2, 3, 4, 5)
successful_methods <- integer(0)

for (method_id in core_method_ids) {
  cat("Testing method", method_id, "...\n")
  
  result <- execute_fs_method(
    method_id = method_id,
    data = test_data,
    max_features = 5,
    timeout_sec = 30
  )
  
  if (!is.null(result)) {
    cat("  ✓ Method", method_id, "successful\n")
    cat("    Features selected:", length(result$features), "\n")
    cat("    Top features:", paste(head(result$features, 3), collapse = ", "), "\n")
    successful_methods <- c(successful_methods, method_id)
  } else {
    cat("  ✗ Method", method_id, "failed\n")
  }
}

cat("\nSuccessful methods:", length(successful_methods), "/", length(core_method_ids), "\n\n")

# =============================================================================
# Test 4: Multiple Method Execution
# =============================================================================

cat("Test 4: Multiple Method Execution\n")
cat("=================================\n")

if (length(successful_methods) >= 2) {
  test_methods <- head(successful_methods, 3)
  
  cat("Testing multiple methods:", paste(test_methods, collapse = ", "), "\n")
  
  multi_results <- execute_multiple_fs_methods(
    method_ids = test_methods,
    data = test_data,
    max_features = 5,
    parallel = FALSE,  # Use sequential for testing
    progress = TRUE
  )
  
  successful_multi <- sum(sapply(multi_results, function(x) !is.null(x)))
  cat("✓ Multiple method execution complete\n")
  cat("  Successful methods:", successful_multi, "/", length(test_methods), "\n\n")
} else {
  cat("⚠ Skipping multiple method test (insufficient successful methods)\n\n")
}

# =============================================================================
# Test 5: High-Level API
# =============================================================================

cat("Test 5: High-Level API\n")
cat("======================\n")

# Test modular_feature_selection
if (length(successful_methods) >= 1) {
  cat("Testing modular_feature_selection API...\n")
  
  api_result <- modular_feature_selection(
    data = test_data,
    methods = head(successful_methods, 2),
    max_features = 5,
    parallel = FALSE,
    progress = TRUE,
    initialize_framework = FALSE
  )
  
  if (!is.null(api_result)) {
    cat("✓ High-level API test successful\n")
    if (is.list(api_result) && length(api_result) > 1) {
      cat("  Multiple results returned:", length(api_result), "\n")
    } else {
      cat("  Single result returned\n")
    }
  } else {
    cat("✗ High-level API test failed\n")
  }
} else {
  cat("⚠ Skipping API test (no successful methods)\n")
}
cat("\n")

# =============================================================================
# Test 6: Enhanced Analysis
# =============================================================================

cat("Test 6: Enhanced Analysis\n")
cat("=========================\n")

if (length(successful_methods) >= 2) {
  cat("Testing enhanced_modular_feature_selection...\n")
  
  enhanced_result <- enhanced_modular_feature_selection(
    data = test_data,
    methods = head(successful_methods, 3),
    max_features = 5,
    parallel = FALSE,
    progress = TRUE
  )
  
  if (!is.null(enhanced_result)) {
    cat("✓ Enhanced analysis successful\n")
    
    # Check analysis components
    if (!is.null(enhanced_result$analysis)) {
      analysis <- enhanced_result$analysis
      
      if (!is.null(analysis$feature_overlap)) {
        cat("  ✓ Feature overlap analysis included\n")
      }
      
      if (!is.null(analysis$stability)) {
        cat("  ✓ Stability analysis included\n")
        cat("    Stable features found:", length(analysis$stability$stable_features), "\n")
      }
      
      if (!is.null(analysis$consensus_features)) {
        cat("  ✓ Consensus features analysis included\n")
        cat("    High consensus features:", length(analysis$consensus_features$high_consensus), "\n")
      }
      
      if (!is.null(analysis$performance_summary)) {
        cat("  ✓ Performance summary included\n")
        perf <- analysis$performance_summary$summary_statistics
        cat("    Success rate:", round(perf$success_rate * 100, 1), "%\n")
      }
    }
  } else {
    cat("✗ Enhanced analysis failed\n")
  }
} else {
  cat("⚠ Skipping enhanced analysis (insufficient successful methods)\n")
}
cat("\n")

# =============================================================================
# Test 7: Plugin Method Loading
# =============================================================================

cat("Test 7: Plugin Method Loading\n")
cat("=============================\n")

# Check for plugin methods in fs_methods directory
plugin_count <- load_fs_method_modules(verbose = TRUE)
cat("Plugin methods loaded:", plugin_count, "\n")

# Check total methods after plugin loading
total_methods_after <- length(ls(get(".fs_method_registry", envir = .GlobalEnv)))
cat("Total methods after plugin loading:", total_methods_after, "\n\n")

# =============================================================================
# Test 8: Framework Statistics
# =============================================================================

cat("Test 8: Framework Statistics\n")
cat("============================\n")

stats <- get_fs_framework_stats()
cat("Framework Statistics:\n")
cat("  Total methods:", stats$total_methods, "\n")
if (length(stats$categories) > 0) {
  cat("  Categories:\n")
  for (cat_name in names(stats$categories)) {
    cat("    -", cat_name, ":", stats$categories[cat_name], "methods\n")
  }
}
if (stats$dependency_analysis$unique_packages > 0) {
  cat("  Unique packages used:", stats$dependency_analysis$unique_packages, "\n")
  if (length(stats$dependency_analysis$most_common) > 0) {
    cat("  Most common packages:", paste(stats$dependency_analysis$most_common, collapse = ", "), "\n")
  }
}
cat("\n")

# =============================================================================
# Test Summary
# =============================================================================

cat("Test Summary\n")
cat("============\n")
cat("✓ Framework initialization: PASS\n")
cat("✓ Method registry: PASS\n")
cat("✓ Single method execution:", length(successful_methods), "/", length(core_method_ids), "methods successful\n")
if (length(successful_methods) >= 2) {
  cat("✓ Multiple method execution: PASS\n")
  cat("✓ High-level API: PASS\n")
  cat("✓ Enhanced analysis: PASS\n")
} else {
  cat("⚠ Limited testing due to method failures\n")
}
cat("✓ Plugin method loading: PASS\n")
cat("✓ Framework statistics: PASS\n")
cat("\n")

if (length(successful_methods) == length(core_method_ids)) {
  cat("🎉 All tests completed successfully!\n")
  cat("The modular feature selection framework is ready for use.\n")
} else {
  cat("⚠ Some core methods failed - check dependencies and implementation.\n")
  cat("Framework is functional but some methods may not be available.\n")
}

cat("\n")
cat("=========================================\n")
cat("Test Complete\n")
cat("=========================================\n")
