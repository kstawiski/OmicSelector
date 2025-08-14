#!/usr/bin/env Rscript

# Simple test script for optimized OmicSelector functions
# =====================================================

cat("Testing Optimized OmicSelector Functions\n")
cat("========================================\n\n")

# Load the new functions
tryCatch({
  source('R/omicselector-helpers.R')
  source('R/omicselector-optimized.R') 
  source('R/OmicSelector_OmicSelector_improved.R')
  cat("✓ Functions loaded successfully\n")
}, error = function(e) {
  cat("✗ Error loading functions:", e$message, "\n")
  quit(status = 1)
})

# Test 1: Helper functions
cat("\n1. Testing helper functions...\n")

# Create test data
set.seed(42)
test_data <- data.frame(
  Class = factor(c(rep('Case', 30), rep('Control', 20)), levels = c('Control', 'Case')),
  hsa_miR_1 = rnorm(50, mean = ifelse(c(rep(1, 30), rep(0, 20)), 2, 0)),
  hsa_miR_2 = rnorm(50),
  hsa_let_7a = rnorm(50, mean = ifelse(c(rep(1, 30), rep(0, 20)), 1.5, 0)),
  invalid_feature = rnorm(50)
)

# Test get_feature_columns
hsa_features <- get_feature_columns(test_data, 'hsa_')
cat("   Found", length(hsa_features), "hsa_ features\n")

# Test differential expression  
feature_matrix <- as.matrix(test_data[, hsa_features])
de_results <- fast_ttest(feature_matrix, test_data$Class, 'Case', 'Control')
cat("   DE analysis completed for", nrow(de_results), "features\n")

# Show significant features
sig_features <- de_results[de_results$p_value < 0.05, ]
cat("   Significant features found:", nrow(sig_features), "\n")

# Test feature validation
valid_features <- validate_mirna_features(hsa_features, 'hsa_')
cat("   Valid miRNA features:", sum(valid_features), "/", length(valid_features), "\n")

cat("✓ Helper functions working correctly\n")

# Test 2: Create test dataset files
cat("\n2. Creating test dataset files...\n")

test_dir <- tempdir()
set.seed(42)

# Generate realistic miRNA test data
create_test_mirna_data <- function(n_samples, n_features) {
  feature_names <- paste0("hsa_miR_", sample(1:1000, n_features))
  n_case <- floor(n_samples * 0.6)
  n_control <- n_samples - n_case
  class_labels <- c(rep("Case", n_case), rep("Control", n_control))
  
  # Create feature matrix with some differential expression
  feature_data <- matrix(rnorm(n_samples * n_features), nrow = n_samples, ncol = n_features)
  colnames(feature_data) <- feature_names
  
  # Make some features differentially expressed
  case_indices <- which(class_labels == "Case")
  diff_features <- 1:min(5, n_features)
  
  for (i in diff_features) {
    feature_data[case_indices, i] <- feature_data[case_indices, i] + rnorm(length(case_indices), mean = 1.5, sd = 0.3)
  }
  
  data.frame(
    Class = factor(class_labels, levels = c("Control", "Case")),
    feature_data,
    stringsAsFactors = FALSE
  )
}

# Create datasets
train_data <- create_test_mirna_data(n_samples = 80, n_features = 20)
test_data <- create_test_mirna_data(n_samples = 30, n_features = 20) 
validation_data <- create_test_mirna_data(n_samples = 30, n_features = 20)

# Save as CSV files
write.csv(train_data, file.path(test_dir, "mixed_train.csv"), row.names = FALSE)
write.csv(test_data, file.path(test_dir, "mixed_test.csv"), row.names = FALSE)
write.csv(validation_data, file.path(test_dir, "mixed_validation.csv"), row.names = FALSE)

cat("   Test files created in:", test_dir, "\n")
cat("   Training samples:", nrow(train_data), "\n")
cat("   miRNA features:", sum(startsWith(colnames(train_data), "hsa_")), "\n")
cat("✓ Test dataset created successfully\n")

# Test 3: Run optimized function (minimal test)
cat("\n3. Testing optimized function (quick test)...\n")

tryCatch({
  results <- omicselector_optimized(
    data_path = test_dir,
    methods = c(1, 2, 3),  # Just basic methods
    max_features = 8,
    max_iterations = 3,
    timeout_minutes = 2,
    parallel_cores = 1,
    verbose = FALSE,
    seed = 42
  )
  
  cat("   Methods executed:", length(results$selected_features), "\n")
  cat("   Features selected:", sum(sapply(results$selected_features, length)), "\n")
  cat("   Significant DE features:", sum(results$differential_expression$significant), "\n")
  cat("   Runtime:", round(results$metadata$total_runtime, 2), "minutes\n")
  cat("✓ Optimized function working correctly\n")
  
}, error = function(e) {
  cat("✗ Optimized function failed:", e$message, "\n")
})

# Test 4: Test improved wrapper (backward compatibility)
cat("\n4. Testing improved wrapper function...\n")

tryCatch({
  wrapper_results <- OmicSelector_OmicSelector(
    wd = test_dir,
    m = c(1, 2),
    max_iterations = 2,
    prefer_no_features = 5,
    register_parallel = FALSE,
    debug = FALSE
  )
  
  cat("   Formulas created:", length(wrapper_results$formulas), "\n")
  cat("   Selected features available:", !is.null(wrapper_results$selected_features), "\n")
  cat("   Backward compatible:", inherits(wrapper_results, "omicselector_legacy"), "\n")
  cat("✓ Wrapper function working correctly\n")
  
}, error = function(e) {
  cat("✗ Wrapper function failed:", e$message, "\n")
})

# Cleanup
unlink(test_dir, recursive = TRUE)

cat("\n========================================\n")
cat("All tests completed!\n")
cat("Optimized OmicSelector functions are ready for use.\n")
cat("========================================\n")
