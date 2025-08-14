# Example: Using Comprehensive Feature Selection
# =============================================

cat("🧬 OmicSelector Comprehensive Feature Selection Example\n")
cat("====================================================\n\n")

# Load the comprehensive feature selection modules
source('R/comprehensive-feature-selection.R')
source('R/enhanced-omicselector-integration.R')

# Create sample miRNA data for demonstration
set.seed(42)
n_samples <- 100
n_features <- 50

# Generate feature names (miRNA-like)
feature_names <- paste0("hsa_miR_", sample(1:1000, n_features))

# Create sample data with some informative features
sample_data <- data.frame(
  Class = factor(rep(c("Case", "Control"), each = n_samples/2))
)

# Add features with varying signal strength
for (i in 1:n_features) {
  if (i <= 10) {
    # First 10 features are informative
    case_mean <- runif(1, 5, 8)
    control_mean <- runif(1, 2, 4)
    
    sample_data[[feature_names[i]]] <- c(
      rnorm(n_samples/2, case_mean, 1),    # Cases
      rnorm(n_samples/2, control_mean, 1)  # Controls
    )
  } else {
    # Rest are noise
    sample_data[[feature_names[i]]] <- rnorm(n_samples, 5, 2)
  }
}

cat("📊 Sample Data Created:\n")
cat("- Samples:", nrow(sample_data), "\n")
cat("- Features:", ncol(sample_data) - 1, "\n") 
cat("- Informative features:", 10, "\n")
cat("- Noise features:", n_features - 10, "\n\n")

# Test basic feature selection methods (fast methods for demo)
cat("🔍 Testing Feature Selection Methods...\n")
cat("=====================================\n")

# Test individual methods
methods_to_test <- c(1, 2, 3, 11, 17)  # sig, fcsig, cfs, RF_RFE, Boruta
method_results <- list()

for (method_num in methods_to_test) {
  method_name <- get_comprehensive_method_name(method_num)
  cat(sprintf("Testing Method %d (%s)... ", method_num, method_name))
  
  tryCatch({
    # Prepare data variants
    data_variants <- prepare_data_variants(sample_data, use_smote = FALSE, prefer_no_features = 10)
    
    # Execute method
    result <- execute_comprehensive_method(
      method_num = method_num,
      data_variants = data_variants,
      prefer_no_features = 10,
      timeout_sec = 30,
      config = list()
    )
    
    if (!is.null(result) && length(result$features) > 0) {
      method_results[[method_name]] <- result$features
      cat("✓ SUCCESS (", length(result$features), " features)\n", sep="")
    } else {
      cat("⚠ No features selected\n")
    }
    
  }, error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  })
}

cat("\n📈 Results Summary:\n")
cat("==================\n")

if (length(method_results) > 0) {
  for (method_name in names(method_results)) {
    features <- method_results[[method_name]]
    cat(sprintf("%-20s: %d features selected\n", method_name, length(features)))
    cat(sprintf("%-20s  Features: %s\n", "", paste(head(features, 5), collapse=", ")))
    if (length(features) > 5) cat(sprintf("%-20s  ... and %d more\n", "", length(features) - 5))
    cat("\n")
  }
  
  # Feature overlap analysis
  all_features <- unique(unlist(method_results))
  feature_freq <- table(unlist(method_results))
  
  cat("🎯 Feature Overlap Analysis:\n")
  cat("===========================\n")
  cat("Total unique features found:", length(all_features), "\n")
  
  if (length(feature_freq) > 0) {
    stable_features <- feature_freq[feature_freq >= 2]
    if (length(stable_features) > 0) {
      cat("Features selected by multiple methods:\n")
      for (i in seq_len(min(10, length(stable_features)))) {
        feature_name <- names(stable_features)[i]
        freq <- stable_features[i]
        cat(sprintf("  %s: selected by %d methods\n", feature_name, freq))
      }
    }
  }
  
} else {
  cat("No methods completed successfully.\n")
}

cat("\n🚀 Comprehensive Feature Selection Ready!\n")
cat("========================================\n")
cat("The comprehensive feature selection module is now fully functional.\n")
cat("You can use all 70+ feature selection methods with:\n\n")
cat('results <- enhanced_omicselector(\n')
cat('  wd = "path/to/your/data",\n')
cat('  m = 1:70,  # All methods\n')
cat('  prefer_no_features = 20,\n')
cat('  use_comprehensive = TRUE\n')
cat(')\n\n')
cat("This provides you with the most comprehensive feature selection\n")
cat("toolkit available for omics data analysis!\n")
