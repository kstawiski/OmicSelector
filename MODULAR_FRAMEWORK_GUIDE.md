# Modular Feature Selection Framework - Developer Guide

## Overview

The OmicSelector package now includes a highly modular, plugin-based feature selection framework that makes adding new methods extremely easy. This guide shows you how to add new feature selection methods in just a few simple steps.

## Quick Start: Adding a New Method

### Step 1: Copy the Template

1. Navigate to `R/fs_methods/`
2. Copy `fs_method_template.R` to `fs_method_[your_method_name].R`
3. Replace the placeholders with your method implementation

### Step 2: Implement Your Method

```r
# Example: Adding a Simple Variance Filter
fs_method_variance_filter <- function(data, config = list(), max_features = 20, use_smote = FALSE, timeout_sec = 300) {
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract features (skip class column)
  features <- data[, -1, drop = FALSE]
  
  # Get configuration
  min_variance <- config$min_variance %||% 0.01
  
  # Calculate variance for each feature
  variances <- apply(features, 2, var, na.rm = TRUE)
  
  # Filter by minimum variance and select top features
  valid_features <- variances >= min_variance
  filtered_variances <- variances[valid_features]
  
  if (length(filtered_variances) == 0) {
    warning("No features meet variance threshold")
    return(NULL)
  }
  
  # Select top features by variance
  n_select <- min(max_features, length(filtered_variances))
  selected_indices <- order(filtered_variances, decreasing = TRUE)[seq_len(n_select)]
  selected_features <- names(filtered_variances)[selected_indices]
  
  return(list(
    features = selected_features,
    scores = filtered_variances[selected_indices],
    metadata = list(
      method = "Variance Filter",
      n_features_input = ncol(features),
      n_features_selected = length(selected_features),
      parameters_used = list(
        min_variance = min_variance,
        max_features = max_features
      )
    )
  ))
}

# Register the method
register_fs_method(
  method_id = 201,  # Use unique ID (200+ for custom methods)
  method_name = "Variance Filter",
  method_function = fs_method_variance_filter,
  category = "filter",
  dependencies = character(0),
  description = "Select features with highest variance above threshold",
  parameters = list(min_variance = 0.01),
  complexity = "low",
  supports_smote = TRUE,
  timeout_default = 60
)
```

### Step 3: Test Your Method

```r
# Load the framework
library(OmicSelector)
initialize_fs_framework()

# Test your method
result <- execute_fs_method(
  method_id = 201,
  data = your_data,
  max_features = 10
)

print(result$features)
```

That's it! Your method is now available throughout the framework.

## Framework Architecture

### Core Components

1. **Method Registry**: Central repository for all feature selection methods
2. **Execution Engine**: Handles method execution with error handling and timeouts
3. **Dependency Management**: Smart loading of optional packages
4. **Analysis Framework**: Provides overlap, stability, and consensus analysis
5. **Plugin System**: Auto-discovers and loads method modules

### Method Categories

- **statistical**: Statistical tests (t-test, ANOVA, chi-square)
- **filter**: Filter methods (correlation, variance, mutual information)
- **wrapper**: Wrapper methods (RFE, forward/backward selection)
- **embedded**: Embedded methods (LASSO, Ridge, Elastic Net)
- **ensemble**: Ensemble methods (Boruta, stability selection)
- **hybrid**: Hybrid approaches combining multiple strategies

### Complexity Levels

- **low**: Fast methods (< 1 minute on typical datasets)
- **medium**: Moderate methods (1-10 minutes on typical datasets)  
- **high**: Slow methods (> 10 minutes on typical datasets)

## Method Implementation Guidelines

### Required Function Signature

```r
fs_method_[name] <- function(data, config = list(), max_features = 20, use_smote = FALSE, timeout_sec = 300)
```

### Required Parameters

- `data`: Data frame with class column as first column
- `config`: Named list of method-specific parameters
- `max_features`: Maximum number of features to select
- `use_smote`: Whether to apply SMOTE balancing
- `timeout_sec`: Timeout in seconds

### Required Return Value

```r
list(
  features = character(),      # Selected feature names
  scores = numeric(),          # Feature scores (optional)
  metadata = list()            # Method execution information
)
```

### Best Practices

1. **Error Handling**: Always use `tryCatch()` for robust error handling
2. **Input Validation**: Check data format and handle edge cases
3. **Dependency Checking**: Check required packages with `requireNamespace()`
4. **Parameter Defaults**: Use `%||%` operator for default values
5. **SMOTE Support**: Use `apply_smote_safely()` if method supports SMOTE
6. **Informative Metadata**: Include comprehensive execution information

### Example Template Structure

```r
#' Method Name Feature Selection
#'
#' @description Brief description of what the method does
#' @param data Training data with class column as first column
#' @param config Configuration list
#' @param max_features Maximum features to select
#' @param use_smote Apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#' @return List with features, scores, and metadata
fs_method_[name] <- function(data, config = list(), max_features = 20, use_smote = FALSE, timeout_sec = 300) {
  
  # 1. Check dependencies
  if (!requireNamespace("required_package", quietly = TRUE)) {
    warning("required_package not available")
    return(NULL)
  }
  
  # 2. Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # 3. Extract data components
  classes <- data[, 1]
  features <- data[, -1, drop = FALSE]
  
  # 4. Get configuration parameters
  param1 <- config$param1 %||% default_value1
  param2 <- config$param2 %||% default_value2
  
  # 5. Implement your algorithm
  tryCatch({
    # Your feature selection logic here
    feature_scores <- calculate_scores(features, classes, param1, param2)
    
    # Handle edge cases
    feature_scores[is.na(feature_scores)] <- 0
    
    # Select top features
    n_select <- min(max_features, length(feature_scores))
    selected_indices <- order(feature_scores, decreasing = TRUE)[seq_len(n_select)]
    selected_features <- names(features)[selected_indices]
    
    # 6. Return results
    return(list(
      features = selected_features,
      scores = feature_scores[selected_indices],
      metadata = list(
        method = "Method Name",
        n_features_input = ncol(features),
        n_features_selected = length(selected_features),
        parameters_used = list(
          param1 = param1,
          param2 = param2
        )
      )
    ))
    
  }, error = function(e) {
    warning("Method execution failed: ", e$message)
    return(NULL)
  })
}

# Register the method
register_fs_method(
  method_id = unique_id,
  method_name = "Method Name",
  method_function = fs_method_[name],
  category = "appropriate_category",
  dependencies = c("required_packages"),
  description = "Brief description",
  parameters = list(param1 = default1, param2 = default2),
  complexity = "low|medium|high",
  supports_smote = TRUE|FALSE,
  timeout_default = seconds
)
```

## Advanced Features

### Custom Configuration

Methods can accept complex configurations:

```r
config <- list(
  algorithm_params = list(
    learning_rate = 0.01,
    max_iterations = 1000
  ),
  preprocessing = list(
    normalize = TRUE,
    remove_correlated = 0.95
  ),
  selection_criteria = list(
    min_score = 0.1,
    max_features = 50
  )
)

result <- execute_fs_method(method_id, data, config = config)
```

### Parallel Execution

The framework automatically handles parallel execution:

```r
# Execute multiple methods in parallel
results <- execute_multiple_fs_methods(
  method_ids = c(1, 2, 3, 4, 5),
  data = data,
  parallel = TRUE,
  timeout_sec = 300
)
```

### Method Dependencies

Declare optional dependencies for graceful degradation:

```r
register_fs_method(
  method_id = 301,
  dependencies = c("randomForest", "caret", "doParallel"),  # Optional packages
  # ... other parameters
)
```

### Performance Monitoring

Methods automatically track execution statistics:

```r
result <- execute_fs_method(method_id, data)
print(result$metadata$execution_time)  # Execution time in seconds
print(result$metadata$execution_start) # Start timestamp
print(result$metadata$config_used)     # Configuration used
```

## Integration with Existing Code

### Backward Compatibility

The modular framework works alongside existing OmicSelector functions:

```r
# Use with existing comprehensive system
comprehensive_results <- run_comprehensive_feature_selection(data)

# Use with new modular system
modular_results <- modular_feature_selection(data, methods = c(1, 2, 3))

# Combine results
all_features <- unique(c(
  unlist(lapply(comprehensive_results, function(x) x$features)),
  unlist(lapply(modular_results, function(x) x$features))
))
```

### Enhanced Analysis

Get comprehensive analysis across all methods:

```r
enhanced_results <- enhanced_modular_feature_selection(
  data = data,
  methods = c(1, 2, 3, 4, 5),
  analysis_options = list(
    feature_overlap = TRUE,
    stability_analysis = TRUE,
    consensus_features = TRUE,
    performance_summary = TRUE
  )
)

# Access consensus features
consensus <- enhanced_results$analysis$consensus_features
print(consensus$high_consensus)  # Features selected by >75% of methods
```

## Method ID Guidelines

- **1-100**: Reserved for core methods
- **101-200**: Package-provided plugin methods
- **201-1000**: User-defined custom methods
- **1001+**: Reserved for future extensions

## Testing Your Methods

### Basic Testing

```r
# Test with sample data
test_data <- data.frame(
  class = factor(rep(c("A", "B"), each = 50)),
  feature1 = rnorm(100, mean = rep(c(0, 1), each = 50)),
  feature2 = rnorm(100),
  feature3 = rnorm(100, mean = rep(c(1, 0), each = 50))
)

result <- execute_fs_method(your_method_id, test_data, max_features = 2)
```

### Comprehensive Testing

```r
# Test with multiple datasets and configurations
test_configs <- list(
  list(param1 = 0.1),
  list(param1 = 0.5),
  list(param1 = 1.0)
)

for (config in test_configs) {
  result <- execute_fs_method(your_method_id, test_data, config = config)
  print(paste("Config:", deparse(config), "Selected:", length(result$features)))
}
```

## Contributing Methods

To contribute methods to the official OmicSelector package:

1. Implement your method following these guidelines
2. Add comprehensive tests and documentation
3. Ensure compatibility across different data types
4. Submit a pull request with your method file

## Troubleshooting

### Common Issues

1. **Method not found**: Ensure `initialize_fs_framework()` was called
2. **Dependency errors**: Check that required packages are installed
3. **Empty results**: Verify data format (class column first)
4. **Timeout errors**: Increase `timeout_sec` for complex methods

### Debug Mode

Enable detailed logging:

```r
# Check method registration
methods_df <- list_fs_methods()
print(methods_df[methods_df$id == your_method_id, ])

# Check dependencies
dep_check <- check_method_dependencies(your_method_id)
print(dep_check)

# Test method execution with error details
result <- tryCatch({
  execute_fs_method(your_method_id, data)
}, error = function(e) {
  print(paste("Error:", e$message))
  return(NULL)
})
```

This modular framework makes OmicSelector infinitely extensible while maintaining simplicity and robustness. Happy method development!
