#' Core Feature Selection Methods
#'
#' @description
#' Core feature selection methods that are included with the OmicSelector package.
#' These methods provide essential functionality and serve as examples for the
#' modular framework.
#'
#' @author OmicSelector Team

# =============================================================================
# CORE METHOD IMPLEMENTATIONS
# =============================================================================

#' Student's t-test Feature Selection
#'
#' @description
#' Select features using Student's t-test between classes.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_ttest <- function(data, config = list(), max_features = 20, use_smote = FALSE, timeout_sec = 300) {
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract classes and features
  classes <- data[, 1]
  features <- data[, -1, drop = FALSE]
  
  # Calculate t-test p-values
  p_values <- apply(features, 2, function(x) {
    tryCatch({
      if (length(unique(classes)) == 2) {
        t.test(x ~ classes)$p.value
      } else {
        # For multi-class, use ANOVA
        summary(aov(x ~ classes))[[1]][1, "Pr(>F)"]
      }
    }, error = function(e) 1.0)
  })
  
  # Select top features
  n_select <- min(max_features, length(p_values))
  selected_indices <- order(p_values)[seq_len(n_select)]
  selected_features <- names(features)[selected_indices]
  
  return(list(
    features = selected_features,
    scores = 1 - p_values[selected_indices],  # Convert p-values to scores
    metadata = list(
      method = "Student's t-test",
      p_values = p_values[selected_indices],
      n_features_tested = ncol(features)
    )
  ))
}

#' Correlation Feature Selection
#'
#' @description
#' Select features based on correlation with class variable.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_correlation <- function(data, config = list(), max_features = 20, use_smote = FALSE, timeout_sec = 300) {
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract classes and features
  classes_numeric <- as.numeric(as.factor(data[, 1]))
  features <- data[, -1, drop = FALSE]
  
  # Calculate correlations
  correlations <- apply(features, 2, function(x) {
    tryCatch({
      abs(cor(x, classes_numeric, use = "complete.obs"))
    }, error = function(e) 0)
  })
  
  # Handle NAs
  correlations[is.na(correlations)] <- 0
  
  # Select top features
  n_select <- min(max_features, length(correlations))
  selected_indices <- order(correlations, decreasing = TRUE)[seq_len(n_select)]
  selected_features <- names(features)[selected_indices]
  
  return(list(
    features = selected_features,
    scores = correlations[selected_indices],
    metadata = list(
      method = "Correlation",
      correlations = correlations[selected_indices],
      n_features_tested = ncol(features)
    )
  ))
}

#' Random Forest Importance Feature Selection
#'
#' @description
#' Select features using Random Forest variable importance.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (ntree, mtry, etc.)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_randomforest <- function(data, config = list(), max_features = 20, use_smote = FALSE, timeout_sec = 300) {
  
  # Check if randomForest is available
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    warning("randomForest package not available")
    return(NULL)
  }
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract classes and features
  classes <- as.factor(data[, 1])
  features <- data[, -1, drop = FALSE]
  
  # Default parameters
  ntree <- config$ntree %||% 500
  mtry <- config$mtry %||% sqrt(ncol(features))
  
  # Train Random Forest
  rf_model <- tryCatch({
    randomForest::randomForest(
      x = features,
      y = classes,
      ntree = ntree,
      mtry = mtry,
      importance = TRUE
    )
  }, error = function(e) {
    warning("Random Forest training failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(rf_model)) {
    return(NULL)
  }
  
  # Get importance scores
  importance_scores <- randomForest::importance(rf_model)[, "MeanDecreaseGini"]
  
  # Select top features
  n_select <- min(max_features, length(importance_scores))
  selected_indices <- order(importance_scores, decreasing = TRUE)[seq_len(n_select)]
  selected_features <- names(features)[selected_indices]
  
  return(list(
    features = selected_features,
    scores = importance_scores[selected_indices],
    metadata = list(
      method = "Random Forest",
      importance_type = "MeanDecreaseGini",
      ntree = ntree,
      mtry = mtry,
      n_features_tested = ncol(features)
    )
  ))
}

#' LASSO Feature Selection
#'
#' @description
#' Select features using LASSO regularization.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (alpha, lambda, etc.)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_lasso <- function(data, config = list(), max_features = 20, use_smote = FALSE, timeout_sec = 300) {
  
  # Check if glmnet is available
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    warning("glmnet package not available")
    return(NULL)
  }
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract classes and features
  classes <- as.factor(data[, 1])
  features <- as.matrix(data[, -1, drop = FALSE])
  
  # Default parameters
  alpha <- config$alpha %||% 1.0  # 1.0 for LASSO
  
  # Train LASSO model with cross-validation
  cv_model <- tryCatch({
    glmnet::cv.glmnet(
      x = features,
      y = classes,
      alpha = alpha,
      family = if (length(unique(classes)) == 2) "binomial" else "multinomial"
    )
  }, error = function(e) {
    warning("LASSO training failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(cv_model)) {
    return(NULL)
  }
  
  # Get coefficients at optimal lambda
  coefs <- coef(cv_model, s = "lambda.min")
  
  # Extract non-zero coefficients
  if (is.list(coefs)) {
    # Multinomial case
    all_coefs <- do.call(cbind, coefs)
    feature_importance <- rowSums(abs(all_coefs))[-1]  # Remove intercept
  } else {
    # Binomial case
    feature_importance <- abs(as.numeric(coefs))[-1]  # Remove intercept
  }
  
  # Get non-zero features
  nonzero_features <- which(feature_importance > 0)
  
  if (length(nonzero_features) == 0) {
    warning("LASSO selected no features")
    return(NULL)
  }
  
  # Select top features (up to max_features)
  importance_values <- feature_importance[nonzero_features]
  n_select <- min(max_features, length(nonzero_features))
  top_indices <- order(importance_values, decreasing = TRUE)[seq_len(n_select)]
  selected_indices <- nonzero_features[top_indices]
  selected_features <- colnames(features)[selected_indices]
  
  return(list(
    features = selected_features,
    scores = feature_importance[selected_indices],
    metadata = list(
      method = "LASSO",
      alpha = alpha,
      lambda_min = cv_model$lambda.min,
      n_nonzero_features = length(nonzero_features),
      n_features_tested = ncol(features)
    )
  ))
}

#' Boruta Feature Selection
#'
#' @description
#' Select features using the Boruta algorithm.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (maxRuns, pValue, etc.)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_boruta <- function(data, config = list(), max_features = 20, use_smote = FALSE, timeout_sec = 300) {
  
  # Check if Boruta is available
  if (!requireNamespace("Boruta", quietly = TRUE)) {
    warning("Boruta package not available")
    return(NULL)
  }
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract classes and features
  classes <- as.factor(data[, 1])
  features <- data[, -1, drop = FALSE]
  
  # Default parameters
  maxRuns <- config$maxRuns %||% 100
  pValue <- config$pValue %||% 0.01
  
  # Run Boruta
  boruta_result <- tryCatch({
    Boruta::Boruta(
      x = features,
      y = classes,
      maxRuns = maxRuns,
      pValue = pValue
    )
  }, error = function(e) {
    warning("Boruta failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(boruta_result)) {
    return(NULL)
  }
  
  # Get confirmed features
  confirmed_features <- names(boruta_result$finalDecision[boruta_result$finalDecision == "Confirmed"])
  
  if (length(confirmed_features) == 0) {
    warning("Boruta confirmed no features")
    return(NULL)
  }
  
  # Get importance scores for confirmed features
  importance_scores <- boruta_result$ImpHistory[nrow(boruta_result$ImpHistory), confirmed_features]
  
  # Select top features (up to max_features)
  n_select <- min(max_features, length(confirmed_features))
  selected_indices <- order(importance_scores, decreasing = TRUE)[seq_len(n_select)]
  selected_features <- confirmed_features[selected_indices]
  
  return(list(
    features = selected_features,
    scores = importance_scores[selected_indices],
    metadata = list(
      method = "Boruta",
      maxRuns = maxRuns,
      pValue = pValue,
      n_confirmed = length(confirmed_features),
      n_tentative = sum(boruta_result$finalDecision == "Tentative"),
      n_rejected = sum(boruta_result$finalDecision == "Rejected"),
      n_features_tested = ncol(features)
    )
  ))
}

# =============================================================================
# CORE METHODS REGISTRATION
# =============================================================================

#' Load Core Feature Selection Methods
#'
#' @description
#' Register all core feature selection methods with the framework.
#'
#' @param verbose Show loading messages
#' @return Number of methods loaded
#' @export
load_core_fs_methods <- function(verbose = TRUE) {
  
  if (verbose) {
    cat("Registering core feature selection methods...\n")
  }
  
  methods_loaded <- 0
  
  # Method 1: Student's t-test
  tryCatch({
    register_fs_method(
      method_id = 1,
      method_name = "Student's t-test",
      method_function = fs_method_ttest,
      category = "statistical",
      dependencies = character(0),
      description = "Select features using Student's t-test between classes",
      parameters = list(),
      complexity = "low",
      supports_smote = TRUE,
      timeout_default = 60
    )
    methods_loaded <- methods_loaded + 1
    if (verbose) cat("  ✓ Registered: Student's t-test (ID: 1)\n")
  }, error = function(e) {
    if (verbose) cat("  ✗ Failed to register t-test:", e$message, "\n")
  })
  
  # Method 2: Correlation
  tryCatch({
    register_fs_method(
      method_id = 2,
      method_name = "Correlation",
      method_function = fs_method_correlation,
      category = "statistical",
      dependencies = character(0),
      description = "Select features based on correlation with class variable",
      parameters = list(),
      complexity = "low",
      supports_smote = TRUE,
      timeout_default = 60
    )
    methods_loaded <- methods_loaded + 1
    if (verbose) cat("  ✓ Registered: Correlation (ID: 2)\n")
  }, error = function(e) {
    if (verbose) cat("  ✗ Failed to register correlation:", e$message, "\n")
  })
  
  # Method 3: Random Forest
  tryCatch({
    register_fs_method(
      method_id = 3,
      method_name = "Random Forest",
      method_function = fs_method_randomforest,
      category = "wrapper",
      dependencies = c("randomForest"),
      description = "Select features using Random Forest variable importance",
      parameters = list(ntree = 500, mtry = "sqrt"),
      complexity = "medium",
      supports_smote = TRUE,
      timeout_default = 300
    )
    methods_loaded <- methods_loaded + 1
    if (verbose) cat("  ✓ Registered: Random Forest (ID: 3)\n")
  }, error = function(e) {
    if (verbose) cat("  ✗ Failed to register Random Forest:", e$message, "\n")
  })
  
  # Method 4: LASSO
  tryCatch({
    register_fs_method(
      method_id = 4,
      method_name = "LASSO",
      method_function = fs_method_lasso,
      category = "embedded",
      dependencies = c("glmnet"),
      description = "Select features using LASSO regularization",
      parameters = list(alpha = 1.0),
      complexity = "medium",
      supports_smote = TRUE,
      timeout_default = 180
    )
    methods_loaded <- methods_loaded + 1
    if (verbose) cat("  ✓ Registered: LASSO (ID: 4)\n")
  }, error = function(e) {
    if (verbose) cat("  ✗ Failed to register LASSO:", e$message, "\n")
  })
  
  # Method 5: Boruta
  tryCatch({
    register_fs_method(
      method_id = 5,
      method_name = "Boruta",
      method_function = fs_method_boruta,
      category = "wrapper",
      dependencies = c("Boruta"),
      description = "Select features using the Boruta algorithm",
      parameters = list(maxRuns = 100, pValue = 0.01),
      complexity = "high",
      supports_smote = TRUE,
      timeout_default = 600
    )
    methods_loaded <- methods_loaded + 1
    if (verbose) cat("  ✓ Registered: Boruta (ID: 5)\n")
  }, error = function(e) {
    if (verbose) cat("  ✗ Failed to register Boruta:", e$message, "\n")
  })
  
  if (verbose) {
    cat("Core methods registration complete:", methods_loaded, "methods loaded\n")
  }
  
  return(methods_loaded)
}

# =============================================================================
# UTILITIES
# =============================================================================

#' Null-coalescing operator
#'
#' @param x First value
#' @param y Second value (used if x is NULL)
#' @return x if not NULL, otherwise y
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Apply SMOTE Safely
#'
#' @description
#' Apply SMOTE balancing with error handling.
#'
#' @param data Data frame with class column as first column
#' @return Balanced data frame or original data if SMOTE fails
apply_smote_safely <- function(data) {
  
  if (!requireNamespace("SMOTE", quietly = TRUE) && 
      !requireNamespace("smotefamily", quietly = TRUE) &&
      !requireNamespace("DMwR", quietly = TRUE)) {
    warning("No SMOTE package available. Using original data.")
    return(data)
  }
  
  tryCatch({
    # Try DMwR first (most common)
    if (requireNamespace("DMwR", quietly = TRUE)) {
      result <- DMwR::SMOTE(as.formula(paste(names(data)[1], "~ .")), data)
      return(result)
    }
    
    # Fallback to other packages
    if (requireNamespace("smotefamily", quietly = TRUE)) {
      classes <- data[, 1]
      features <- data[, -1, drop = FALSE]
      result <- smotefamily::SMOTE(features, classes)
      balanced_data <- data.frame(class = result$data$class, result$data[, -ncol(result$data)])
      names(balanced_data)[1] <- names(data)[1]
      return(balanced_data)
    }
    
    # If all else fails, return original data
    warning("SMOTE balancing failed. Using original data.")
    return(data)
    
  }, error = function(e) {
    warning("SMOTE balancing failed: ", e$message, ". Using original data.")
    return(data)
  })
}
