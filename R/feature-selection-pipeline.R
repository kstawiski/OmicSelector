#' Feature Selection Pipeline
#'
#' @description
#' Modular feature selection pipeline that runs multiple methods
#'
#' @param results OmicSelector results object
#' @param methods Integer vector with method numbers
#' @param config Configuration list
#' @return Updated results object
#' @keywords internal
run_feature_selection_pipeline <- function(results, methods, config) {
  
  # Create progress bar
  pb <- create_progress_bar(length(methods), "Running feature selection methods")
  
  # Set up method execution environment
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  
  for (i in seq_along(methods)) {
    method_num <- methods[i]
    method_name <- get_method_name(method_num)
    
    log_info("Running method {method_num}: {method_name}")
    
    tryCatch({
      # Run individual method with timeout
      method_result <- run_single_method(
        method_num, 
        config,
        timeout = config$timeout_sec
      )
      
      if (!is.null(method_result$formula)) {
        results$formulas[[length(results$formulas) + 1]] <- method_result$formula
        results$method_results[[as.character(method_num)]] <- method_result
        log_info("Method {method_num} completed successfully")
      } else {
        log_warn("Method {method_num} returned no formula")
      }
      
    }, error = function(e) {
      error_msg <- paste("Method", method_num, "failed:", e$message)
      results$errors[[as.character(method_num)]] <- error_msg
      log_error(error_msg)
    })
    
    update_progress_bar(pb, i)
  }
  
  finish_progress_bar(pb)
  
  log_info("Feature selection pipeline completed")
  log_info("Generated {length(results$formulas)} formulas from {length(methods)} methods")
  
  return(results)
}

#' Get Method Name
#'
#' @param method_num Integer method number
#' @return Character string with method name
#' @keywords internal
get_method_name <- function(method_num) {
  
  method_names <- c(
    "1" = "Sig (t-test p<0.05)",
    "2" = "Fcsig (t-test + FC filter)",
    "3" = "CFS (Correlation-based Feature Selection)",
    "4" = "Classloop",
    "5" = "FCFS (Forward CFS)",
    "6" = "MDL_AUC",
    "7" = "MDL_SU", 
    "8" = "MDL_CorrSF",
    "9" = "bounceR",
    "10" = "RandomForestRFE_CV",
    "11" = "RandomForestRFE",
    "12" = "GeneticAlgorithmRF",
    "13" = "SimulatedAnnealing",
    "14" = "Boruta_Confirmed",
    "15" = "Boruta_Tentative",
    "16" = "Boruta_Confirmed_Tentative",
    "17" = "Boruta",
    "18" = "spFSR",
    "19" = "varSelRF",
    "20" = "varSelRF",
    "21" = "WxNet"
  )
  
  # Add more method names as needed
  if (method_num <= 70 && !as.character(method_num) %in% names(method_names)) {
    return(paste("Method", method_num))
  }
  
  return(method_names[as.character(method_num)] %||% paste("Method", method_num))
}

#' Run Single Feature Selection Method
#'
#' @param method_num Integer method number
#' @param config Configuration list
#' @param timeout Timeout in seconds
#' @return Method result list
#' @keywords internal
run_single_method <- function(method_num, config, timeout = 3600) {
  
  # Implement timeout mechanism
  result <- NULL
  
  if (requireNamespace("future", quietly = TRUE)) {
    # Use future for timeout
    f <- future::future({
      execute_method(method_num, config)
    })
    
    result <- tryCatch({
      future::value(f, timeout = timeout)
    }, error = function(e) {
      if (grepl("timeout", e$message, ignore.case = TRUE)) {
        stop("Method ", method_num, " timed out after ", timeout, " seconds")
      } else {
        stop(e$message)
      }
    })
    
  } else {
    # Fallback without timeout
    result <- execute_method(method_num, config)
  }
  
  return(result)
}

#' Execute Feature Selection Method
#'
#' @param method_num Integer method number
#' @param config Configuration list
#' @return Method result list
#' @keywords internal
execute_method <- function(method_num, config) {
  
  # Load required data
  train_data <- safe_read_csv("mixed_train.csv")
  if (is.null(train_data)) {
    stop("Failed to load training data")
  }
  
  # Execute method based on number
  result <- switch(as.character(method_num),
    "1" = method_sig(train_data, config),
    "2" = method_fcsig(train_data, config),
    "3" = method_cfs(train_data, config),
    "11" = method_rf_rfe(train_data, config),
    "17" = method_boruta(train_data, config),
    "20" = method_varselrf(train_data, config),
    {
      # For methods not yet implemented, return a placeholder
      log_warn("Method {method_num} not implemented in refactored version")
      return(list(formula = NULL, features = NULL))
    }
  )
  
  return(result)
}

#' Method 1: Significant Features (t-test)
#'
#' @param data Training data
#' @param config Configuration
#' @return Method result
#' @keywords internal
method_sig <- function(data, config) {
  
  log_debug("Running Sig method (t-test p<0.05)")
  
  # Assume first column is class, rest are features
  class_col <- data[, 1]
  feature_cols <- data[, -1]
  
  # Perform t-tests
  p_values <- numeric(ncol(feature_cols))
  
  for (i in seq_along(p_values)) {
    tryCatch({
      test_result <- t.test(feature_cols[, i] ~ class_col)
      p_values[i] <- test_result$p.value
    }, error = function(e) {
      p_values[i] <- 1.0
    })
  }
  
  # Apply multiple testing correction
  p_adjusted <- p.adjust(p_values, method = "BH")
  
  # Select significant features
  significant <- which(p_adjusted < 0.05)
  
  if (length(significant) == 0) {
    return(list(formula = NULL, features = NULL))
  }
  
  # Limit to preferred number of features
  if (length(significant) > config$prefer_no_features) {
    # Order by p-value and take top features
    top_indices <- order(p_adjusted[significant])[1:config$prefer_no_features]
    significant <- significant[top_indices]
  }
  
  selected_features <- colnames(feature_cols)[significant]
  formula <- paste("class ~", paste(selected_features, collapse = " + "))
  
  log_debug("Sig method selected {length(selected_features)} features")
  
  return(list(
    formula = formula,
    features = selected_features,
    p_values = p_adjusted[significant]
  ))
}

#' Method 2: Fold Change + Significance
#'
#' @param data Training data
#' @param config Configuration
#' @return Method result
#' @keywords internal
method_fcsig <- function(data, config) {
  
  log_debug("Running Fcsig method (t-test + FC filter)")
  
  # Get significant features first
  sig_result <- method_sig(data, config)
  
  if (is.null(sig_result$features)) {
    return(list(formula = NULL, features = NULL))
  }
  
  # Calculate fold changes
  class_col <- data[, 1]
  feature_cols <- data[, -1]
  
  unique_classes <- unique(class_col)
  if (length(unique_classes) != 2) {
    log_warn("Fold change calculation requires exactly 2 classes")
    return(sig_result)
  }
  
  class1_data <- feature_cols[class_col == unique_classes[1], ]
  class2_data <- feature_cols[class_col == unique_classes[2], ]
  
  mean1 <- colMeans(class1_data, na.rm = TRUE)
  mean2 <- colMeans(class2_data, na.rm = TRUE)
  
  # Calculate log2 fold change
  log2fc <- log2((mean2 + 1e-6) / (mean1 + 1e-6))
  
  # Filter by fold change threshold
  fc_threshold <- 1.0  # |log2FC| > 1
  fc_filter <- abs(log2fc) > fc_threshold
  
  # Get features that pass both significance and FC filter
  sig_indices <- match(sig_result$features, colnames(feature_cols))
  final_indices <- sig_indices[fc_filter[sig_indices]]
  final_indices <- final_indices[!is.na(final_indices)]
  
  if (length(final_indices) == 0) {
    return(list(formula = NULL, features = NULL))
  }
  
  selected_features <- colnames(feature_cols)[final_indices]
  formula <- paste("class ~", paste(selected_features, collapse = " + "))
  
  log_debug("Fcsig method selected {length(selected_features)} features")
  
  return(list(
    formula = formula,
    features = selected_features,
    log2fc = log2fc[final_indices]
  ))
}

#' Method 3: Correlation-based Feature Selection
#'
#' @param data Training data  
#' @param config Configuration
#' @return Method result
#' @keywords internal
method_cfs <- function(data, config) {
  
  log_debug("Running CFS method")
  
  # Check if required package is available
  if (!requireNamespace("FSelector", quietly = TRUE)) {
    log_warn("FSelector package not available for CFS method")
    return(list(formula = NULL, features = NULL))
  }
  
  # This is a simplified version - full implementation would use FSelector
  # For now, return top correlated features with class
  class_col <- as.numeric(as.factor(data[, 1]))
  feature_cols <- data[, -1]
  
  correlations <- numeric(ncol(feature_cols))
  for (i in seq_along(correlations)) {
    correlations[i] <- abs(cor(feature_cols[, i], class_col, use = "complete.obs"))
  }
  
  # Remove NAs
  correlations[is.na(correlations)] <- 0
  
  # Select top correlated features
  n_select <- min(config$prefer_no_features, sum(correlations > 0))
  if (n_select == 0) {
    return(list(formula = NULL, features = NULL))
  }
  
  top_indices <- order(correlations, decreasing = TRUE)[1:n_select]
  selected_features <- colnames(feature_cols)[top_indices]
  formula <- paste("class ~", paste(selected_features, collapse = " + "))
  
  log_debug("CFS method selected {length(selected_features)} features")
  
  return(list(
    formula = formula,
    features = selected_features,
    correlations = correlations[top_indices]
  ))
}

#' Method 11: Random Forest RFE
#'
#' @param data Training data
#' @param config Configuration
#' @return Method result
#' @keywords internal
method_rf_rfe <- function(data, config) {
  
  log_debug("Running Random Forest RFE method")
  
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    log_warn("randomForest package not available")
    return(list(formula = NULL, features = NULL))
  }
  
  class_col <- as.factor(data[, 1])
  feature_cols <- data[, -1]
  
  # Build initial random forest
  rf_model <- randomForest::randomForest(
    x = feature_cols,
    y = class_col,
    importance = TRUE,
    ntree = 100
  )
  
  # Get feature importance
  importance_scores <- randomForest::importance(rf_model)[, "MeanDecreaseGini"]
  
  # Select top features
  n_select <- min(config$prefer_no_features, length(importance_scores))
  top_indices <- order(importance_scores, decreasing = TRUE)[1:n_select]
  
  selected_features <- colnames(feature_cols)[top_indices]
  formula <- paste("class ~", paste(selected_features, collapse = " + "))
  
  log_debug("RF RFE method selected {length(selected_features)} features")
  
  return(list(
    formula = formula,
    features = selected_features,
    importance = importance_scores[top_indices]
  ))
}

#' Method 17: Boruta
#'
#' @param data Training data
#' @param config Configuration
#' @return Method result
#' @keywords internal
method_boruta <- function(data, config) {
  
  log_debug("Running Boruta method")
  
  if (!requireNamespace("Boruta", quietly = TRUE)) {
    log_warn("Boruta package not available")
    return(list(formula = NULL, features = NULL))
  }
  
  class_col <- as.factor(data[, 1])
  feature_cols <- data[, -1]
  
  # Run Boruta
  boruta_result <- Boruta::Boruta(
    x = feature_cols,
    y = class_col,
    maxRuns = config$max_iterations * 10
  )
  
  # Get confirmed features
  confirmed_features <- names(boruta_result$finalDecision[boruta_result$finalDecision == "Confirmed"])
  
  if (length(confirmed_features) == 0) {
    return(list(formula = NULL, features = NULL))
  }
  
  # Limit to preferred number
  if (length(confirmed_features) > config$prefer_no_features) {
    # Order by importance and take top
    importance_stats <- boruta_result$ImpHistory
    final_importance <- apply(importance_stats, 2, median, na.rm = TRUE)
    confirmed_importance <- final_importance[confirmed_features]
    top_indices <- order(confirmed_importance, decreasing = TRUE)[1:config$prefer_no_features]
    confirmed_features <- confirmed_features[top_indices]
  }
  
  formula <- paste("class ~", paste(confirmed_features, collapse = " + "))
  
  log_debug("Boruta method selected {length(confirmed_features)} features")
  
  return(list(
    formula = formula,
    features = confirmed_features,
    decision = boruta_result$finalDecision[confirmed_features]
  ))
}

#' Method 20: Variable Selection Random Forest
#'
#' @param data Training data
#' @param config Configuration
#' @return Method result
#' @keywords internal
method_varselrf <- function(data, config) {
  
  log_debug("Running varSelRF method")
  
  if (!requireNamespace("varSelRF", quietly = TRUE)) {
    log_warn("varSelRF package not available")
    return(list(formula = NULL, features = NULL))
  }
  
  class_col <- as.factor(data[, 1])
  feature_cols <- data[, -1]
  
  # Run variable selection RF
  varsel_result <- varSelRF::varSelRF(
    xdata = feature_cols,
    Class = class_col,
    ntree = 100,
    ntreeIterat = 50,
    vars.drop.frac = 0.2
  )
  
  selected_features <- varsel_result$selected.vars
  
  if (length(selected_features) == 0) {
    return(list(formula = NULL, features = NULL))
  }
  
  formula <- paste("class ~", paste(selected_features, collapse = " + "))
  
  log_debug("varSelRF method selected {length(selected_features)} features")
  
  return(list(
    formula = formula,
    features = selected_features,
    oob_error = varsel_result$selec.history$OOB
  ))
}
