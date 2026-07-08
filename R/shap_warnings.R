#' @title Correlation-Aware SHAP Interpretation
#'
#' @description
#' Functions to compute SHAP values with warnings about correlated features.
#' High correlations between features can make SHAP values misleading because:
#' 1. Importance can be arbitrarily split between correlated features
#' 2. SHAP values may not reflect true causal relationships
#' 3. Removing one correlated feature could dramatically change others' SHAP values
#'
#' @references
#' Hooker et al. (2021). Unrestricted Permutation Forces Extrapolation.
#' Lundberg & Lee (2017). SHAP: A Unified Approach to Interpreting Model Predictions.
#'
#' @name shap_warnings
NULL


#' @title Compute SHAP Values with Correlation Warnings
#'
#' @description
#' Computes SHAP values for model predictions and flags features that have
#' high correlations with other features, which can make SHAP interpretations
#' unreliable.
#'
#' @param model A trained model (mlr3 Learner, caret model, or predict function)
#' @param data Data.frame of features used for SHAP computation
#' @param n_samples Number of samples for SHAP approximation (default: 100)
#' @param correlation_threshold Correlation threshold for warnings (default: 0.7)
#' @param method SHAP method: "permutation" or "kernel" (default: "permutation")
#' @param verbose Print progress messages
#'
#' @return A list with:
#'   - shap_values: Matrix of SHAP values (samples x features)
#'   - feature_importance: Mean absolute SHAP per feature
#'   - correlations: Correlation matrix of features
#'   - warnings: Data.frame of correlated feature pairs with warnings
#'   - reliable_features: Features with no high correlations (safe to interpret)
#'
#' @examples
#' \dontrun{
#' # Compute SHAP with correlation warnings
#' result <- compute_shap_with_warnings(model, test_data)
#'
#' # Check which features have reliable SHAP interpretations
#' print(result$reliable_features)
#'
#' # View correlation warnings
#' print(result$warnings)
#' }
#'
#' @export
compute_shap_with_warnings <- function(model,
                                        data,
                                        n_samples = 100,
                                        correlation_threshold = 0.7,
                                        method = "permutation",
                                        verbose = TRUE) {

  # Input validation
  stopifnot(
    "data must be a data.frame" = is.data.frame(data),
    "n_samples must be positive" = n_samples > 0,
    "correlation_threshold must be between 0 and 1" =
      correlation_threshold > 0 && correlation_threshold < 1,
    "method must be 'permutation' or 'kernel'" =
      method %in% c("permutation", "kernel")
  )

  if (!requireNamespace("DALEX", quietly = TRUE)) {
    stop("Package 'DALEX' is required for SHAP computation. Install with: install.packages('DALEX')",
         call. = FALSE)
  }

  # Step 1: Compute feature correlations
  if (verbose) message("Computing feature correlations...")

  numeric_cols <- sapply(data, is.numeric)
  numeric_data <- data[, numeric_cols, drop = FALSE]

  if (ncol(numeric_data) < 2) {
    stop("Need at least 2 numeric features for correlation analysis", call. = FALSE)
  }

  cor_matrix <- cor(numeric_data, use = "pairwise.complete.obs")

  # Step 2: Identify highly correlated pairs
  high_cor_pairs <- find_high_correlations(cor_matrix, correlation_threshold)

  # Step 3: Create DALEX explainer
  if (verbose) message("Creating explainer...")

  # Handle different model types
  predict_fn <- create_predict_function(model)

  explainer <- DALEX::explain(
    model = model,
    data = numeric_data,
    predict_function = predict_fn,
    verbose = FALSE
  )

  # Step 4: Compute SHAP values
  if (verbose) message(sprintf("Computing SHAP values (%d samples)...", n_samples))

  # Sample observations for SHAP
  sample_idx <- sample(nrow(data), min(n_samples, nrow(data)))
  sample_data <- numeric_data[sample_idx, , drop = FALSE]

  shap_results <- list()

  for (i in seq_len(nrow(sample_data))) {
    shap_i <- DALEX::predict_parts(
      explainer,
      new_observation = sample_data[i, , drop = FALSE],
      type = "shap",
      B = 25  # Number of permutations per observation
    )

    # Extract SHAP values
    shap_df <- as.data.frame(shap_i)
    shap_results[[i]] <- shap_df
  }

  # Aggregate SHAP values
  shap_matrix <- aggregate_shap_values(shap_results, colnames(numeric_data))

  # Step 5: Compute feature importance (mean |SHAP|)
  feature_importance <- colMeans(abs(shap_matrix), na.rm = TRUE)
  feature_importance <- sort(feature_importance, decreasing = TRUE)

  # Step 6: Generate correlation warnings
  warnings_df <- generate_correlation_warnings(
    high_cor_pairs,
    feature_importance,
    correlation_threshold
  )

  # Step 7: Identify reliable features
  correlated_features <- unique(c(high_cor_pairs$feature1, high_cor_pairs$feature2))
  reliable_features <- setdiff(names(feature_importance), correlated_features)

  if (verbose) {
    n_warnings <- nrow(warnings_df)
    n_reliable <- length(reliable_features)

    message(sprintf("\nResults:"))
    message(sprintf("  Total features: %d", ncol(numeric_data)))
    message(sprintf("  Features with high correlations: %d", length(correlated_features)))
    message(sprintf("  Reliable features (safe to interpret): %d", n_reliable))

    if (n_warnings > 0) {
      message(sprintf("\n  WARNING: %d feature pairs have correlation > %.2f",
                      n_warnings, correlation_threshold))
      message("  SHAP values for these features may be unreliable!")
    }
  }

  list(
    shap_values = shap_matrix,
    feature_importance = feature_importance,
    correlations = cor_matrix,
    warnings = warnings_df,
    reliable_features = reliable_features,
    correlated_features = correlated_features,
    threshold = correlation_threshold
  )
}


#' @title Find High Correlation Pairs
#'
#' @description
#' Identifies pairs of features with correlation above a threshold.
#'
#' @param cor_matrix Correlation matrix
#' @param threshold Correlation threshold
#'
#' @return Data.frame with feature1, feature2, correlation columns
#'
#' @keywords internal
find_high_correlations <- function(cor_matrix, threshold) {
  # Get upper triangle pairs
  pairs <- which(abs(cor_matrix) > threshold & upper.tri(cor_matrix), arr.ind = TRUE)

  if (nrow(pairs) == 0) {
    return(data.frame(
      feature1 = character(0),
      feature2 = character(0),
      correlation = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    feature1 = rownames(cor_matrix)[pairs[, 1]],
    feature2 = colnames(cor_matrix)[pairs[, 2]],
    correlation = cor_matrix[pairs],
    stringsAsFactors = FALSE
  )
}


#' @title Create Predict Function for Various Model Types
#'
#' @description
#' Creates a standardized predict function for different model types.
#'
#' @param model A model object
#'
#' @return A predict function
#'
#' @keywords internal
create_predict_function <- function(model) {
  # mlr3 Learner
  if (inherits(model, "Learner")) {
    return(function(model, newdata) {
      task <- mlr3::as_task_classif(newdata, target = "dummy")
      pred <- model$predict(task)
      pred$prob[, 1]
    })
  }

  # caret model
  if (inherits(model, "train")) {
    return(function(model, newdata) {
      predict(model, newdata, type = "prob")[, 1]
    })
  }

  # Default: assume standard predict
  function(model, newdata) {
    pred <- predict(model, newdata, type = "response")
    if (is.matrix(pred)) pred[, 1] else pred
  }
}


#' @title Aggregate SHAP Values from List
#'
#' @description
#' Aggregates SHAP results into a matrix.
#'
#' @param shap_results List of SHAP data.frames
#' @param feature_names Character vector of feature names
#'
#' @return Matrix of SHAP values (samples x features)
#'
#' @keywords internal
aggregate_shap_values <- function(shap_results, feature_names) {
  n_samples <- length(shap_results)
  n_features <- length(feature_names)

  shap_matrix <- matrix(0, nrow = n_samples, ncol = n_features)
  colnames(shap_matrix) <- feature_names

  for (i in seq_len(n_samples)) {
    shap_df <- shap_results[[i]]

    # Extract contributions by variable
    if ("variable" %in% names(shap_df) && "contribution" %in% names(shap_df)) {
      for (j in seq_len(nrow(shap_df))) {
        var_name <- shap_df$variable[j]
        if (var_name %in% feature_names) {
          shap_matrix[i, var_name] <- shap_df$contribution[j]
        }
      }
    }
  }

  shap_matrix
}


#' @title Generate Correlation Warnings
#'
#' @description
#' Creates warning messages for correlated feature pairs.
#'
#' @param high_cor_pairs Data.frame of correlated pairs
#' @param feature_importance Named vector of feature importance
#' @param threshold Correlation threshold
#'
#' @return Data.frame with warnings
#'
#' @keywords internal
generate_correlation_warnings <- function(high_cor_pairs, feature_importance, threshold) {
  if (nrow(high_cor_pairs) == 0) {
    return(data.frame(
      feature1 = character(0),
      feature2 = character(0),
      correlation = numeric(0),
      importance1 = numeric(0),
      importance2 = numeric(0),
      warning = character(0),
      severity = character(0),
      stringsAsFactors = FALSE
    ))
  }

  # Add importance and generate warnings
  high_cor_pairs$importance1 <- feature_importance[high_cor_pairs$feature1]
  high_cor_pairs$importance2 <- feature_importance[high_cor_pairs$feature2]

  high_cor_pairs$warning <- sapply(seq_len(nrow(high_cor_pairs)), function(i) {
    f1 <- high_cor_pairs$feature1[i]
    f2 <- high_cor_pairs$feature2[i]
    cor <- high_cor_pairs$correlation[i]

    sprintf(
      "%s and %s have correlation %.2f. SHAP importance may be arbitrarily split between them.",
      f1, f2, cor
    )
  })

  # Severity based on correlation strength
  high_cor_pairs$severity <- ifelse(
    abs(high_cor_pairs$correlation) >= 0.9, "HIGH",
    ifelse(abs(high_cor_pairs$correlation) >= 0.8, "MEDIUM", "LOW")
  )

  # Sort by severity and correlation
  high_cor_pairs <- high_cor_pairs[order(
    factor(high_cor_pairs$severity, levels = c("HIGH", "MEDIUM", "LOW")),
    -abs(high_cor_pairs$correlation)
  ), ]

  high_cor_pairs
}


#' @title Print SHAP Warnings Report
#'
#' @description
#' Prints a formatted report of SHAP interpretation warnings.
#'
#' @param shap_result Result from compute_shap_with_warnings()
#'
#' @export
print_shap_warnings <- function(shap_result) {
  cat("SHAP Interpretation Report\n")
  cat("==========================\n\n")

  cat(sprintf("Correlation Threshold: %.2f\n", shap_result$threshold))
  cat(sprintf("Total Features: %d\n", ncol(shap_result$shap_values)))
  cat(sprintf("Reliable Features: %d\n", length(shap_result$reliable_features)))
  cat(sprintf("Correlated Features: %d\n\n", length(shap_result$correlated_features)))

  if (nrow(shap_result$warnings) > 0) {
    cat("CORRELATION WARNINGS:\n")
    cat("---------------------\n")

    for (i in seq_len(nrow(shap_result$warnings))) {
      w <- shap_result$warnings[i, ]
      severity_marker <- switch(w$severity,
        "HIGH" = "[!!!]",
        "MEDIUM" = "[!!]",
        "LOW" = "[!]"
      )
      cat(sprintf("%s %s\n", severity_marker, w$warning))
    }
    cat("\n")
  }

  if (length(shap_result$reliable_features) > 0) {
    cat("RELIABLE FEATURES (safe to interpret):\n")
    cat("--------------------------------------\n")

    # Show top reliable features by importance
    reliable_importance <- shap_result$feature_importance[shap_result$reliable_features]
    reliable_importance <- sort(reliable_importance, decreasing = TRUE)

    top_n <- min(10, length(reliable_importance))
    for (i in seq_len(top_n)) {
      cat(sprintf("  %d. %s (importance: %.4f)\n",
                  i, names(reliable_importance)[i], reliable_importance[i]))
    }
  }

  invisible(shap_result)
}


#' @title Get Reliable SHAP Features
#'
#' @description
#' Returns only the features whose SHAP values are reliable (no high correlations).
#'
#' @param shap_result Result from compute_shap_with_warnings()
#' @param top_n Number of top features to return (default: all)
#'
#' @return Named vector of feature importance for reliable features only
#'
#' @export
get_reliable_shap_features <- function(shap_result, top_n = NULL) {
  reliable_importance <- shap_result$feature_importance[shap_result$reliable_features]
  reliable_importance <- sort(reliable_importance, decreasing = TRUE)

  if (!is.null(top_n) && top_n < length(reliable_importance)) {
    reliable_importance <- reliable_importance[1:top_n]
  }

  reliable_importance
}
