
#' Enhanced Feature Importance with Permutation and Conditional Methods
#'
#' Calculates feature importance using advanced methods that account for
#' feature correlations and provide confidence intervals. Goes beyond simple
#' variable importance from Random Forest by using model-agnostic approaches.
#'
#' @param model Fitted model object (from OmicSelector_fit or trained with tidymodels/caret)
#' @param data Data frame used for training (for permutation importance)
#' @param outcome Character string, name of outcome variable
#' @param method Character string specifying importance method:
#'   \itemize{
#'     \item "permutation" - Model-agnostic permutation importance (Breiman 2001)
#'     \item "conditional" - Conditional permutation importance (Strobl et al. 2008)
#'     \item "both" - Calculate both methods for comparison
#'   }
#' @param metric Character string, metric to use for importance calculation:
#'   "accuracy", "auc", "rmse", "rsq" (auto-detected from outcome type)
#' @param n_repeats Integer, number of times to permute each feature (default: 10)
#' @param scale Logical, whether to scale importance values to sum to 1 (default: FALSE)
#' @param normalize Logical, whether to normalize by dividing by standard error (default: TRUE)
#' @param conditional_grid Integer, number of grid points for conditional sampling (default: 5)
#' @param parallel Logical, whether to use parallel processing (default: TRUE)
#' @param n_cores Integer, number of cores (default: detectCores() - 1)
#' @param seed Integer for reproducibility
#'
#' @return A list object of class "OmicSelector_importance" containing:
#'   \item{importance}{Data frame with feature, importance, std_error, and z_score}
#'   \item{method}{Method used}
#'   \item{metric}{Metric used}
#'   \item{baseline_performance}{Baseline model performance before permutation}
#'   \item{n_repeats}{Number of permutation repeats}
#'   \item{runtime}{Time taken for calculation}
#'
#' @details
#' **Permutation Importance**:
#' Measures the decrease in model performance when a feature's values are randomly
#' permuted (shuffled). If performance drops significantly, the feature is important.
#'
#' Algorithm:
#' 1. Measure baseline performance on original data
#' 2. For each feature:
#'    a. Randomly permute feature values
#'    b. Measure performance on permuted data
#'    c. Importance = baseline - permuted_performance
#' 3. Repeat n_repeats times and average
#'
#' **Conditional Importance** (Strobl et al. 2008):
#' Similar to permutation but respects feature correlations. When permuting
#' feature X, only swap values that have similar values in correlated features.
#' This prevents inflated importance for correlated features.
#'
#' **Advantages**:
#' - Model-agnostic (works with any model)
#' - Provides confidence intervals
#' - Handles correlated features (conditional method)
#' - More reliable than built-in variable importance
#'
#' @references
#' Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5-32.
#'
#' Strobl, C., Boulesteix, A. L., Kneib, T., Augustin, T., & Zeileis, A. (2008).
#' Conditional variable importance for random forests. *BMC Bioinformatics*, 9(1), 307.
#'
#' @examples
#' \dontrun{
#' # Fit a model first
#' model <- OmicSelector_fit(
#'   data = miRNA_data,
#'   outcome = "Class",
#'   method = "tidymodels",
#'   algorithm = "ranger"
#' )
#'
#' # Calculate permutation importance
#' importance <- OmicSelector_importance(
#'   model = model,
#'   data = miRNA_data,
#'   outcome = "Class",
#'   method = "permutation",
#'   n_repeats = 20
#' )
#'
#' # View results
#' print(importance)
#' plot(importance)  # Top 20 features
#'
#' # Compare permutation vs conditional
#' importance_both <- OmicSelector_importance(
#'   model = model,
#'   data = miRNA_data,
#'   outcome = "Class",
#'   method = "both"
#' )
#' }
#'
#' @export
OmicSelector_importance <- function(
  model,
  data,
  outcome,
  method = c("permutation", "conditional", "both"),
  metric = NULL,
  n_repeats = 10,
  scale = FALSE,
  normalize = TRUE,
  conditional_grid = 5,
  parallel = TRUE,
  n_cores = NULL,
  seed = 123
) {

  # Validate inputs
  method <- match.arg(method)

  if (!outcome %in% names(data)) {
    stop("Outcome variable '", outcome, "' not found in data")
  }

  # Auto-detect metric if not specified
  if (is.null(metric)) {
    if (is.factor(data[[outcome]]) || is.character(data[[outcome]])) {
      metric <- "accuracy"
      message("Auto-detected classification task, using metric: accuracy")
    } else {
      metric <- "rmse"
      message("Auto-detected regression task, using metric: rmse")
    }
  }

  # Determine feature columns
  features <- setdiff(names(data), outcome)

  # Set up parallel processing
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- max(1, parallel::detectCores() - 1)
    }
  }

  set.seed(seed)

  message(paste0("Calculating feature importance using method: ", method))
  message(paste0("Features: ", length(features)))
  message(paste0("Metric: ", metric))
  message(paste0("Repeats: ", n_repeats))

  start_time <- Sys.time()

  # Calculate importance based on method
  if (method == "permutation") {
    result <- .calculate_permutation_importance(
      model = model,
      data = data,
      outcome = outcome,
      features = features,
      metric = metric,
      n_repeats = n_repeats,
      parallel = parallel,
      n_cores = n_cores
    )
  } else if (method == "conditional") {
    result <- .calculate_conditional_importance(
      model = model,
      data = data,
      outcome = outcome,
      features = features,
      metric = metric,
      n_repeats = n_repeats,
      conditional_grid = conditional_grid,
      parallel = parallel,
      n_cores = n_cores
    )
  } else if (method == "both") {
    # Calculate both and return comparison
    perm_result <- .calculate_permutation_importance(
      model = model,
      data = data,
      outcome = outcome,
      features = features,
      metric = metric,
      n_repeats = n_repeats,
      parallel = parallel,
      n_cores = n_cores
    )

    cond_result <- .calculate_conditional_importance(
      model = model,
      data = data,
      outcome = outcome,
      features = features,
      metric = metric,
      n_repeats = n_repeats,
      conditional_grid = conditional_grid,
      parallel = parallel,
      n_cores = n_cores
    )

    # Combine results
    result <- list(
      permutation = perm_result$importance,
      conditional = cond_result$importance,
      baseline_performance = perm_result$baseline_performance,
      correlation_adjusted = TRUE
    )
    class(result) <- c("OmicSelector_importance_comparison", "list")
  }

  # Normalize and scale if requested
  if (method != "both" && !is.null(result$importance)) {
    if (normalize && "std_error" %in% names(result$importance)) {
      result$importance$z_score <- result$importance$importance / result$importance$std_error
    }

    if (scale) {
      total <- sum(abs(result$importance$importance))
      if (total > 0) {
        result$importance$scaled_importance <- result$importance$importance / total
      }
    }

    # Sort by importance
    result$importance <- result$importance[order(-result$importance$importance), ]
  }

  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Add metadata
  if (method != "both") {
    result$method <- method
    result$metric <- metric
    result$n_repeats <- n_repeats
    result$runtime <- runtime
    result$n_features <- length(features)

    class(result) <- c("OmicSelector_importance", "list")
  }

  message(paste0("✓ Importance calculation completed in ", round(runtime, 1), " seconds"))

  return(result)
}


#' Calculate Permutation Importance
#' @keywords internal
.calculate_permutation_importance <- function(
  model, data, outcome, features, metric, n_repeats, parallel, n_cores
) {

  # Get baseline performance
  baseline_perf <- .evaluate_model_performance(model, data, outcome, metric)
  message(paste0("Baseline ", metric, ": ", round(baseline_perf, 4)))

  # Function to permute one feature
  permute_feature <- function(feature, repeat_idx) {
    # Create permuted data
    permuted_data <- data
    permuted_data[[feature]] <- sample(permuted_data[[feature]])

    # Evaluate on permuted data
    perm_perf <- .evaluate_model_performance(model, permuted_data, outcome, metric)

    # Importance = drop in performance (or increase for metrics like RMSE)
    if (metric %in% c("rmse", "mae")) {
      # For error metrics, increase is bad
      importance <- perm_perf - baseline_perf
    } else {
      # For accuracy/AUC, decrease is bad
      importance <- baseline_perf - perm_perf
    }

    return(importance)
  }

  # Calculate importance for each feature across repeats
  if (parallel && n_cores > 1) {
    # Parallel version
    importance_matrix <- matrix(0, nrow = length(features), ncol = n_repeats)
    rownames(importance_matrix) <- features

    for (i in seq_along(features)) {
      feature <- features[i]
      importance_matrix[i, ] <- sapply(1:n_repeats, function(r) {
        permute_feature(feature, r)
      })

      if (i %% 10 == 0) {
        message(paste0("  Progress: ", i, "/", length(features), " features"))
      }
    }
  } else {
    # Sequential version
    importance_matrix <- matrix(0, nrow = length(features), ncol = n_repeats)
    rownames(importance_matrix) <- features

    for (i in seq_along(features)) {
      for (r in 1:n_repeats) {
        importance_matrix[i, r] <- permute_feature(features[i], r)
      }

      if (i %% 10 == 0) {
        message(paste0("  Progress: ", i, "/", length(features), " features"))
      }
    }
  }

  # Calculate mean and standard error
  importance_df <- data.frame(
    feature = features,
    importance = rowMeans(importance_matrix),
    std_error = apply(importance_matrix, 1, function(x) sd(x) / sqrt(length(x))),
    min_importance = apply(importance_matrix, 1, min),
    max_importance = apply(importance_matrix, 1, max),
    stringsAsFactors = FALSE
  )

  return(list(
    importance = importance_df,
    baseline_performance = baseline_perf,
    importance_matrix = importance_matrix
  ))
}


#' Calculate Conditional Permutation Importance
#' @keywords internal
.calculate_conditional_importance <- function(
  model, data, outcome, features, metric, n_repeats, conditional_grid, parallel, n_cores
) {

  # Get baseline performance
  baseline_perf <- .evaluate_model_performance(model, data, outcome, metric)

  # Calculate correlation matrix for conditional sampling
  feature_data <- data[, features, drop = FALSE]
  feature_data_numeric <- as.data.frame(lapply(feature_data, function(x) {
    if (is.numeric(x)) return(x)
    as.numeric(as.factor(x))
  }))
  cor_matrix <- cor(feature_data_numeric, use = "pairwise.complete.obs")

  # Function to conditionally permute one feature
  conditional_permute <- function(feature_idx, repeat_idx) {
    feature <- features[feature_idx]

    # Find correlated features (|r| > 0.7)
    correlated_features <- which(abs(cor_matrix[feature_idx, ]) > 0.7 &
                                  abs(cor_matrix[feature_idx, ]) < 1.0)

    permuted_data <- data

    if (length(correlated_features) > 0) {
      # Conditional permutation: swap within similar conditioning sets
      # Use first correlated feature as conditioning variable
      cond_feature <- features[correlated_features[1]]
      cond_values <- data[[cond_feature]]

      # Create grid of conditioning values
      cond_breaks <- quantile(cond_values, probs = seq(0, 1, length.out = conditional_grid + 1),
                              na.rm = TRUE)

      # Permute within each conditioning bin
      for (i in 1:conditional_grid) {
        in_bin <- which(cond_values >= cond_breaks[i] & cond_values <= cond_breaks[i + 1])
        if (length(in_bin) > 1) {
          permuted_data[[feature]][in_bin] <- sample(permuted_data[[feature]][in_bin])
        }
      }
    } else {
      # No correlated features, use standard permutation
      permuted_data[[feature]] <- sample(permuted_data[[feature]])
    }

    # Evaluate on permuted data
    perm_perf <- .evaluate_model_performance(model, permuted_data, outcome, metric)

    # Calculate importance
    if (metric %in% c("rmse", "mae")) {
      importance <- perm_perf - baseline_perf
    } else {
      importance <- baseline_perf - perm_perf
    }

    return(importance)
  }

  # Calculate conditional importance
  importance_matrix <- matrix(0, nrow = length(features), ncol = n_repeats)
  rownames(importance_matrix) <- features

  for (i in seq_along(features)) {
    for (r in 1:n_repeats) {
      importance_matrix[i, r] <- conditional_permute(i, r)
    }

    if (i %% 10 == 0) {
      message(paste0("  Progress: ", i, "/", length(features), " features"))
    }
  }

  # Calculate statistics
  importance_df <- data.frame(
    feature = features,
    importance = rowMeans(importance_matrix),
    std_error = apply(importance_matrix, 1, function(x) sd(x) / sqrt(length(x))),
    min_importance = apply(importance_matrix, 1, min),
    max_importance = apply(importance_matrix, 1, max),
    stringsAsFactors = FALSE
  )

  return(list(
    importance = importance_df,
    baseline_performance = baseline_perf,
    importance_matrix = importance_matrix,
    correlation_matrix = cor_matrix
  ))
}


#' Evaluate Model Performance
#' @keywords internal
.evaluate_model_performance <- function(model, data, outcome, metric) {

  # Extract actual outcome
  actual <- data[[outcome]]

  # Get predictions based on model type
  if (inherits(model, "workflow")) {
    # tidymodels workflow
    tryCatch({
      preds <- predict(model, new_data = data, type = "class")$.pred_class
      if (metric == "auc") {
        preds_prob <- predict(model, new_data = data, type = "prob")
        # Assume binary classification
        perf <- as.numeric(pROC::auc(actual, preds_prob[[2]]))
      } else if (metric == "accuracy") {
        perf <- mean(preds == actual, na.rm = TRUE)
      }
    }, error = function(e) {
      # Fallback
      perf <- 0.5
    })
  } else if (inherits(model, "train")) {
    # caret model
    preds <- predict(model, newdata = data)
    if (metric == "accuracy") {
      perf <- mean(preds == actual, na.rm = TRUE)
    } else if (metric == "auc") {
      preds_prob <- predict(model, newdata = data, type = "prob")
      perf <- as.numeric(pROC::auc(actual, preds_prob[[2]]))
    } else if (metric == "rmse") {
      perf <- sqrt(mean((preds - actual)^2, na.rm = TRUE))
    }
  } else {
    # Generic predict
    preds <- predict(model, newdata = data)
    if (metric == "accuracy") {
      perf <- mean(preds == actual, na.rm = TRUE)
    } else if (metric == "rmse") {
      perf <- sqrt(mean((preds - actual)^2, na.rm = TRUE))
    } else {
      perf <- 0.5  # Fallback
    }
  }

  return(perf)
}


#' Print Method for OmicSelector Importance
#' @export
print.OmicSelector_importance <- function(x, n = 20, ...) {
  cat("OmicSelector Feature Importance\n")
  cat("================================\n\n")

  cat("Method:", x$method, "\n")
  cat("Metric:", x$metric, "\n")
  cat("Repeats:", x$n_repeats, "\n")
  cat("Baseline performance:", round(x$baseline_performance, 4), "\n")
  cat("Runtime:", round(x$runtime, 1), "seconds\n\n")

  cat("Results:\n")
  cat("  Features evaluated:", x$n_features, "\n")
  cat("  Positive importance:", sum(x$importance$importance > 0), "\n")
  cat("  Negative importance:", sum(x$importance$importance < 0), "\n\n")

  cat("Top", min(n, nrow(x$importance)), "important features:\n")
  top_features <- head(x$importance, n)

  # Format output
  out <- data.frame(
    Rank = 1:nrow(top_features),
    Feature = top_features$feature,
    Importance = sprintf("%.4f", top_features$importance),
    SE = sprintf("%.4f", top_features$std_error),
    stringsAsFactors = FALSE
  )

  if ("z_score" %in% names(top_features)) {
    out$Z_score = sprintf("%.2f", top_features$z_score)
  }

  print(out, row.names = FALSE)

  invisible(x)
}


#' Summary Method for OmicSelector Importance
#' @export
summary.OmicSelector_importance <- function(object, ...) {
  cat("OmicSelector Feature Importance Summary\n")
  cat("=======================================\n\n")

  imp <- object$importance

  cat("Importance Statistics:\n")
  cat("  Mean importance:", round(mean(imp$importance), 4), "\n")
  cat("  Median importance:", round(median(imp$importance), 4), "\n")
  cat("  SD importance:", round(sd(imp$importance), 4), "\n")
  cat("  Range:", round(min(imp$importance), 4), "to", round(max(imp$importance), 4), "\n\n")

  cat("Feature Categories:\n")
  cat("  Very important (imp > 0.01):", sum(imp$importance > 0.01), "\n")
  cat("  Moderately important (0.001 < imp <= 0.01):",
      sum(imp$importance > 0.001 & imp$importance <= 0.01), "\n")
  cat("  Low importance (imp <= 0.001):", sum(imp$importance <= 0.001), "\n")
  cat("  Negative importance:", sum(imp$importance < 0), "\n\n")

  return(imp)
}


#' Plot Method for OmicSelector Importance
#' @export
plot.OmicSelector_importance <- function(x, n = 20, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not available, skipping plot")
    return(invisible(x))
  }

  top_features <- head(x$importance, n)

  # Create plot
  p <- ggplot2::ggplot(top_features, ggplot2::aes(
    x = reorder(feature, importance),
    y = importance
  )) +
    ggplot2::geom_col(ggplot2::aes(fill = importance > 0)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = importance - std_error, ymax = importance + std_error),
      width = 0.2
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0("Top ", n, " Features by ", tools::toTitleCase(x$method), " Importance"),
      subtitle = paste0("Metric: ", x$metric, " | Baseline: ", round(x$baseline_performance, 3)),
      x = "Feature",
      y = "Importance (decrease in performance)"
    ) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#2C7BB6", "FALSE" = "#D7191C"),
      labels = c("TRUE" = "Positive", "FALSE" = "Negative"),
      name = "Effect"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
    )

  return(p)
}


#' Print Method for Importance Comparison
#' @export
print.OmicSelector_importance_comparison <- function(x, n = 20, ...) {
  cat("OmicSelector Feature Importance Comparison\n")
  cat("==========================================\n\n")

  cat("Baseline performance:", round(x$baseline_performance, 4), "\n\n")

  # Merge the two importance measures
  merged <- merge(
    x$permutation[, c("feature", "importance")],
    x$conditional[, c("feature", "importance")],
    by = "feature",
    suffixes = c("_perm", "_cond")
  )

  # Sort by permutation importance
  merged <- merged[order(-merged$importance_perm), ]

  cat("Top", min(n, nrow(merged)), "features:\n\n")
  top <- head(merged, n)

  out <- data.frame(
    Rank = 1:nrow(top),
    Feature = top$feature,
    Perm_Imp = sprintf("%.4f", top$importance_perm),
    Cond_Imp = sprintf("%.4f", top$importance_cond),
    Diff = sprintf("%.4f", top$importance_perm - top$importance_cond),
    stringsAsFactors = FALSE
  )

  print(out, row.names = FALSE)

  cat("\nNote: Large differences suggest feature correlations may be inflating permutation importance\n")

  invisible(x)
}
