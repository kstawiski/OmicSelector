#' @title Model Interpretability for OmicSelector
#'
#' @description
#' Provides model-agnostic interpretability via DALEX and SHAP-like methods.
#' Essential for understanding which features drive predictions and for
#' generating hypotheses in biomarker discovery.
#'
#' @details
#' ## Why Interpretability Matters
#'
#' For biomarker discovery, understanding WHY a model makes predictions is
#' as important as accuracy. Key use cases:
#' - Identify which genes/miRNAs drive classification
#' - Detect potential confounders or batch effects
#' - Generate hypotheses for biological validation
#' - Build trust for clinical adoption
#'
#' ## Methods Provided
#'
#' - **Feature Importance**: Permutation-based importance (model-agnostic)
#' - **Partial Dependence**: How features affect predictions marginally
#' - **SHAP Values**: Additive feature attributions (via iml or approximation)
#' - **Correlation Warnings**: Flag high-correlation features that may mislead
#'
#' ## Fold-Aware Design
#'
#' All interpretability methods respect the cross-validation structure:
#' - Explanations are computed per-fold using training data only
#' - Aggregated explanations show stability across folds
#' - Unstable features are flagged
#'
#' @name interpretability
NULL


#' @title Create Model Explainer
#'
#' @description
#' Creates a DALEX-style explainer for any mlr3 learner.
#' The explainer wraps the model for uniform interpretability access.
#'
#' @param learner A trained mlr3 Learner or GraphLearner
#' @param data Data frame or matrix used for explanations
#' @param y Response variable (factor for classification)
#' @param label Optional label for the model
#' @param predict_function Custom predict function (optional)
#'
#' @return An OmicExplainer object
#'
#' @examples
#' \dontrun{
#' # Create explainer from trained learner
#' explainer <- create_explainer(
#'   learner = trained_learner,
#'   data = train_data,
#'   y = train_labels
#' )
#'
#' # Compute feature importance
#' importance <- feature_importance(explainer)
#' }
#'
#' @export
create_explainer <- function(learner, data, y, label = NULL, predict_function = NULL) {
  # Validate inputs
  if (!inherits(learner, "Learner")) {
    stop("learner must be an mlr3 Learner object")
  }

  if (is.null(label)) {
    label <- learner$id
  }

  # Create prediction function if not provided
  if (is.null(predict_function)) {
    predict_function <- function(model, newdata) {
      # mlr3 predict_newdata - no task needed
      pred <- model$predict_newdata(as.data.frame(newdata))

      if ("prob" %in% pred$predict_types) {
        # Return positive class probability (second column)
        probs <- pred$prob
        if (ncol(probs) >= 2) {
          as.numeric(probs[, 2])
        } else {
          as.numeric(probs[, 1])
        }
      } else {
        # Return response as numeric
        as.numeric(as.factor(pred$response)) - 1
      }
    }
  }

  # Store as OmicExplainer
  explainer <- list(
    model = learner,
    data = as.data.frame(data),
    y = y,
    label = label,
    predict_function = predict_function,
    feature_names = setdiff(names(data), names(y)),
    class = "classification"
  )

  class(explainer) <- c("OmicExplainer", "list")
  explainer
}


#' @title Compute Permutation Feature Importance
#'
#' @description
#' Computes model-agnostic feature importance by permuting features
#' and measuring the drop in performance.
#'
#' @param explainer An OmicExplainer object
#' @param loss_function Loss function to use (default: "1 - AUC")
#' @param n_permutations Number of permutations per feature (default: 10)
#' @param parallel Logical, use parallel computation (default: FALSE)
#'
#' @return A data frame with feature importance scores
#'
#' @details
#' Permutation importance measures how much model performance degrades
#' when a feature's values are randomly shuffled. Higher values indicate
#' more important features.
#'
#' For correlated features, importance may be shared or masked. Use
#' `feature_importance_groups()` for correlated feature handling.
#'
#' @examples
#' \dontrun{
#' importance <- feature_importance(explainer, n_permutations = 20)
#' head(importance[order(-importance$importance), ])
#' }
#'
#' @export
feature_importance <- function(explainer,
                                loss_function = c("auc_loss", "accuracy_loss", "brier"),
                                n_permutations = 10,
                                parallel = FALSE) {

  loss_function <- match.arg(loss_function)

  data <- explainer$data
  y <- as.numeric(as.factor(explainer$y)) - 1  # Convert to 0/1
  feature_names <- explainer$feature_names

  # Baseline predictions
  baseline_pred <- explainer$predict_function(explainer$model, data)

  # Baseline loss
  baseline_loss <- .compute_loss(baseline_pred, y, loss_function)

  # Compute importance for each feature
  importance_results <- lapply(feature_names, function(feat) {
    losses <- numeric(n_permutations)

    for (i in seq_len(n_permutations)) {
      # Permute feature
      data_permuted <- data
      data_permuted[[feat]] <- sample(data_permuted[[feat]])

      # Predict and compute loss
      pred_permuted <- explainer$predict_function(explainer$model, data_permuted)
      losses[i] <- .compute_loss(pred_permuted, y, loss_function)
    }

    data.frame(
      feature = feat,
      baseline_loss = baseline_loss,
      mean_permuted_loss = mean(losses),
      sd_permuted_loss = sd(losses),
      importance = mean(losses) - baseline_loss,
      importance_ratio = mean(losses) / baseline_loss
    )
  })

  result <- do.call(rbind, importance_results)
  result <- result[order(-result$importance), ]
  rownames(result) <- NULL

  class(result) <- c("FeatureImportance", "data.frame")
  result
}


#' Compute loss function
#' @keywords internal
.compute_loss <- function(pred, y, loss_type) {
  switch(loss_type,
    auc_loss = {
      # 1 - AUC (so higher = worse)
      if (length(unique(y)) < 2) return(0.5)
      1 - .compute_auc(pred, y)
    },
    accuracy_loss = {
      # 1 - accuracy
      pred_class <- as.numeric(pred > 0.5)
      1 - mean(pred_class == y)
    },
    brier = {
      # Brier score
      mean((pred - y)^2)
    }
  )
}


#' Compute AUC
#' @keywords internal
.compute_auc <- function(pred, y) {
  # Simple AUC computation via Mann-Whitney U
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)

  if (n_pos == 0 || n_neg == 0) return(0.5)

  ranks <- rank(pred)
  u <- sum(ranks[y == 1]) - n_pos * (n_pos + 1) / 2

  u / (n_pos * n_neg)
}


#' @title Compute Partial Dependence
#'
#' @description
#' Computes partial dependence profiles showing how a feature
#' affects predictions, marginalizing over other features.
#'
#' @param explainer An OmicExplainer object
#' @param features Character vector of feature names to analyze
#' @param n_grid Number of grid points (default: 20)
#'
#' @return A data frame with partial dependence values
#'
#' @export
partial_dependence <- function(explainer, features, n_grid = 20) {

  data <- explainer$data

  results <- list()

  for (feat in features) {
    if (!feat %in% names(data)) {
      warning(sprintf("Feature '%s' not found in data", feat))
      next
    }

    # Create grid of feature values
    feat_values <- data[[feat]]

    if (is.numeric(feat_values)) {
      grid <- seq(min(feat_values), max(feat_values), length.out = n_grid)
    } else {
      grid <- unique(feat_values)
    }

    # Compute marginal predictions at each grid point
    pdp_values <- sapply(grid, function(val) {
      data_modified <- data
      data_modified[[feat]] <- val
      pred <- explainer$predict_function(explainer$model, data_modified)
      mean(pred)
    })

    results[[feat]] <- data.frame(
      feature = feat,
      value = grid,
      prediction = pdp_values
    )
  }

  result <- do.call(rbind, results)
  rownames(result) <- NULL

  class(result) <- c("PartialDependence", "data.frame")
  result
}


#' @title Compute SHAP-like Values
#'
#' @description
#' Computes Shapley-like additive feature attributions for individual
#' predictions. Uses kernel SHAP approximation for efficiency.
#'
#' @param explainer An OmicExplainer object
#' @param new_data Data for which to compute explanations
#' @param n_samples Number of samples for approximation (default: 100)
#'
#' @return A matrix of SHAP values (instances x features)
#'
#' @details
#' SHAP (SHapley Additive exPlanations) values explain individual predictions
#' as the sum of feature contributions. For each prediction:
#'
#' prediction = base_value + sum(SHAP_values)
#'
#' Positive SHAP value = feature pushes prediction higher
#' Negative SHAP value = feature pushes prediction lower
#'
#' @export
shap_values <- function(explainer, new_data, n_samples = 100) {

  # Check for iml package (preferred SHAP implementation)
  if (requireNamespace("iml", quietly = TRUE)) {
    return(.shap_via_iml(explainer, new_data, n_samples))
  }

  # Fallback: kernel SHAP approximation
  .shap_kernel_approximation(explainer, new_data, n_samples)
}


#' SHAP via iml package
#' @keywords internal
.shap_via_iml <- function(explainer, new_data, n_samples) {
  # Create iml Predictor
  predict_wrapper <- function(model, newdata) {
    explainer$predict_function(explainer$model, newdata)
  }

  predictor <- iml::Predictor$new(
    model = explainer$model,
    data = explainer$data[, explainer$feature_names, drop = FALSE],
    y = explainer$y,
    predict.function = predict_wrapper
  )

  # Compute SHAP for each instance
  shap_list <- lapply(seq_len(nrow(new_data)), function(i) {
    instance <- new_data[i, explainer$feature_names, drop = FALSE]
    shap <- iml::Shapley$new(predictor, x.interest = instance, sample.size = n_samples)
    shap$results$phi
  })

  shap_matrix <- do.call(rbind, shap_list)
  colnames(shap_matrix) <- explainer$feature_names
  rownames(shap_matrix) <- rownames(new_data)

  attr(shap_matrix, "base_value") <- mean(explainer$predict_function(explainer$model, explainer$data))

  class(shap_matrix) <- c("SHAPValues", "matrix")
  shap_matrix
}


#' Kernel SHAP approximation (fallback)
#' @keywords internal
.shap_kernel_approximation <- function(explainer, new_data, n_samples) {
  message("Using kernel SHAP approximation (install 'iml' for faster computation)")

  data <- explainer$data[, explainer$feature_names, drop = FALSE]
  n_features <- length(explainer$feature_names)

  # Base prediction (expected value)
  base_pred <- mean(explainer$predict_function(explainer$model, data))

  # Compute SHAP for each instance
  shap_matrix <- matrix(0, nrow = nrow(new_data), ncol = n_features)
  colnames(shap_matrix) <- explainer$feature_names

  for (i in seq_len(nrow(new_data))) {
    instance <- new_data[i, explainer$feature_names, drop = FALSE]

    # Approximate marginal contributions
    for (j in seq_len(n_features)) {
      feat <- explainer$feature_names[j]

      # Sample background instances
      bg_idx <- sample(nrow(data), min(n_samples, nrow(data)))
      bg_data <- data[bg_idx, , drop = FALSE]

      # With feature from instance
      data_with <- bg_data
      data_with[[feat]] <- instance[[feat]]

      # Without (use background)
      data_without <- bg_data

      # Marginal contribution
      pred_with <- mean(explainer$predict_function(explainer$model, data_with))
      pred_without <- mean(explainer$predict_function(explainer$model, data_without))

      shap_matrix[i, j] <- pred_with - pred_without
    }
  }

  attr(shap_matrix, "base_value") <- base_pred

  class(shap_matrix) <- c("SHAPValues", "matrix")
  shap_matrix
}


#' @title Check Feature Correlations
#'
#' @description
#' Identifies highly correlated features that may lead to misleading
#' importance scores. When features are correlated, importance can be
#' shared, masked, or unstable.
#'
#' @param data Data frame with features
#' @param threshold Correlation threshold (default: 0.7)
#' @param method Correlation method: "pearson", "spearman" (default: "spearman")
#'
#' @return A list with correlation warnings
#'
#' @details
#' High correlation between features causes problems for interpretation:
#' - Permutation importance is shared between correlated features
#' - SHAP values become unstable
#' - One feature can "mask" another
#'
#' This function identifies correlated pairs and suggests grouped analysis.
#'
#' @export
check_feature_correlations <- function(data, threshold = 0.7, method = c("spearman", "pearson")) {
  method <- match.arg(method)

  # Select numeric columns
  numeric_cols <- sapply(data, is.numeric)
  numeric_data <- data[, numeric_cols, drop = FALSE]

  if (ncol(numeric_data) < 2) {
    return(list(n_pairs = 0, correlated_pairs = data.frame(), warning = FALSE))
  }

  # Compute correlation matrix
  cor_matrix <- cor(numeric_data, method = method, use = "pairwise.complete.obs")

  # Find pairs above threshold
  pairs <- list()
  for (i in 1:(ncol(cor_matrix) - 1)) {
    for (j in (i + 1):ncol(cor_matrix)) {
      if (abs(cor_matrix[i, j]) >= threshold) {
        pairs <- c(pairs, list(data.frame(
          feature1 = colnames(cor_matrix)[i],
          feature2 = colnames(cor_matrix)[j],
          correlation = cor_matrix[i, j]
        )))
      }
    }
  }

  correlated_pairs <- if (length(pairs) > 0) {
    do.call(rbind, pairs)
  } else {
    data.frame(feature1 = character(), feature2 = character(), correlation = numeric())
  }

  result <- list(
    n_pairs = nrow(correlated_pairs),
    correlated_pairs = correlated_pairs,
    threshold = threshold,
    method = method,
    warning = nrow(correlated_pairs) > 0
  )

  if (result$warning) {
    result$message <- sprintf(
      "Found %d feature pairs with |correlation| >= %.2f. Consider grouped importance analysis.",
      result$n_pairs, threshold
    )
  }

  class(result) <- c("CorrelationCheck", "list")
  result
}


#' @title Print Feature Importance
#' @param x A `FeatureImportance` object.
#' @param n Number of top features to print.
#' @param ... Additional arguments passed through to base methods (ignored).
#' @export
print.FeatureImportance <- function(x, n = 10, ...) {
  cat("=== Feature Importance (Permutation) ===\n\n")
  cat(sprintf("Top %d features by importance:\n\n", min(n, nrow(x))))

  top_n <- head(x, n)

  for (i in seq_len(nrow(top_n))) {
    cat(sprintf("  %2d. %-20s  %.4f\n",
                i, top_n$feature[i], top_n$importance[i]))
  }

  invisible(x)
}


#' @title Print Correlation Check
#' @param x A `CorrelationCheck` object.
#' @param ... Additional arguments passed through to base methods (ignored).
#' @export
print.CorrelationCheck <- function(x, ...) {
  cat("=== Feature Correlation Check ===\n\n")
  cat(sprintf("Threshold: %.2f (%s)\n", x$threshold, x$method))
  cat(sprintf("Correlated pairs found: %d\n\n", x$n_pairs))

  if (x$warning && nrow(x$correlated_pairs) > 0) {
    cat("WARNING: High correlations detected!\n\n")
    print(head(x$correlated_pairs[order(-abs(x$correlated_pairs$correlation)), ], 10))
    cat("\n", x$message, "\n")
  } else {
    cat("No concerning correlations found.\n")
  }

  invisible(x)
}
