#' @title Auto XAI: Automatic Explainability for Biomarker Models
#'
#' @description
#' Provides automatic interpretability outputs for mlr3 models using DALEX.
#' Designed for high-dimensional omics data with correlation-aware warnings.
#'
#' @details
#' Key features:
#' - mlr3 to DALEX adapter for probabilistic classifiers
#' - Permutation-based feature importance
#' - Partial Dependence Plots (PDPs) for top features
#' - SHAP values for individual predictions
#' - Correlation diagnostics and warnings
#' - Automated report generation
#'
#' @name AutoXAI
NULL


#' @title Create DALEX Explainer from mlr3 Learner
#'
#' @description
#' Converts a trained mlr3 learner into a DALEX explainer for interpretability.
#'
#' @param learner A trained mlr3 Learner (must have predict_type = "prob")
#' @param task The mlr3 task used for training
#' @param data Optional: explicit data for explainer (if NULL, uses task data)
#' @param label Label for the explainer (default: learner ID)
#' @param features Optional: subset of features to include (speeds up explanations)
#'
#' @return A DALEX explainer object
#'
#' @examples
#' \dontrun{
#' learner <- lrn("classif.ranger", predict_type = "prob")
#' learner$train(task)
#' explainer <- xai_explainer_mlr3(learner, task)
#' }
#'
#' @export
xai_explainer_mlr3 <- function(learner,
                                task,
                                data = NULL,
                                label = NULL,
                                features = NULL) {

  if (!requireNamespace("DALEX", quietly = TRUE)) {
    stop("Package 'DALEX' is required. Install with: install.packages('DALEX')",
         call. = FALSE)
  }

  # Validate learner
  if (!inherits(learner, "Learner")) {
    stop("learner must be an mlr3 Learner object", call. = FALSE)
  }

  if (learner$predict_type != "prob") {
    stop("Learner must have predict_type = 'prob' for XAI", call. = FALSE)
  }

  if (is.null(learner$model)) {
    stop("Learner must be trained before creating explainer", call. = FALSE)
  }

  # Get data
  if (is.null(data)) {
    data <- as.data.frame(task$data())
  }

  # Get target
  target_name <- task$target_names
  positive_class <- task$positive

  # Extract features (use learner's feature names to ensure consistency)
  learner_features <- learner$state$feature_names
  X <- data[, learner_features, drop = FALSE]
  y <- data[[target_name]]

  # Subset to specific features if requested (for speed)
  # Note: We still pass full feature set to predict, but limit analysis
  analysis_features <- learner_features
  if (!is.null(features)) {
    analysis_features <- intersect(features, learner_features)
    if (length(analysis_features) == 0) {
      analysis_features <- learner_features
    }
  }

  # Convert target to numeric (1 = positive class, 0 = negative)
  y_numeric <- as.numeric(y == positive_class)

  # Create predict function for DALEX
  # Returns probability of positive class
  predict_fun <- function(model, newdata) {
    # Predict using mlr3 learner
    pred <- learner$predict_newdata(newdata)
    probs <- pred$prob

    # Return probability of positive class
    if (positive_class %in% colnames(probs)) {
      as.numeric(probs[, positive_class])
    } else {
      as.numeric(probs[, 2])  # Fallback to second column
    }
  }

  # Set label
  if (is.null(label)) {
    label <- learner$id
  }

  # Create explainer with full feature set
  explainer <- DALEX::explain(
    model = learner,
    data = X,
    y = y_numeric,
    predict_function = predict_fun,
    label = label,
    verbose = FALSE
  )

  # Store analysis features as attribute
  attr(explainer, "analysis_features") <- analysis_features

  explainer
}


#' @title Compute Correlation Diagnostics for Features
#'
#' @description
#' Identifies highly correlated feature pairs that may cause misleading
#' SHAP values or permutation importance.
#'
#' @param data Data frame of features
#' @param features Optional: subset of features to analyze
#' @param method Correlation method ("spearman", "pearson")
#' @param threshold Correlation threshold for warnings (default: 0.7)
#'
#' @return A list with correlation matrix, clusters, and warnings
#'
#' @export
xai_correlations <- function(data,
                              features = NULL,
                              method = "spearman",
                              threshold = 0.7) {

  # Subset features
  if (!is.null(features)) {
    data <- data[, intersect(features, names(data)), drop = FALSE]
  }

  # Remove non-numeric columns
  numeric_cols <- sapply(data, is.numeric)
  data <- data[, numeric_cols, drop = FALSE]

  if (ncol(data) < 2) {
    return(list(
      matrix = NULL,
      high_cor_pairs = data.frame(),
      clusters = list(),
      n_warnings = 0,
      threshold = threshold
    ))
  }

  # Compute correlation matrix
  cor_matrix <- cor(data, use = "pairwise.complete.obs", method = method)

  # Find high correlation pairs
  high_cor_pairs <- data.frame(
    feature1 = character(0),
    feature2 = character(0),
    correlation = numeric(0),
    stringsAsFactors = FALSE
  )

  n <- ncol(cor_matrix)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (abs(cor_matrix[i, j]) >= threshold) {
        high_cor_pairs <- rbind(high_cor_pairs, data.frame(
          feature1 = colnames(cor_matrix)[i],
          feature2 = colnames(cor_matrix)[j],
          correlation = cor_matrix[i, j],
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Sort by absolute correlation
  if (nrow(high_cor_pairs) > 0) {
    high_cor_pairs <- high_cor_pairs[order(-abs(high_cor_pairs$correlation)), ]
  }

  # Identify correlation clusters using simple connected components
  clusters <- identify_correlation_clusters(high_cor_pairs, names(data))

  list(
    matrix = cor_matrix,
    high_cor_pairs = high_cor_pairs,
    clusters = clusters,
    n_warnings = nrow(high_cor_pairs),
    threshold = threshold,
    method = method
  )
}


#' @title Identify Correlation Clusters
#'
#' @description
#' Groups highly correlated features into clusters.
#'
#' @param high_cor_pairs Data frame from xai_correlations
#' @param all_features All feature names
#'
#' @return List of character vectors (clusters)
#'
#' @keywords internal
identify_correlation_clusters <- function(high_cor_pairs, all_features) {

  if (nrow(high_cor_pairs) == 0) {
    return(list())
  }

  # Build adjacency list
  adj <- list()
  for (i in seq_len(nrow(high_cor_pairs))) {
    f1 <- high_cor_pairs$feature1[i]
    f2 <- high_cor_pairs$feature2[i]

    adj[[f1]] <- unique(c(adj[[f1]], f2))
    adj[[f2]] <- unique(c(adj[[f2]], f1))
  }

  # Find connected components using BFS
  visited <- character(0)
  clusters <- list()

  for (feature in names(adj)) {
    if (!(feature %in% visited)) {
      # BFS to find component
      queue <- feature
      component <- character(0)

      while (length(queue) > 0) {
        current <- queue[1]
        queue <- queue[-1]

        if (!(current %in% visited)) {
          visited <- c(visited, current)
          component <- c(component, current)

          neighbors <- adj[[current]]
          for (neighbor in neighbors) {
            if (!(neighbor %in% visited)) {
              queue <- c(queue, neighbor)
            }
          }
        }
      }

      if (length(component) > 1) {
        clusters[[length(clusters) + 1]] <- sort(component)
      }
    }
  }

  clusters
}


#' @title Compute Permutation Feature Importance
#'
#' @description
#' Calculates permutation-based feature importance using DALEX.
#'
#' @param explainer A DALEX explainer
#' @param n_perm Number of permutations (default: 10)
#' @param loss_function Loss function ("1-auc" for classification)
#' @param features Optional: subset of features to evaluate
#'
#' @return A model_parts object with importance scores
#'
#' @export
xai_importance <- function(explainer,
                            n_perm = 10,
                            loss_function = "1 - auc",
                            features = NULL) {

  if (!requireNamespace("DALEX", quietly = TRUE)) {
    stop("Package 'DALEX' is required", call. = FALSE)
  }

  # Subset features if specified
  if (!is.null(features)) {
    available <- intersect(features, colnames(explainer$data))
    if (length(available) == 0) {
      stop("No specified features found in explainer", call. = FALSE)
    }
  } else {
    features <- colnames(explainer$data)
  }

  # Compute permutation importance
  mp <- DALEX::model_parts(
    explainer,
    variables = features,
    B = n_perm,
    type = "difference",
    loss_function = DALEX::loss_one_minus_auc
  )

  mp
}


#' @title Compute Partial Dependence Plots
#'
#' @description
#' Computes PDPs for specified features using DALEX.
#'
#' @param explainer A DALEX explainer
#' @param features Features to compute PDPs for
#' @param grid_n Number of grid points (default: 20)
#' @param grid_type "quantile" or "uniform"
#'
#' @return A model_profile object
#'
#' @export
xai_pdp <- function(explainer,
                     features,
                     grid_n = 20,
                     grid_type = "quantile") {

  if (!requireNamespace("DALEX", quietly = TRUE)) {
    stop("Package 'DALEX' is required", call. = FALSE)
  }

  # Validate features exist
  available <- intersect(features, colnames(explainer$data))
  if (length(available) == 0) {
    stop("No specified features found in explainer", call. = FALSE)
  }

  # Compute profiles
  DALEX::model_profile(
    explainer,
    variables = available,
    N = min(300, nrow(explainer$data)),  # Limit samples for speed
    type = "partial"
  )
}


#' @title Compute SHAP Values for Observations
#'
#' @description
#' Computes SHAP values for individual predictions.
#'
#' @param explainer A DALEX explainer
#' @param new_observations Data frame of observations to explain
#' @param B Number of random orderings (default: 25)
#'
#' @return A list of predict_parts objects (one per observation)
#'
#' @export
xai_shap <- function(explainer,
                      new_observations,
                      B = 25) {

  if (!requireNamespace("DALEX", quietly = TRUE)) {
    stop("Package 'DALEX' is required", call. = FALSE)
  }

  # Ensure correct columns
  common_cols <- intersect(names(new_observations), colnames(explainer$data))
  if (length(common_cols) == 0) {
    stop("No common features between observations and explainer", call. = FALSE)
  }

  # Compute SHAP for each observation
  shap_list <- list()
  for (i in seq_len(nrow(new_observations))) {
    obs <- new_observations[i, common_cols, drop = FALSE]

    shap_list[[i]] <- tryCatch({
      DALEX::predict_parts(
        explainer,
        new_observation = obs,
        type = "shap",
        B = B
      )
    }, error = function(e) {
      warning(sprintf("SHAP failed for observation %d: %s", i, e$message))
      NULL
    })
  }

  shap_list
}


#' @title Run Complete XAI Pipeline
#'
#' @description
#' Runs the full explainability pipeline: importance, PDPs, SHAP, correlations.
#'
#' @param learner A trained mlr3 learner
#' @param task The mlr3 task
#' @param top_k Number of top features to analyze in detail (default: 10)
#' @param n_shap_obs Number of observations for SHAP (default: 5)
#' @param cor_threshold Correlation threshold for warnings (default: 0.7)
#' @param features Optional: pre-selected feature subset
#' @param verbose Print progress (default: TRUE)
#'
#' @return A list containing all XAI results
#'
#' @examples
#' \dontrun{
#' results <- xai_pipeline(learner, task, top_k = 10)
#' plot(results$importance)
#' print(results$correlations$high_cor_pairs)
#' }
#'
#' @export
xai_pipeline <- function(learner,
                          task,
                          top_k = 10,
                          n_shap_obs = 5,
                          cor_threshold = 0.7,
                          features = NULL,
                          verbose = TRUE) {

  if (verbose) cat("XAI Pipeline Starting...\n")

  # Get training data
  data <- as.data.frame(task$data())
  target_name <- task$target_names

  # Get learner's feature names
  learner_features <- learner$state$feature_names
  X <- data[, learner_features, drop = FALSE]

  # Determine analysis features (subset if specified)
  analysis_features <- learner_features
  if (!is.null(features)) {
    analysis_features <- intersect(features, learner_features)
    if (length(analysis_features) == 0) {
      analysis_features <- learner_features
    }
  }

  # Step 1: Create explainer
  if (verbose) cat("  [1/5] Creating DALEX explainer...\n")
  explainer <- xai_explainer_mlr3(learner, task)

  # Step 2: Correlation diagnostics
  if (verbose) cat("  [2/5] Computing correlation diagnostics...\n")
  correlations <- xai_correlations(X, threshold = cor_threshold)

  if (correlations$n_warnings > 0 && verbose) {
    cat(sprintf("    WARNING: %d highly correlated feature pairs (|r| > %.2f)\n",
                correlations$n_warnings, cor_threshold))
  }

  # Step 3: Permutation importance
  if (verbose) cat("  [3/5] Computing permutation importance...\n")
  importance <- tryCatch({
    xai_importance(explainer, n_perm = 10, features = analysis_features)
  }, error = function(e) {
    warning("Permutation importance failed: ", e$message)
    NULL
  })

  # Get top features by importance
  if (!is.null(importance)) {
    imp_df <- as.data.frame(importance)
    imp_df <- imp_df[!imp_df$variable %in% c("_full_model_", "_baseline_"), ]
    imp_agg <- aggregate(dropout_loss ~ variable, data = imp_df, FUN = mean)
    imp_agg <- imp_agg[order(-imp_agg$dropout_loss), ]
    top_features <- head(as.character(imp_agg$variable), top_k)
  } else {
    # Fallback to variance-based selection from analysis features
    X_sub <- X[, analysis_features, drop = FALSE]
    vars <- sapply(X_sub, var, na.rm = TRUE)
    top_features <- names(sort(vars, decreasing = TRUE))[1:min(top_k, length(vars))]
  }

  # Step 4: PDPs for top features
  if (verbose) cat("  [4/5] Computing PDPs for top features...\n")
  pdps <- tryCatch({
    xai_pdp(explainer, features = top_features)
  }, error = function(e) {
    warning("PDP computation failed: ", e$message)
    NULL
  })

  # Step 5: SHAP for selected observations
  if (verbose) cat("  [5/5] Computing SHAP values...\n")

  # Select diverse observations for SHAP
  # Get predictions to select interesting cases
  pred <- learner$predict_newdata(X)
  probs <- pred$prob[, task$positive]

  # Select: high confidence positive, high confidence negative, borderline
  obs_indices <- c(
    which.max(probs),  # Most confident positive
    which.min(probs),  # Most confident negative
    order(abs(probs - 0.5))[1:min(3, n_shap_obs - 2)]  # Borderline cases
  )
  obs_indices <- unique(obs_indices[1:min(length(obs_indices), n_shap_obs)])

  shap <- tryCatch({
    xai_shap(explainer, X[obs_indices, , drop = FALSE], B = 15)
  }, error = function(e) {
    warning("SHAP computation failed: ", e$message)
    NULL
  })

  if (verbose) cat("XAI Pipeline Complete.\n")

  # Compile results
  list(
    explainer = explainer,
    correlations = correlations,
    importance = importance,
    top_features = top_features,
    pdps = pdps,
    shap = shap,
    shap_observations = obs_indices,
    config = list(
      top_k = top_k,
      n_shap_obs = n_shap_obs,
      cor_threshold = cor_threshold
    )
  )
}


#' @title Print XAI Summary
#'
#' @description
#' Prints a summary of XAI results with correlation warnings.
#'
#' @param xai_results Output from xai_pipeline
#'
#' @export
print_xai_summary <- function(xai_results) {

  cat("\n========================================\n")
  cat("         XAI SUMMARY REPORT\n")
  cat("========================================\n\n")

  # Correlations
  cat("CORRELATION DIAGNOSTICS\n")
  cat("------------------------\n")
  cor_info <- xai_results$correlations
  cat(sprintf("Threshold: |r| > %.2f\n", cor_info$threshold))
  cat(sprintf("High correlation pairs: %d\n", cor_info$n_warnings))

  if (cor_info$n_warnings > 0) {
    cat("\nWARNING: SHAP values may be unreliable for correlated features!\n")
    cat("Top correlated pairs:\n")
    n_show <- min(5, nrow(cor_info$high_cor_pairs))
    for (i in seq_len(n_show)) {
      cat(sprintf("  - %s <-> %s (r = %.3f)\n",
                  cor_info$high_cor_pairs$feature1[i],
                  cor_info$high_cor_pairs$feature2[i],
                  cor_info$high_cor_pairs$correlation[i]))
    }
  }

  # Feature importance
  cat("\nFEATURE IMPORTANCE\n")
  cat("------------------\n")
  if (!is.null(xai_results$importance)) {
    cat(sprintf("Top %d features by permutation importance:\n",
                length(xai_results$top_features)))
    for (i in seq_along(xai_results$top_features)) {
      feat <- xai_results$top_features[i]
      # Check if correlated
      is_corr <- feat %in% cor_info$high_cor_pairs$feature1 |
                 feat %in% cor_info$high_cor_pairs$feature2
      corr_marker <- if (is_corr) " *" else ""
      cat(sprintf("  %2d. %s%s\n", i, feat, corr_marker))
    }
    if (any(cor_info$high_cor_pairs$feature1 %in% xai_results$top_features |
            cor_info$high_cor_pairs$feature2 %in% xai_results$top_features)) {
      cat("  (* = highly correlated with another feature)\n")
    }
  } else {
    cat("  (Importance computation failed)\n")
  }

  # SHAP summary
  cat("\nSHAP ANALYSIS\n")
  cat("-------------\n")
  if (!is.null(xai_results$shap)) {
    cat(sprintf("SHAP computed for %d observations\n",
                length(xai_results$shap_observations)))
    cat("Observation indices: ", paste(xai_results$shap_observations, collapse = ", "), "\n")
  } else {
    cat("  (SHAP computation failed)\n")
  }

  cat("\n========================================\n")
  invisible(xai_results)
}


#' @title Plot XAI Feature Importance
#'
#' @description
#' Creates a feature importance plot with correlation warnings.
#'
#' @param xai_results Output from xai_pipeline
#' @param top_n Number of top features to show (default: 15)
#' @param show_cor_warning Highlight correlated features (default: TRUE)
#'
#' @return A ggplot object
#'
#' @export
plot_xai_importance <- function(xai_results,
                                 top_n = 15,
                                 show_cor_warning = TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required", call. = FALSE)
  }

  if (is.null(xai_results$importance)) {
    stop("No importance data available", call. = FALSE)
  }

  # Get importance data
  imp_df <- as.data.frame(xai_results$importance)
  imp_df <- imp_df[imp_df$variable != "_full_model_" & imp_df$variable != "_baseline_", ]

  # Aggregate by variable
  imp_agg <- aggregate(dropout_loss ~ variable, data = imp_df, FUN = mean)
  imp_agg <- imp_agg[order(-imp_agg$dropout_loss), ]
  imp_agg <- head(imp_agg, top_n)

  # Mark correlated features
  cor_pairs <- xai_results$correlations$high_cor_pairs
  correlated_features <- unique(c(cor_pairs$feature1, cor_pairs$feature2))
  imp_agg$correlated <- imp_agg$variable %in% correlated_features

  # Create plot
  p <- ggplot2::ggplot(imp_agg,
                       ggplot2::aes(x = stats::reorder(variable, dropout_loss),
                                    y = dropout_loss,
                                    fill = correlated)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Permutation Feature Importance",
      subtitle = if (show_cor_warning && any(imp_agg$correlated))
        "Orange bars indicate highly correlated features (interpret with caution)"
      else NULL,
      x = "Feature",
      y = "Dropout Loss (higher = more important)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  if (show_cor_warning) {
    p <- p + ggplot2::scale_fill_manual(
      values = c("FALSE" = "#3498db", "TRUE" = "#e67e22")
    )
  } else {
    p <- p + ggplot2::scale_fill_manual(values = c("FALSE" = "#3498db", "TRUE" = "#3498db"))
  }

  p
}
