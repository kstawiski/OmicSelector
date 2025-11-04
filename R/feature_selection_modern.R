#' @title Advanced Feature Selection with Stability Metrics
#' @description Modern feature selection methods with stability assessment
#'
#' This module implements advanced feature selection techniques that go beyond
#' traditional methods by incorporating stability metrics and ensuring
#' reproducibility across different data subsamples.
#'
#' Key features:
#' - Stability selection (Meinshausen & Bühlmann, 2010)
#' - Nogueira stability metrics (Nogueira & Brown, 2016)
#' - Feature clustering for biomarker replaceability
#' - Model-X Knockoffs for FDR control (Phase 2b)
#'
#' @author Konrad Stawiski
#' @name feature_selection_modern
NULL


#' Stability-Aware Feature Selection
#'
#' Performs feature selection with stability assessment using multiple approaches.
#' This function addresses the problem that many feature selection methods
#' produce unstable results that vary significantly across different data samples.
#'
#' @param data A data frame containing features and outcome
#' @param outcome Character string specifying the outcome variable name
#' @param method Character string specifying the selection method:
#'   \itemize{
#'     \item "stability_selection" - Subsampling-based stability selection
#'     \item "boruta_stable" - Boruta with stability assessment
#'     \item "lasso_stable" - LASSO with multiple runs and stability
#'     \item "rfe_stable" - RFE with stability assessment
#'   }
#' @param n_iterations Integer, number of iterations for stability assessment. Default: 100
#' @param selection_threshold Numeric (0-1), minimum selection frequency to consider
#'   a feature stable. Default: 0.6 (selected in 60% of iterations)
#' @param subsample_rate Numeric (0-1), proportion of data to use in each iteration.
#'   Default: 0.8
#' @param max_features Integer, maximum number of features to select. NULL = no limit.
#' @param parallel Logical, whether to use parallel processing. Default: TRUE
#' @param cores Integer, number of cores for parallel processing
#' @param seed Integer for reproducibility
#' @param verbose Logical, whether to print progress messages. Default: TRUE
#'
#' @return A list object of class "OmicSelector_stable_features" containing:
#'   \item{selected_features}{Character vector of selected stable features}
#'   \item{stability_scores}{Named numeric vector of stability scores for each feature}
#'   \item{selection_frequency}{How often each feature was selected across iterations}
#'   \item{nogueira_metrics}{List of Nogueira stability metrics}
#'   \item{selection_path}{Matrix showing which features were selected in each iteration}
#'   \item{method}{Method used}
#'   \item{parameters}{List of parameters used}
#'
#' @examples
#' \dontrun{
#' # Load data
#' data(miR_Asakura)
#'
#' # Stability selection
#' stable_features <- OmicSelector_stable_features(
#'   data = miR_Asakura,
#'   outcome = "Class",
#'   method = "stability_selection",
#'   n_iterations = 100,
#'   selection_threshold = 0.7
#' )
#'
#' # View results
#' print(stable_features)
#' plot(stable_features)
#'
#' # Get most stable features
#' top_features <- head(stable_features$selected_features, 20)
#' }
#'
#' @references
#' Meinshausen, N., & Bühlmann, P. (2010). Stability selection.
#'   Journal of the Royal Statistical Society: Series B, 72(4), 417-473.
#'
#' Nogueira, S., & Brown, G. (2016). Measuring the stability of feature selection.
#'   In Joint European Conference on Machine Learning and Knowledge Discovery
#'   in Databases (pp. 442-457). Springer.
#'
#' @export
OmicSelector_stable_features <- function(
  data,
  outcome,
  method = c("stability_selection", "boruta_stable", "lasso_stable", "rfe_stable"),
  n_iterations = 100,
  selection_threshold = 0.6,
  subsample_rate = 0.8,
  max_features = NULL,
  parallel = TRUE,
  cores = parallel::detectCores() - 1,
  seed = 123,
  verbose = TRUE
) {

  # Validate inputs
  method <- match.arg(method)

  if (!outcome %in% names(data)) {
    stop(paste0("Outcome variable '", outcome, "' not found in data"))
  }

  if (selection_threshold < 0 || selection_threshold > 1) {
    stop("selection_threshold must be between 0 and 1")
  }

  if (subsample_rate <= 0 || subsample_rate >= 1) {
    stop("subsample_rate must be between 0 and 1")
  }

  if (verbose) {
    message("Starting stability-aware feature selection...")
    message(paste0("Method: ", method))
    message(paste0("Iterations: ", n_iterations))
    message(paste0("Selection threshold: ", selection_threshold))
  }

  # Set seed
  set.seed(seed)

  # Get feature names (exclude outcome)
  feature_names <- setdiff(names(data), outcome)
  n_features <- length(feature_names)

  if (verbose) {
    message(paste0("Features to evaluate: ", n_features))
  }

  # Set up parallel processing
  if (parallel && n_iterations > 10) {
    if (requireNamespace("doParallel", quietly = TRUE)) {
      library(doParallel)
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      on.exit(stopCluster(cl), add = TRUE)
      if (verbose) message(paste0("Using ", cores, " cores for parallel processing"))
    } else {
      parallel <- FALSE
      if (verbose) message("doParallel not available, using sequential processing")
    }
  } else {
    parallel <- FALSE
  }

  # Perform stability-based feature selection
  if (method == "stability_selection") {
    result <- .stability_selection_lasso(
      data = data,
      outcome = outcome,
      feature_names = feature_names,
      n_iterations = n_iterations,
      subsample_rate = subsample_rate,
      selection_threshold = selection_threshold,
      max_features = max_features,
      seed = seed,
      verbose = verbose
    )
  } else if (method == "boruta_stable") {
    result <- .stability_selection_boruta(
      data = data,
      outcome = outcome,
      feature_names = feature_names,
      n_iterations = n_iterations,
      subsample_rate = subsample_rate,
      selection_threshold = selection_threshold,
      max_features = max_features,
      seed = seed,
      verbose = verbose
    )
  } else if (method == "lasso_stable") {
    result <- .stability_selection_lasso(
      data = data,
      outcome = outcome,
      feature_names = feature_names,
      n_iterations = n_iterations,
      subsample_rate = subsample_rate,
      selection_threshold = selection_threshold,
      max_features = max_features,
      seed = seed,
      verbose = verbose
    )
  } else if (method == "rfe_stable") {
    result <- .stability_selection_rfe(
      data = data,
      outcome = outcome,
      feature_names = feature_names,
      n_iterations = n_iterations,
      subsample_rate = subsample_rate,
      selection_threshold = selection_threshold,
      max_features = max_features,
      seed = seed,
      verbose = verbose
    )
  }

  # Calculate Nogueira stability metrics
  if (verbose) message("Calculating stability metrics...")
  nogueira_metrics <- .calculate_nogueira_stability(result$selection_path)

  # Create final result object
  final_result <- list(
    selected_features = result$selected_features,
    stability_scores = result$stability_scores,
    selection_frequency = result$selection_frequency,
    nogueira_metrics = nogueira_metrics,
    selection_path = result$selection_path,
    method = method,
    parameters = list(
      n_iterations = n_iterations,
      selection_threshold = selection_threshold,
      subsample_rate = subsample_rate,
      max_features = max_features,
      n_features_evaluated = n_features,
      n_features_selected = length(result$selected_features)
    )
  )

  class(final_result) <- c("OmicSelector_stable_features", "list")

  if (verbose) {
    message(paste0("✓ Selected ", length(result$selected_features), " stable features"))
    message(paste0("✓ Stability score: ", round(nogueira_metrics$stability, 3)))
  }

  return(final_result)
}


#' Internal: Stability Selection with LASSO
#'
#' @keywords internal
#' @noRd
.stability_selection_lasso <- function(
  data, outcome, feature_names, n_iterations,
  subsample_rate, selection_threshold, max_features, seed, verbose
) {

  # Check for glmnet
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("glmnet package required. Install with: install.packages('glmnet')")
  }

  library(glmnet)

  n_features <- length(feature_names)
  n_samples <- nrow(data)
  subsample_size <- floor(subsample_rate * n_samples)

  # Matrix to store selection results
  selection_matrix <- matrix(0, nrow = n_iterations, ncol = n_features)
  colnames(selection_matrix) <- feature_names

  # Determine outcome type
  is_classification <- is.factor(data[[outcome]]) || is.character(data[[outcome]])
  family <- ifelse(is_classification, "binomial", "gaussian")

  # Prepare data
  x_full <- as.matrix(data[, feature_names, drop = FALSE])
  y_full <- if (is_classification) as.factor(data[[outcome]]) else data[[outcome]]

  # Run iterations
  if (verbose) pb <- txtProgressBar(min = 0, max = n_iterations, style = 3)

  for (i in 1:n_iterations) {
    # Subsample
    subsample_idx <- sample(1:n_samples, size = subsample_size, replace = FALSE)
    x_sub <- x_full[subsample_idx, , drop = FALSE]
    y_sub <- y_full[subsample_idx]

    # Fit LASSO with cross-validation
    tryCatch({
      cv_fit <- cv.glmnet(x_sub, y_sub, family = family, alpha = 1, nfolds = 5)

      # Extract selected features at lambda.min
      coef_matrix <- coef(cv_fit, s = "lambda.min")
      selected_idx <- which(coef_matrix != 0)

      # Remove intercept
      selected_idx <- selected_idx[selected_idx > 1] - 1

      if (length(selected_idx) > 0) {
        selection_matrix[i, selected_idx] <- 1
      }
    }, error = function(e) {
      # Skip this iteration if error
      if (verbose) message(paste0("\nIteration ", i, " failed: ", conditionMessage(e)))
    })

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) close(pb)

  # Calculate selection frequencies
  selection_frequency <- colSums(selection_matrix) / n_iterations
  names(selection_frequency) <- feature_names

  # Select stable features
  stable_features <- names(selection_frequency[selection_frequency >= selection_threshold])

  # Apply max_features if specified
  if (!is.null(max_features) && length(stable_features) > max_features) {
    # Sort by selection frequency and keep top max_features
    sorted_features <- sort(selection_frequency[stable_features], decreasing = TRUE)
    stable_features <- names(sorted_features)[1:max_features]
  }

  result <- list(
    selected_features = stable_features,
    stability_scores = selection_frequency,
    selection_frequency = selection_frequency,
    selection_path = selection_matrix
  )

  return(result)
}


#' Internal: Stability Selection with Boruta
#'
#' @keywords internal
#' @noRd
.stability_selection_boruta <- function(
  data, outcome, feature_names, n_iterations,
  subsample_rate, selection_threshold, max_features, seed, verbose
) {

  # Check for Boruta
  if (!requireNamespace("Boruta", quietly = TRUE)) {
    stop("Boruta package required. Install with: install.packages('Boruta')")
  }

  library(Boruta)

  n_features <- length(feature_names)
  n_samples <- nrow(data)
  subsample_size <- floor(subsample_rate * n_samples)

  # Matrix to store selection results
  selection_matrix <- matrix(0, nrow = n_iterations, ncol = n_features)
  colnames(selection_matrix) <- feature_names

  # Prepare formula
  formula_obj <- as.formula(paste0(outcome, " ~ ."))

  # Run iterations
  if (verbose) pb <- txtProgressBar(min = 0, max = n_iterations, style = 3)

  for (i in 1:n_iterations) {
    # Subsample
    subsample_idx <- sample(1:n_samples, size = subsample_size, replace = FALSE)
    data_sub <- data[subsample_idx, ]

    # Run Boruta
    tryCatch({
      boruta_result <- Boruta(formula_obj, data = data_sub, doTrace = 0, maxRuns = 50)

      # Get confirmed features
      decision <- boruta_result$finalDecision
      confirmed_features <- names(decision[decision == "Confirmed"])

      # Mark as selected
      if (length(confirmed_features) > 0) {
        confirmed_idx <- match(confirmed_features, feature_names)
        confirmed_idx <- confirmed_idx[!is.na(confirmed_idx)]
        if (length(confirmed_idx) > 0) {
          selection_matrix[i, confirmed_idx] <- 1
        }
      }
    }, error = function(e) {
      if (verbose) message(paste0("\nIteration ", i, " failed: ", conditionMessage(e)))
    })

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) close(pb)

  # Calculate selection frequencies
  selection_frequency <- colSums(selection_matrix) / n_iterations
  names(selection_frequency) <- feature_names

  # Select stable features
  stable_features <- names(selection_frequency[selection_frequency >= selection_threshold])

  # Apply max_features if specified
  if (!is.null(max_features) && length(stable_features) > max_features) {
    sorted_features <- sort(selection_frequency[stable_features], decreasing = TRUE)
    stable_features <- names(sorted_features)[1:max_features]
  }

  result <- list(
    selected_features = stable_features,
    stability_scores = selection_frequency,
    selection_frequency = selection_frequency,
    selection_path = selection_matrix
  )

  return(result)
}


#' Internal: Stability Selection with RFE
#'
#' @keywords internal
#' @noRd
.stability_selection_rfe <- function(
  data, outcome, feature_names, n_iterations,
  subsample_rate, selection_threshold, max_features, seed, verbose
) {

  # Check for caret
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("caret package required. Install with: install.packages('caret')")
  }

  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("randomForest package required. Install with: install.packages('randomForest')")
  }

  library(caret)
  library(randomForest)

  n_features <- length(feature_names)
  n_samples <- nrow(data)
  subsample_size <- floor(subsample_rate * n_samples)

  # Matrix to store selection results
  selection_matrix <- matrix(0, nrow = n_iterations, ncol = n_features)
  colnames(selection_matrix) <- feature_names

  # RFE control
  rfe_control <- rfeControl(functions = rfFuncs, method = "cv", number = 3, verbose = FALSE)

  # Prepare data
  x_full <- data[, feature_names, drop = FALSE]
  y_full <- data[[outcome]]

  # Run iterations
  if (verbose) pb <- txtProgressBar(min = 0, max = n_iterations, style = 3)

  for (i in 1:n_iterations) {
    # Subsample
    subsample_idx <- sample(1:n_samples, size = subsample_size, replace = FALSE)
    x_sub <- x_full[subsample_idx, , drop = FALSE]
    y_sub <- y_full[subsample_idx]

    # Run RFE
    tryCatch({
      # Use smaller subset of sizes for speed
      sizes <- c(5, 10, 15, 20, min(30, n_features))
      sizes <- unique(sort(sizes[sizes <= n_features]))

      rfe_result <- rfe(x_sub, y_sub, sizes = sizes, rfeControl = rfe_control)

      # Get selected features
      selected_features <- rfe_result$optVariables

      if (length(selected_features) > 0) {
        selected_idx <- match(selected_features, feature_names)
        selected_idx <- selected_idx[!is.na(selected_idx)]
        if (length(selected_idx) > 0) {
          selection_matrix[i, selected_idx] <- 1
        }
      }
    }, error = function(e) {
      if (verbose) message(paste0("\nIteration ", i, " failed: ", conditionMessage(e)))
    })

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) close(pb)

  # Calculate selection frequencies
  selection_frequency <- colSums(selection_matrix) / n_iterations
  names(selection_frequency) <- feature_names

  # Select stable features
  stable_features <- names(selection_frequency[selection_frequency >= selection_threshold])

  # Apply max_features if specified
  if (!is.null(max_features) && length(stable_features) > max_features) {
    sorted_features <- sort(selection_frequency[stable_features], decreasing = TRUE)
    stable_features <- names(sorted_features)[1:max_features]
  }

  result <- list(
    selected_features = stable_features,
    stability_scores = selection_frequency,
    selection_frequency = selection_frequency,
    selection_path = selection_matrix
  )

  return(result)
}


#' Internal: Calculate Nogueira Stability Metrics
#'
#' Calculates various stability metrics as described in Nogueira & Brown (2016)
#'
#' @keywords internal
#' @noRd
.calculate_nogueira_stability <- function(selection_matrix) {

  n_iterations <- nrow(selection_matrix)
  n_features <- ncol(selection_matrix)

  # Average number of features selected per iteration
  k_bar <- mean(rowSums(selection_matrix))

  # Stability metric (Kuncheva index)
  # Measures pairwise similarity between feature sets
  if (n_iterations < 2) {
    stability <- NA
  } else {
    pairwise_similarities <- numeric()

    for (i in 1:(n_iterations - 1)) {
      for (j in (i + 1):n_iterations) {
        set_i <- which(selection_matrix[i, ] == 1)
        set_j <- which(selection_matrix[j, ] == 1)

        n_i <- length(set_i)
        n_j <- length(set_j)
        intersection <- length(intersect(set_i, set_j))

        if (n_i > 0 && n_j > 0) {
          # Kuncheva index
          kuncheva <- (intersection - (n_i * n_j / n_features)) /
                      (min(n_i, n_j) - (n_i * n_j / n_features))
          pairwise_similarities <- c(pairwise_similarities, kuncheva)
        }
      }
    }

    stability <- mean(pairwise_similarities, na.rm = TRUE)
  }

  # Jaccard similarity
  jaccard_sims <- numeric()
  if (n_iterations > 1) {
    for (i in 1:(n_iterations - 1)) {
      for (j in (i + 1):n_iterations) {
        set_i <- which(selection_matrix[i, ] == 1)
        set_j <- which(selection_matrix[j, ] == 1)

        if (length(set_i) > 0 || length(set_j) > 0) {
          intersection <- length(intersect(set_i, set_j))
          union <- length(union(set_i, set_j))
          jaccard_sims <- c(jaccard_sims, intersection / union)
        }
      }
    }
  }

  result <- list(
    stability = stability,
    kuncheva_index = stability,
    mean_jaccard = mean(jaccard_sims, na.rm = TRUE),
    median_jaccard = median(jaccard_sims, na.rm = TRUE),
    avg_n_features = k_bar,
    variance_n_features = var(rowSums(selection_matrix))
  )

  return(result)
}


#' Print Method for OmicSelector_stable_features
#'
#' @param x An OmicSelector_stable_features object
#' @param ... Additional arguments (not used)
#' @export
print.OmicSelector_stable_features <- function(x, ...) {
  cat("OmicSelector Stable Feature Selection\n")
  cat("======================================\n\n")

  cat("Method:", x$method, "\n")
  cat("Iterations:", x$parameters$n_iterations, "\n")
  cat("Selection threshold:", x$parameters$selection_threshold, "\n\n")

  cat("Results:\n")
  cat("  Features evaluated:", x$parameters$n_features_evaluated, "\n")
  cat("  Features selected:", x$parameters$n_features_selected, "\n")
  cat("  Selection rate:", round(x$parameters$n_features_selected / x$parameters$n_features_evaluated * 100, 1), "%\n\n")

  cat("Stability Metrics:\n")
  cat("  Kuncheva index:", round(x$nogueira_metrics$kuncheva_index, 3), "\n")
  cat("  Mean Jaccard similarity:", round(x$nogueira_metrics$mean_jaccard, 3), "\n")
  cat("  Avg features per iteration:", round(x$nogueira_metrics$avg_n_features, 1), "\n\n")

  if (length(x$selected_features) > 0) {
    cat("Top 10 stable features:\n")
    top_10 <- head(sort(x$stability_scores[x$selected_features], decreasing = TRUE), 10)
    for (i in seq_along(top_10)) {
      cat(sprintf("  %2d. %-30s (%.3f)\n", i, names(top_10)[i], top_10[i]))
    }
  } else {
    cat("No features met the stability threshold.\n")
  }

  cat("\nUse summary() for detailed results.\n")
}


#' Summary Method for OmicSelector_stable_features
#'
#' @param object An OmicSelector_stable_features object
#' @param ... Additional arguments (not used)
#' @export
summary.OmicSelector_stable_features <- function(object, ...) {
  print(object)

  cat("\n")
  cat("All Selected Features (", length(object$selected_features), "):\n", sep = "")
  cat("=====================================\n")

  if (length(object$selected_features) > 0) {
    selected_scores <- sort(object$stability_scores[object$selected_features], decreasing = TRUE)
    for (i in seq_along(selected_scores)) {
      cat(sprintf("%3d. %-30s %.3f\n", i, names(selected_scores)[i], selected_scores[i]))
    }
  }
}


#' Plot Method for OmicSelector_stable_features
#'
#' @param x An OmicSelector_stable_features object
#' @param ... Additional arguments (not used)
#' @export
plot.OmicSelector_stable_features <- function(x, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting")
  }

  library(ggplot2)

  # Plot 1: Stability scores
  top_50 <- head(sort(x$stability_scores, decreasing = TRUE), 50)

  p1 <- data.frame(
    Feature = factor(names(top_50), levels = rev(names(top_50))),
    Stability = as.numeric(top_50),
    Selected = names(top_50) %in% x$selected_features
  ) %>%
    ggplot(aes(x = Feature, y = Stability, fill = Selected)) +
    geom_col() +
    geom_hline(yintercept = x$parameters$selection_threshold,
               linetype = "dashed", color = "red") +
    coord_flip() +
    scale_fill_manual(values = c("gray70", "steelblue")) +
    labs(title = "Top 50 Features by Stability Score",
         subtitle = paste0("Method: ", x$method, ", Threshold: ", x$parameters$selection_threshold),
         x = "Feature", y = "Selection Frequency") +
    theme_minimal() +
    theme(legend.position = "bottom")

  print(p1)
}
