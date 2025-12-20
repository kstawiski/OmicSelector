#' @title Nogueira Stability Index for Feature Selection
#'
#' @description
#' Functions for computing the Nogueira Stability Index, which measures
#' the consistency of feature selection across cross-validation folds.
#'
#' @details
#' The Nogueira Stability Index is preferred over Jaccard/Kuncheva indices
#' because it is invariant to the number of selected features. A high stability
#' (close to 1) indicates that the same features are consistently selected
#' across different training subsets, suggesting robust signal rather than
#' overfitting to specific data splits.
#'
#' Formula:
#' \deqn{S = 1 - \frac{V}{V_0}}
#'
#' Where:
#' - V = (1/p) * sum(pi_j * (1 - pi_j)) [observed variance]
#' - V_0 = (k_bar/p) * (1 - k_bar/p) [expected variance under random selection]
#' - pi_j = selection frequency for feature j
#' - k_bar = average number of selected features
#' - p = total number of candidate features
#'
#' @name stability
NULL


#' @title Compute Nogueira Stability Index
#'
#' @description
#' Computes the Nogueira Stability Index from a list of selected feature sets
#' across multiple cross-validation folds.
#'
#' @param feature_sets A list where each element is a character vector of
#'   selected feature names from one CV fold
#' @param all_features Character vector of all candidate feature names
#'   (the universe of features before selection)
#'
#' @return A list containing:
#'   - nogueira_index: The stability index (0 to 1, higher = more stable)
#'   - selection_frequency: Named vector of per-feature selection frequencies
#'   - avg_subset_size: Average number of features selected per fold
#'   - subset_sizes: Vector of subset sizes per fold
#'   - consensus_features: Features selected in ALL folds
#'   - n_folds: Number of folds analyzed
#'   - n_features: Total number of candidate features
#'   - warning: Any warning messages (e.g., degenerate cases)
#'
#' @examples
#' \dontrun{
#' # Feature sets from 5-fold CV
#' sets <- list(
#'   c("gene1", "gene2", "gene3"),
#'   c("gene1", "gene2", "gene4"),
#'   c("gene1", "gene2", "gene3"),
#'   c("gene1", "gene3", "gene5"),
#'   c("gene1", "gene2", "gene3")
#' )
#' all_genes <- paste0("gene", 1:100)
#'
#' stability <- compute_nogueira_stability(sets, all_genes)
#' print(stability$nogueira_index)
#' }
#'
#' @export
compute_nogueira_stability <- function(feature_sets, all_features) {

  # Input validation
  if (!is.list(feature_sets)) {
    stop("feature_sets must be a list of character vectors", call. = FALSE)
  }

  m <- length(feature_sets)  # Number of folds
  p <- length(all_features)  # Total features

  if (m < 2) {
    stop("At least 2 folds required to compute stability", call. = FALSE)
  }

  if (p == 0) {
    stop("all_features cannot be empty", call. = FALSE)
  }

  # Initialize result
  result <- list(
    nogueira_index = NA_real_,
    selection_frequency = NULL,
    avg_subset_size = NA_real_,
    subset_sizes = integer(m),
    consensus_features = character(0),
    n_folds = m,
    n_features = p,
    warning = NULL
  )

  # Build binary selection matrix Z (m x p)
  # Z[i,j] = 1 if feature j selected in fold i
  Z <- matrix(0, nrow = m, ncol = p)
  colnames(Z) <- all_features

  for (i in seq_len(m)) {
    selected <- feature_sets[[i]]
    if (length(selected) > 0) {
      # Handle features that might not be in all_features
      valid_selected <- intersect(selected, all_features)
      if (length(valid_selected) > 0) {
        Z[i, valid_selected] <- 1
      }
    }
    result$subset_sizes[i] <- sum(Z[i, ])
  }

  # Compute selection frequency for each feature
  # pi_j = (1/m) * sum_i(Z_ij)
  pi_j <- colMeans(Z)
  result$selection_frequency <- pi_j

  # Average subset size
  # k_bar = sum_j(pi_j) = mean(subset_sizes) * p / p = mean(subset_sizes)
  k_bar <- mean(result$subset_sizes)
  result$avg_subset_size <- k_bar

  # Consensus features (selected in ALL folds)
  result$consensus_features <- names(pi_j)[pi_j == 1]

  # Observed variance V
  # V = (1/p) * sum_j(pi_j * (1 - pi_j))
  V <- mean(pi_j * (1 - pi_j))

  # Expected variance under random selection V_0
  # V_0 = (k_bar/p) * (1 - k_bar/p)
  k_ratio <- k_bar / p
  V_0 <- k_ratio * (1 - k_ratio)

  # Handle edge cases
  if (V_0 == 0) {
    # Degenerate case: all features selected or none selected
    if (k_bar == 0) {
      result$warning <- "Degenerate case: no features selected in any fold"
      result$nogueira_index <- 1.0  # Perfectly stable (all agree on empty set)
    } else if (k_bar == p) {
      result$warning <- "Degenerate case: all features selected in every fold"
      result$nogueira_index <- 1.0  # Perfectly stable (all agree on full set)
    } else {
      result$warning <- "Unexpected V_0 = 0 condition"
      result$nogueira_index <- NA_real_
    }
  } else {
    # Standard case: compute Nogueira index
    result$nogueira_index <- 1 - (V / V_0)
  }

  # Add interpretation
  result$interpretation <- interpret_stability(result$nogueira_index)

  class(result) <- c("NogueiraStability", "list")
  return(result)
}


#' @title Interpret Stability Index Value
#'
#' @description
#' Provides a qualitative interpretation of a stability index value.
#'
#' @param stability_index Numeric stability value (0 to 1)
#'
#' @return Character string with interpretation
#' @keywords internal
interpret_stability <- function(stability_index) {
  if (is.na(stability_index)) {
    return("Unable to compute (degenerate case)")
  }

  if (stability_index >= 0.9) {
    return("Excellent: Very stable feature selection")
  } else if (stability_index >= 0.7) {
    return("Good: Reasonably stable feature selection")
  } else if (stability_index >= 0.5) {
    return("Moderate: Some instability in feature selection")
  } else if (stability_index >= 0.3) {
    return("Poor: Unstable feature selection - possible overfitting")
  } else {
    return("Very Poor: Highly unstable - likely random/overfit")
  }
}


#' @title Print method for NogueiraStability
#' @param x A NogueiraStability object
#' @param ... Additional arguments (ignored)
#' @export
print.NogueiraStability <- function(x, ...) {
  cat("=== Nogueira Stability Index ===\n\n")

  cat(sprintf("Stability Index: %.4f\n", x$nogueira_index))
  cat(sprintf("Interpretation: %s\n\n", x$interpretation))

  cat(sprintf("Folds analyzed: %d\n", x$n_folds))
  cat(sprintf("Total features: %d\n", x$n_features))
  cat(sprintf("Avg selected per fold: %.1f\n", x$avg_subset_size))
  cat(sprintf("Subset size range: %d - %d\n",
              min(x$subset_sizes), max(x$subset_sizes)))

  cat(sprintf("\nConsensus features (in all folds): %d\n",
              length(x$consensus_features)))
  if (length(x$consensus_features) > 0 && length(x$consensus_features) <= 10) {
    cat("  ", paste(x$consensus_features, collapse = ", "), "\n")
  } else if (length(x$consensus_features) > 10) {
    cat("  ", paste(x$consensus_features[1:10], collapse = ", "), ", ...\n")
  }

  if (!is.null(x$warning)) {
    cat(sprintf("\nWarning: %s\n", x$warning))
  }

  invisible(x)
}


#' @title Extract Selected Features from Trained GraphLearner
#'
#' @description
#' Extracts the selected feature names from a trained GraphLearner by
#' inspecting the filter PipeOp's state.
#'
#' @param learner A trained GraphLearner object
#' @param filter_id The ID of the filter PipeOp in the graph (default: "filter")
#'
#' @return Character vector of selected feature names, or NULL if extraction fails
#'
#' @export
extract_selected_features <- function(learner, filter_id = "filter") {

  # Check if learner is trained
  if (is.null(learner$model)) {
    warning("Learner is not trained. Cannot extract features.")
    return(NULL)
  }

  # Handle GraphLearner

  if (inherits(learner, "GraphLearner")) {
    graph <- learner$graph

    # Find the filter PipeOp by ID
    if (!filter_id %in% graph$ids()) {
      # Try common filter ID patterns
      filter_ids <- grep("filter|select|flt", graph$ids(), value = TRUE, ignore.case = TRUE)
      if (length(filter_ids) > 0) {
        filter_id <- filter_ids[1]
        message(sprintf("Using filter PipeOp: %s", filter_id))
      } else {
        warning("Could not find filter PipeOp in graph")
        return(NULL)
      }
    }

    # Get the trained PipeOp
    pipeop <- graph$pipeops[[filter_id]]

    # Extract features from state
    if (!is.null(pipeop$state)) {
      # Different PipeOps store features differently
      # Try common patterns
      if (!is.null(pipeop$state$features)) {
        return(pipeop$state$features)
      }
      if (!is.null(pipeop$state$selected_features)) {
        return(pipeop$state$selected_features)
      }
      if (!is.null(pipeop$state$outtasklayout)) {
        # For PipeOpFilter, features are in outtasklayout
        return(pipeop$state$outtasklayout[[1]]$id)
      }
    }

    warning("Could not extract features from filter PipeOp state")
    return(NULL)
  }

  # Handle AutoFSelector results
  if (inherits(learner, "AutoFSelector") || !is.null(learner$fselect_result)) {
    if (!is.null(learner$fselect_result)) {
      return(learner$fselect_result$features)
    }
  }

  warning("Unsupported learner type for feature extraction")
  return(NULL)
}


#' @title Extract Features from All Folds in ResampleResult
#'
#' @description
#' Extracts selected features from each fold of a ResampleResult.
#'
#' @param resample_result An mlr3 ResampleResult object
#' @param filter_id The ID of the filter PipeOp (default: "filter")
#'
#' @return A list of character vectors, one per fold
#'
#' @export
extract_features_from_resample <- function(resample_result, filter_id = "filter") {

  n_iters <- resample_result$iters

  feature_sets <- vector("list", n_iters)

  for (i in seq_len(n_iters)) {
    learner <- resample_result$learners[[i]]
    features <- extract_selected_features(learner, filter_id)
    feature_sets[[i]] <- if (is.null(features)) character(0) else features
  }

  return(feature_sets)
}


#' @title Compute Stability from ResampleResult
#'
#' @description
#' Convenience function to compute Nogueira Stability directly from
#' an mlr3 ResampleResult.
#'
#' @param resample_result An mlr3 ResampleResult object
#' @param all_features Character vector of all candidate features
#' @param filter_id The ID of the filter PipeOp (default: "filter")
#'
#' @return A NogueiraStability object
#'
#' @export
compute_stability_from_resample <- function(resample_result, all_features,
                                             filter_id = "filter") {

  # Extract feature sets from all folds
  feature_sets <- extract_features_from_resample(resample_result, filter_id)

  # Check if extraction succeeded
  if (all(sapply(feature_sets, length) == 0)) {
    warning("No features extracted from any fold. Check filter_id and learner state.")
    return(NULL)
  }

  # Compute stability
  compute_nogueira_stability(feature_sets, all_features)
}
