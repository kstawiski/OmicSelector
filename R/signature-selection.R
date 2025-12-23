#' @title Signature Selection: Multi-Objective Best Biomarker Selection
#'
#' @description
#' Functions for selecting the best biomarker signature from nested CV results.
#' Implements multiple selection strategies: constrained 1SE rule, weighted scoring,
#' and Pareto optimality for transparent multi-objective trade-offs.
#'
#' @details
#' ## Clinical Context
#'
#' In biomarker discovery, selecting the "best" signature requires balancing:
#' - **Performance** (AUC, accuracy): Predictive power
#' - **Stability**: Reproducibility across resampling (Nogueira Index)
#' - **Parsimony**: Fewer features = cheaper, more practical clinical panels
#'
#' This module provides three selection modes:
#' 1. `constrained_1se`: Conservative selection within 1 standard error of best
#' 2. `weighted`: Explicit weighted combination of objectives
#' 3. `pareto`: Return non-dominated solutions for human inspection
#'
#' @name signature-selection
NULL


#' @title Select Best Biomarker Signature from Nested CV Results
#'
#' @description
#' Multi-objective signature selection with three modes: constrained 1SE rule,
#' weighted scoring, or Pareto frontier. Designed for clinical biomarker panels
#' where parsimony and stability are as important as raw performance.
#'
#' @param nested_result A NestedCVResult object from BenchmarkService$run()
#' @param mode Selection mode: "constrained_1se" (default), "weighted", or "pareto"
#' @param metric Character, performance metric column name (default: "classif.auc")
#' @param metric_higher_better Logical, TRUE if higher metric is better (default: TRUE)
#' @param auc_min Numeric, minimum acceptable AUC (hard constraint)
#' @param stability_min Numeric, minimum acceptable stability index (hard constraint)
#' @param k_max Integer, maximum acceptable feature count (hard constraint for clinical panels)
#' @param weights Named numeric vector for weighted mode: c(performance=, stability=, parsimony=)
#' @param na_stability How to handle NA stability: "exclude_if_constrained" (default),
#'   "worst" (treat as worst), or "neutral" (treat as middle)
#'
#' @return A data.table with candidate signatures ranked by the selection criteria.
#'   Contains columns: learner_id, mean_metric, se_metric, stability, mean_k,
#'   score (for weighted), pareto (for pareto mode), selected (best candidate flag)
#'
#' @examples
#' \dontrun{
#' # Run nested CV
#' result <- benchmark_service$run()
#'
#' # Select best with 1SE rule (conservative)
#' best <- select_best_signature(result, mode = "constrained_1se")
#'
#' # Select with explicit weights
#' best <- select_best_signature(result, mode = "weighted",
#'   weights = c(performance = 0.5, stability = 0.3, parsimony = 0.2))
#'
#' # Get Pareto frontier for manual inspection
#' pareto_set <- select_best_signature(result, mode = "pareto")
#' }
#'
#' @export
select_best_signature <- function(
    nested_result,
    mode = c("constrained_1se", "weighted", "pareto"),
    metric = "classif.auc",
    metric_higher_better = TRUE,
    auc_min = NULL,
    stability_min = NULL,
    k_max = NULL,
    weights = c(performance = 0.5, stability = 0.3, parsimony = 0.2),
    na_stability = c("exclude_if_constrained", "worst", "neutral")
) {

  mode <- match.arg(mode)
  na_stability <- match.arg(na_stability)

  # Validate input
  if (!inherits(nested_result, "NestedCVResult")) {
    stop("nested_result must be a NestedCVResult object from BenchmarkService$run()")
  }

  if (is.null(nested_result$benchmark_result)) {
    stop("nested_result must contain benchmark_result")
  }

  # Extract candidate-level summary with per-fold statistics
  candidates <- .extract_candidate_summary(
    nested_result = nested_result,
    metric = metric,
    metric_higher_better = metric_higher_better
  )

  if (nrow(candidates) == 0) {
    warning("No candidates found in nested_result")
    return(candidates)
  }

  # Apply hard constraints
  candidates <- .apply_hard_constraints(
    candidates = candidates,
    auc_min = auc_min,
    stability_min = stability_min,
    k_max = k_max,
    na_stability = na_stability
  )

  if (nrow(candidates) == 0) {
    warning("No candidates remaining after applying constraints")
    return(candidates)
  }

  # Apply mode-specific selection
  result <- switch(mode,
    "constrained_1se" = .select_1se(candidates, na_stability),
    "weighted" = .select_weighted(candidates, weights, na_stability),
    "pareto" = .select_pareto(candidates)
  )

  # Add selection metadata
  attr(result, "mode") <- mode
  attr(result, "metric") <- metric
  attr(result, "constraints") <- list(auc_min = auc_min, stability_min = stability_min, k_max = k_max)

  return(result)
}


#' @title Extract Candidate Summary from Nested CV Results
#'
#' @description
#' Internal function to extract per-candidate statistics across outer folds.
#'
#' @param nested_result NestedCVResult object
#' @param metric Performance metric column name
#' @param metric_higher_better Whether higher metric is better
#'
#' @return data.table with candidate-level summary
#' @keywords internal
.extract_candidate_summary <- function(nested_result, metric, metric_higher_better) {

  bmr <- nested_result$benchmark_result
  stab <- nested_result$stability

  # Get per-fold scores from benchmark result
  # bmr$score() returns per-iteration (fold) results
  # First try without measures, then with the requested measure if needed
  score_dt <- data.table::as.data.table(bmr$score())

  # Check if metric exists, if not try to add it via measures

  if (!metric %in% names(score_dt)) {
    tryCatch({
      # Try to create and pass the measure
      if (requireNamespace("mlr3", quietly = TRUE)) {
        msr_obj <- mlr3::msr(metric)
        score_dt <- data.table::as.data.table(bmr$score(list(msr_obj)))
      }
    }, error = function(e) NULL)
  }

  # Final check
  if (!metric %in% names(score_dt)) {
    available <- names(score_dt)[grepl("^classif\\.|^regr\\.", names(score_dt))]
    stop(sprintf(
      "Metric '%s' not found. Available metrics: %s\nHint: Ensure the metric was computed during benchmark run.",
      metric, paste(available, collapse = ", ")
    ))
  }

  # Aggregate across outer folds
  # Keys: learner_id identifies the candidate
  candidates <- score_dt[, .(
    mean_metric = mean(get(metric), na.rm = TRUE),
    sd_metric = stats::sd(get(metric), na.rm = TRUE),
    n_folds = sum(!is.na(get(metric)))
  ), by = .(learner_id)]

  # Compute SE

  candidates[, se_metric := fifelse(n_folds > 1, sd_metric / sqrt(n_folds), NA_real_)]

  # Store original metric value for constraint checking
  # Constraints should be specified in original metric units (e.g., auc_min=0.7 for AUC)
  candidates[, original_metric := mean_metric]

  # If higher is not better (e.g., brier score), flip sign for internal ranking
  # This allows all selection logic to use "higher is better" convention
  if (!metric_higher_better) {
    candidates[, mean_metric := -mean_metric]
  }

  # Store metric direction for constraint interpretation
  attr(candidates, "metric_higher_better") <- metric_higher_better

  # Attach stability information
  candidates[, stability := NA_real_]

  if (!is.null(stab$nogueira_index)) {
    ni <- stab$nogueira_index

    if (is.numeric(ni) && length(ni) == 1 && !is.na(ni)) {
      # Single global stability value
      candidates[, stability := ni]
    } else if (is.data.frame(ni) || data.table::is.data.table(ni)) {
      # Per-learner stability
      ni_dt <- data.table::as.data.table(ni)
      if ("learner_id" %in% names(ni_dt) && "nogueira_index" %in% names(ni_dt)) {
        candidates <- merge(candidates, ni_dt[, .(learner_id, stability = nogueira_index)],
                            by = "learner_id", all.x = TRUE)
      }
    }
  }

  # Attach feature count (parsimony) information
  candidates[, mean_k := NA_real_]

  if (!is.null(stab$selected_features_per_fold)) {
    sf <- stab$selected_features_per_fold

    if (length(sf) > 0) {
      sf_summary <- .summarize_feature_counts(sf)

      if (!is.null(sf_summary) && "learner_id" %in% names(sf_summary)) {
        candidates <- merge(candidates, sf_summary[, .(learner_id, mean_k)],
                            by = "learner_id", all.x = TRUE)
      }
    }
  }

  # Try to extract feature counts from graph learner states if not available
  if (all(is.na(candidates$mean_k))) {
    candidates <- .try_extract_feature_counts(candidates, bmr)

    # If we extracted features, compute stability from them
    features_per_fold <- attr(candidates, "features_per_fold")
    if (!is.null(features_per_fold) && length(features_per_fold) > 0) {
      # Get total candidate feature count from original task (NOT union of selected)
      # This is critical for correct Nogueira stability computation
      total_features <- tryCatch({
        # mlr3 benchmark result stores tasks
        if (!is.null(bmr$tasks) && length(bmr$tasks$task) > 0) {
          bmr$tasks$task[[1]]$n_features
        } else {
          NULL
        }
      }, error = function(e) NULL)

      for (lid in names(features_per_fold)) {
        fold_features <- features_per_fold[[lid]]
        if (length(fold_features) >= 2) {
          # Use total candidate features for p, not union of selected features
          # If we couldn't get from task, fall back to union (suboptimal but better than nothing)
          p <- if (!is.null(total_features) && total_features > 0) {
            total_features
          } else {
            # Fallback: use union of selected (with warning on first occurrence)
            length(unique(unlist(fold_features)))
          }

          # Compute Nogueira Stability Index with correct p
          stability_val <- .compute_nogueira_index(fold_features, p)
          candidates[learner_id == lid, stability := stability_val]
        }
      }
    }
  }

  candidates[]
}


#' @title Compute Nogueira Stability Index
#'
#' @description
#' Computes the Nogueira Stability Index for feature selection consistency.
#' The index ranges from 0 (completely unstable) to 1 (perfectly stable).
#'
#' @details
#' The Nogueira Stability Index is defined as:
#' SI = 1 - (observed_variance / max_variance)
#'
#' Where observed variance is the average pairwise disagreement between
#' feature sets, and max variance is what we'd expect from random selection.
#'
#' @param fold_features List of character vectors, each containing features
#'   selected in one fold
#' @param p Total number of features in the dataset
#'
#' @return Numeric stability index between 0 and 1
#' @keywords internal
#'
#' @references
#' Nogueira, S., Sechidis, K., & Brown, G. (2018). On the Stability of
#' Feature Selection Algorithms. Journal of Machine Learning Research, 18(174), 1-54.
.compute_nogueira_index <- function(fold_features, p) {
  n_folds <- length(fold_features)

  if (n_folds < 2 || p == 0) {
    return(NA_real_)
  }

  # Create binary selection matrix: rows = folds, cols = features
  all_features <- unique(unlist(fold_features))
  if (length(all_features) == 0) {
    return(NA_real_)
  }

  # Build selection matrix
  Z <- matrix(0, nrow = n_folds, ncol = length(all_features))
  colnames(Z) <- all_features

  for (i in seq_len(n_folds)) {
    if (length(fold_features[[i]]) > 0) {
      Z[i, fold_features[[i]]] <- 1
    }
  }

  # Compute feature-wise selection frequencies
  pf <- colMeans(Z)  # Proportion of times each feature was selected

  # Mean number of features selected per fold
  k_bar <- mean(rowSums(Z))

  if (k_bar == 0 || k_bar == p) {
    # Edge case: all or none selected
    return(1.0)
  }

  # Nogueira's formula:
  # SI = 1 - (1/k_bar) * sum(pf * (1 - pf)) / (1 - k_bar/p)
  #
  # Simplified: using variance formula
  # Observed variance = (1/p) * sum(pf * (1-pf)) where p = TOTAL candidate features
  # Maximum variance under random selection = k_bar/p * (1 - k_bar/p)
  #
  # NOTE: Features never selected have pf=0, contributing 0 to the sum,
  # so we only need to sum over the union of selected features
  # but MUST divide by total candidates p, not the union size

  observed_var <- sum(pf * (1 - pf)) / p  # Use p (total candidates), not length(all_features)
  max_var <- (k_bar / p) * (1 - k_bar / p)

  if (max_var == 0) {
    return(1.0)
  }

  stability <- 1 - observed_var / max_var

  # Clamp to [0, 1] (numerical precision issues)
  stability <- max(0, min(1, stability))

  return(stability)
}


#' @title Summarize Feature Counts from Selected Features Per Fold
#' @keywords internal
.summarize_feature_counts <- function(sf) {
  if (is.data.frame(sf) || data.table::is.data.table(sf)) {
    sf_dt <- data.table::as.data.table(sf)

    if ("selected_features" %in% names(sf_dt)) {
      # List column of feature vectors
      sf_dt[, k := vapply(selected_features, length, integer(1))]
    } else if ("features" %in% names(sf_dt)) {
      sf_dt[, k := vapply(features, length, integer(1))]
    } else if ("n_features" %in% names(sf_dt)) {
      sf_dt[, k := n_features]
    } else {
      return(NULL)
    }

    if ("learner_id" %in% names(sf_dt)) {
      return(sf_dt[, .(mean_k = mean(k, na.rm = TRUE)), by = .(learner_id)])
    }
  }

  return(NULL)
}


#' @title Try to Extract Feature Counts from Benchmark Result
#' @keywords internal
.try_extract_feature_counts <- function(candidates, bmr) {
  # Attempt to extract feature counts from trained learners
  # This works for AutoFSelector and GraphLearner with feature selection

  # Also collect features per fold for stability computation
  features_per_fold <- list()

  tryCatch({
    for (i in seq_len(bmr$n_resample_results)) {
      rr <- bmr$resample_result(i)
      current_learner_id <- rr$learner$id

      fold_features <- list()

      # Try to get selected features from each iteration
      k_vals <- vapply(seq_len(rr$iters), function(iter) {
        lrn <- tryCatch(rr$learners[[iter]], error = function(e) NULL)
        if (is.null(lrn)) return(NA_integer_)

        features <- .extract_features_from_learner(lrn)

        if (!is.null(features) && length(features) > 0) {
          fold_features[[iter]] <<- features
          return(length(features))
        }

        return(NA_integer_)
      }, integer(1))

      mean_k <- mean(k_vals, na.rm = TRUE)
      if (!is.na(mean_k)) {
        candidates[learner_id == current_learner_id, mean_k := mean_k]
      }

      if (length(fold_features) > 0) {
        features_per_fold[[current_learner_id]] <- fold_features
      }
    }
  }, error = function(e) {
    # Silently fail - feature counts are optional
  })

  # Store features_per_fold as attribute for stability computation
  attr(candidates, "features_per_fold") <- features_per_fold

  candidates
}


#' @title Extract Selected Features from a Trained Learner
#'
#' @description
#' Internal function to extract selected feature names from various learner types.
#' Handles AutoFSelector, GraphLearner with filter PipeOps, and other patterns.
#'
#' @param lrn A trained learner object
#' @return Character vector of selected feature names, or NULL if not extractable
#' @keywords internal
.extract_features_from_learner <- function(lrn) {
  features <- NULL


  # Method 1: AutoFSelector
  if (inherits(lrn, "AutoFSelector") && !is.null(lrn$fselect_result)) {
    features <- lrn$fselect_result$features
    if (!is.null(features)) return(features)
  }

  # Method 2: GraphLearner with filter/selector PipeOps
  if (inherits(lrn, "GraphLearner") && !is.null(lrn$state$model)) {
    state_model <- lrn$state$model

    for (nm in names(state_model)) {
      po_state <- state_model[[nm]]

      if (!is.list(po_state)) next

      # Check various field names used by different PipeOps
      # PipeOpFilter stores in $features
      if ("features" %in% names(po_state)) {
        features <- po_state$features
        if (!is.null(features) && length(features) > 0) return(features)
      }

      # Some PipeOps use $selected_features
      if ("selected_features" %in% names(po_state)) {
        features <- po_state$selected_features
        if (!is.null(features) && length(features) > 0) return(features)
      }

      # PipeOpSelect stores in $selection
      if ("selection" %in% names(po_state)) {
        sel <- po_state$selection
        if (is.character(sel)) {
          features <- sel
          if (length(features) > 0) return(features)
        }
      }
    }
  }

  # Method 3: Check model$features directly (some wrappers)
  if (!is.null(lrn$model$features)) {
    features <- lrn$model$features
    if (!is.null(features) && length(features) > 0) return(features)
  }

  # Method 4: For learners with feature_names in state
  if (!is.null(lrn$state$feature_names)) {
    features <- lrn$state$feature_names
    if (!is.null(features) && length(features) > 0) return(features)
  }

  return(NULL)
}


#' @title Apply Hard Constraints to Candidates
#' @keywords internal
.apply_hard_constraints <- function(candidates, auc_min, stability_min, k_max, na_stability) {

  dt <- data.table::copy(candidates)

  # AUC constraint - use original_metric (not sign-flipped) for constraint checking

  # The auc_min parameter should be specified in original metric units
  # (e.g., auc_min=0.7 for AUC, where higher is better)
  if (!is.null(auc_min)) {
    # Use original_metric if available (preserves original scale before sign-flip)
    if ("original_metric" %in% names(dt)) {
      dt <- dt[original_metric >= auc_min]
    } else {
      # Fallback for backwards compatibility
      dt <- dt[mean_metric >= auc_min]
    }
  }

  # Stability constraint
  if (!is.null(stability_min)) {
    if (na_stability == "exclude_if_constrained") {
      dt <- dt[!is.na(stability) & stability >= stability_min]
    } else {
      # Even with "worst" or "neutral", we must exclude NA when there's a hard constraint
      dt <- dt[!is.na(stability) & stability >= stability_min]
    }
  }

  # Feature count constraint (critical for clinical panels)
  if (!is.null(k_max)) {
    # If mean_k is NA, we cannot verify the constraint
    dt <- dt[!is.na(mean_k) & mean_k <= k_max]
  }

  dt
}


#' @title 1SE Selection Rule
#' @keywords internal
.select_1se <- function(candidates, na_stability) {

  dt <- data.table::copy(candidates)

  # Find best candidate by mean metric
  best_idx <- which.max(dt$mean_metric)
  best <- dt[best_idx]

  # Calculate 1SE threshold
  threshold <- best$mean_metric - best$se_metric

  # Handle NA SE (e.g., single fold)
  if (is.na(threshold)) {
    threshold <- best$mean_metric
  }

  # Keep candidates within 1SE of best
  pool <- dt[mean_metric >= threshold]

  # Tie-break: stability (desc) -> parsimony (asc) -> metric (desc)
  # Handle NA values for tie-breaking
  pool[, stability_tb := stability]

  if (na_stability == "worst") {
    pool[is.na(stability_tb), stability_tb := -Inf]
  } else if (na_stability == "neutral") {
    med_stab <- stats::median(pool$stability, na.rm = TRUE)
    if (is.na(med_stab)) med_stab <- 0
    pool[is.na(stability_tb), stability_tb := med_stab]
  }

  # Handle NA mean_k (treat as worst = largest)
  pool[, k_tb := mean_k]
  max_k <- max(pool$mean_k, na.rm = TRUE)
  if (is.finite(max_k)) {
    pool[is.na(k_tb), k_tb := max_k + 1]
  } else {
    pool[is.na(k_tb), k_tb := Inf]
  }

  # Order: best stability, smallest k, highest metric
  pool <- pool[order(-stability_tb, k_tb, -mean_metric)]

  # Mark best selection
  pool[, selected := FALSE]
  pool[1L, selected := TRUE]

  # Remove temporary columns
  pool[, c("stability_tb", "k_tb") := NULL]

  pool[]
}


#' @title Weighted Scoring Selection
#' @keywords internal
.select_weighted <- function(candidates, weights, na_stability) {

  dt <- data.table::copy(candidates)

  # Validate and normalize weights
  w <- weights
  required_names <- c("performance", "stability", "parsimony")

  if (!all(required_names %in% names(w))) {
    stop("weights must have names: ", paste(required_names, collapse = ", "))
  }

  w <- w[required_names]
  w <- w / sum(w)

  # Min-max normalization helper
  minmax01 <- function(x, higher_better = TRUE) {
    x_clean <- x[is.finite(x)]
    if (length(x_clean) == 0) return(rep(0.5, length(x)))

    rng <- range(x_clean, na.rm = TRUE)
    if (rng[1] == rng[2]) {
      return(rep(1, length(x)))
    }

    normalized <- (x - rng[1]) / (rng[2] - rng[1])
    if (!higher_better) normalized <- 1 - normalized

    normalized
  }

  # Normalize performance (higher is better)
  dt[, perf_norm := minmax01(mean_metric, higher_better = TRUE)]

  # Normalize stability (higher is better)
  stab_vec <- dt$stability
  if (all(is.na(stab_vec))) {
    dt[, stab_norm := 0.5]  # No info, neutral contribution
  } else {
    if (na_stability == "worst") {
      min_stab <- min(stab_vec, na.rm = TRUE)
      stab_vec[is.na(stab_vec)] <- min_stab - 0.1
    } else if (na_stability == "neutral") {
      med_stab <- stats::median(stab_vec, na.rm = TRUE)
      stab_vec[is.na(stab_vec)] <- med_stab
    }
    dt[, stab_norm := minmax01(stab_vec, higher_better = TRUE)]
    dt[is.na(stab_norm), stab_norm := 0]
  }

  # Normalize parsimony (smaller k is better)
  k_vec <- dt$mean_k
  if (all(is.na(k_vec))) {
    dt[, pars_norm := 0.5]  # No info, neutral contribution
  } else {
    # NA k treated as worst (largest)
    max_k <- max(k_vec, na.rm = TRUE)
    k_vec[is.na(k_vec)] <- max_k + 1

    dt[, pars_norm := minmax01(k_vec, higher_better = FALSE)]
    dt[is.na(pars_norm), pars_norm := 0]
  }

  # Compute weighted score
  dt[, score := w["performance"] * perf_norm +
                w["stability"] * stab_norm +
                w["parsimony"] * pars_norm]

  # Rank by score, then by mean_metric as tiebreaker
  dt <- dt[order(-score, -mean_metric)]

  # Mark best selection
  dt[, selected := FALSE]
  dt[1L, selected := TRUE]

  dt[]
}


#' @title Pareto Frontier Selection
#' @keywords internal
.select_pareto <- function(candidates) {

  dt <- data.table::copy(candidates)

  # Exclude candidates with missing objectives
  # For Pareto, we need all three objectives
  dt_complete <- dt[!is.na(mean_metric) & !is.na(stability) & !is.na(mean_k)]

  if (nrow(dt_complete) == 0) {
    # Fall back to candidates with at least metric available
    dt_complete <- dt[!is.na(mean_metric)]

    if (nrow(dt_complete) > 0) {
      warning("Stability and/or feature counts unavailable. ",
              "Pareto analysis limited to available objectives.")

      # If only metric available, all are on frontier (single objective)
      dt_complete[, pareto := TRUE]
      dt_complete[, selected := FALSE]

      # Select best by metric
      best_idx <- which.max(dt_complete$mean_metric)
      dt_complete[best_idx, selected := TRUE]

      return(dt_complete)
    }

    return(dt_complete)
  }

  # Compute Pareto frontier
  # Objectives: maximize mean_metric, maximize stability, minimize mean_k
  pareto_idx <- .pareto_frontier_idx(
    dt_complete,
    cols = c("mean_metric", "stability", "mean_k"),
    maximize = c(TRUE, TRUE, FALSE)
  )

  dt_complete[, pareto := FALSE]
  dt_complete[pareto_idx, pareto := TRUE]

  # Compute distance to Utopia for recommended selection
  # Utopia: max metric, max stability, min k (all normalized)
  dt_pareto <- dt_complete[pareto == TRUE]

  if (nrow(dt_pareto) > 0) {
    # Normalize within Pareto set
    dt_pareto[, metric_norm := .minmax_normalize(mean_metric)]
    dt_pareto[, stab_norm := .minmax_normalize(stability)]
    dt_pareto[, k_norm := .minmax_normalize(mean_k)]

    # Distance to Utopia (1, 1, 0)
    dt_pareto[, dist_utopia := sqrt(
      (1 - metric_norm)^2 +
      (1 - stab_norm)^2 +
      (k_norm - 0)^2
    )]

    # Mark recommended (closest to Utopia)
    dt_pareto[, selected := FALSE]
    best_idx <- which.min(dt_pareto$dist_utopia)
    dt_pareto[best_idx, selected := TRUE]

    # Merge back
    dt_complete <- dt_pareto
  } else {
    dt_complete[, selected := FALSE]
  }

  dt_complete[]
}


#' @title Compute Pareto Frontier Indices
#'
#' @description
#' Simple O(n^2) dominance check for Pareto optimality.
#' Dependency-free implementation suitable for typical candidate set sizes.
#'
#' @param dt data.table with objective columns
#' @param cols Character vector of objective column names
#' @param maximize Logical vector (same length as cols): TRUE if higher is better
#'
#' @return Integer vector of row indices that are Pareto-optimal
#' @keywords internal
.pareto_frontier_idx <- function(dt, cols, maximize) {

  stopifnot(length(cols) == length(maximize))

  X <- as.matrix(dt[, ..cols])

  # Transform to "all maximize" by flipping sign for minimization objectives
  for (j in seq_along(cols)) {
    if (!maximize[j]) {
      X[, j] <- -X[, j]
    }
  }

  n <- nrow(X)
  is_dominated <- rep(FALSE, n)

  # Check dominance: A dominates B if all(A >= B) and any(A > B)
  for (i in seq_len(n)) {
    if (is_dominated[i]) next

    for (k in seq_len(n)) {
      if (i == k || is_dominated[k]) next

      # Check if i dominates k
      ge_all <- all(X[i, ] >= X[k, ])
      gt_any <- any(X[i, ] > X[k, ])

      if (ge_all && gt_any) {
        is_dominated[k] <- TRUE
      }
    }
  }

  which(!is_dominated)
}


#' @title Min-Max Normalize Vector to [0, 1]
#' @keywords internal
.minmax_normalize <- function(x) {
  x_clean <- x[is.finite(x)]
  if (length(x_clean) == 0) return(rep(0.5, length(x)))

  rng <- range(x_clean, na.rm = TRUE)
  if (rng[1] == rng[2]) {
    return(rep(1, length(x)))
  }

  (x - rng[1]) / (rng[2] - rng[1])
}


#' @title Get Consensus Features from Best Signature
#'
#' @description
#' Extract the consensus feature set from the selected best signature(s).
#' Aggregates features selected across outer folds.
#'
#' @param nested_result A NestedCVResult object
#' @param best_signature Result from select_best_signature() or learner_id string
#' @param min_frequency Minimum proportion of folds where feature must be selected (default: 0.5)
#'
#' @return A data.table with columns: feature, frequency, selected (TRUE if >= min_frequency)
#'
#' @examples
#' \dontrun{
#' best <- select_best_signature(result)
#' consensus <- get_consensus_features(result, best, min_frequency = 0.8)
#' final_features <- consensus[selected == TRUE]$feature
#' }
#'
#' @export
get_consensus_features <- function(nested_result, best_signature, min_frequency = 0.5) {

  checkmate::assert_number(min_frequency, lower = 0, upper = 1)

  # Extract learner_id
  if (is.character(best_signature)) {
    learner_id <- best_signature
  } else if (is.data.frame(best_signature) || data.table::is.data.table(best_signature)) {
    # Get selected row or first row
    best_dt <- data.table::as.data.table(best_signature)
    if ("selected" %in% names(best_dt)) {
      learner_id <- best_dt[selected == TRUE]$learner_id[1]
    } else {
      learner_id <- best_dt$learner_id[1]
    }
  } else {
    stop("best_signature must be a learner_id string or data.table from select_best_signature()")
  }

  if (is.null(learner_id) || is.na(learner_id)) {
    stop("Could not extract learner_id from best_signature")
  }

  # Try to extract features from stability info
  stab <- nested_result$stability
  sf <- stab$selected_features_per_fold

  if (!is.null(sf) && length(sf) > 0) {
    return(.extract_consensus_from_sf(sf, learner_id, min_frequency))
  }

  # Try to extract from benchmark result
  bmr <- nested_result$benchmark_result

  return(.extract_consensus_from_bmr(bmr, learner_id, min_frequency))
}


#' @title Extract Consensus Features from Selected Features Per Fold
#' @keywords internal
.extract_consensus_from_sf <- function(sf, learner_id, min_frequency) {

  if (data.table::is.data.table(sf) || is.data.frame(sf)) {
    sf_dt <- data.table::as.data.table(sf)

    if (!"learner_id" %in% names(sf_dt)) {
      stop("selected_features_per_fold must have 'learner_id' column")
    }

    # Use a local variable to avoid data.table scoping issues
    # (learner_id parameter would be shadowed by column name)
    target_learner <- learner_id
    sf_learner <- sf_dt[learner_id == target_learner]

    if (nrow(sf_learner) == 0) {
      warning("No feature information found for learner: ", learner_id)
      return(data.table::data.table(feature = character(), frequency = numeric(), selected = logical()))
    }

    # Get feature column
    feature_col <- intersect(c("selected_features", "features"), names(sf_learner))[1]

    if (is.na(feature_col)) {
      stop("selected_features_per_fold must have 'selected_features' or 'features' column")
    }

    # Count feature frequency across folds
    all_features <- unlist(sf_learner[[feature_col]])
    n_folds <- nrow(sf_learner)

    freq_table <- table(all_features) / n_folds
    result <- data.table::data.table(
      feature = names(freq_table),
      frequency = as.numeric(freq_table)
    )

    result[, selected := frequency >= min_frequency]
    result <- result[order(-frequency)]

    return(result)
  }

  # Handle list format
  if (is.list(sf)) {
    if (learner_id %in% names(sf)) {
      features_list <- sf[[learner_id]]
    } else {
      warning("No feature information found for learner: ", learner_id)
      return(data.table::data.table(feature = character(), frequency = numeric(), selected = logical()))
    }

    if (length(features_list) == 0) {
      return(data.table::data.table(feature = character(), frequency = numeric(), selected = logical()))
    }

    all_features <- unlist(features_list)
    n_folds <- length(features_list)

    freq_table <- table(all_features) / n_folds
    result <- data.table::data.table(
      feature = names(freq_table),
      frequency = as.numeric(freq_table)
    )

    result[, selected := frequency >= min_frequency]
    result <- result[order(-frequency)]

    return(result)
  }

  data.table::data.table(feature = character(), frequency = numeric(), selected = logical())
}


#' @title Extract Consensus Features from Benchmark Result
#' @keywords internal
.extract_consensus_from_bmr <- function(bmr, learner_id, min_frequency) {

  features_per_fold <- list()

  tryCatch({
    for (i in seq_len(bmr$n_resample_results)) {
      rr <- bmr$resample_result(i)

      if (rr$learner$id != learner_id) next

      for (iter in seq_len(rr$iters)) {
        # Try different methods to access the trained learner
        lrn <- NULL
        tryCatch({
          lrn <- rr$learners[[iter]]
        }, error = function(e) {
          tryCatch({
            lrn <<- rr$learner(iter)
          }, error = function(e2) NULL)
        })

        if (is.null(lrn)) next

        # Use the unified feature extraction function
        features <- .extract_features_from_learner(lrn)

        if (!is.null(features) && length(features) > 0) {
          features_per_fold[[length(features_per_fold) + 1]] <- features
        }
      }
    }
  }, error = function(e) {
    warning("Error extracting features from benchmark result: ", e$message)
  })

  if (length(features_per_fold) == 0) {
    warning("Could not extract feature information for learner: ", learner_id)
    return(data.table::data.table(feature = character(), frequency = numeric(), selected = logical()))
  }

  all_features <- unlist(features_per_fold)
  n_folds <- length(features_per_fold)

  freq_table <- table(all_features) / n_folds
  result <- data.table::data.table(
    feature = names(freq_table),
    frequency = as.numeric(freq_table)
  )

  result[, selected := frequency >= min_frequency]
  result <- result[order(-frequency)]

  result
}


#' @title Get Selected Features Per Fold
#'
#' @description
#' Extract the raw list of selected features for each fold from a nested CV result.
#' Useful for debugging feature selection and examining fold-to-fold variability.
#'
#' @param nested_result A NestedCVResult object from BenchmarkService$run()
#' @param learner_id Optional learner ID to filter (default: all learners)
#'
#' @return A data.table with columns: learner_id, fold, n_features, features (list column)
#'
#' @examples
#' \dontrun{
#' result <- benchmark_service$run()
#' features_df <- get_selected_features_per_fold(result)
#' print(features_df)
#'
#' # View features for a specific learner/fold
#' features_df[learner_id == "my_learner" & fold == 1]$features[[1]]
#' }
#'
#' @export
get_selected_features_per_fold <- function(nested_result, learner_id = NULL) {

  if (!inherits(nested_result, "NestedCVResult")) {
    stop("nested_result must be a NestedCVResult object")
  }

  bmr <- nested_result$benchmark_result

  results <- list()

  for (i in seq_len(bmr$n_resample_results)) {
    rr <- bmr$resample_result(i)
    current_learner_id <- rr$learner$id

    # Skip if filtering by learner_id
    if (!is.null(learner_id) && current_learner_id != learner_id) next

    for (iter in seq_len(rr$iters)) {
      lrn <- tryCatch(rr$learners[[iter]], error = function(e) NULL)
      if (is.null(lrn)) next

      features <- .extract_features_from_learner(lrn)

      results[[length(results) + 1]] <- data.table::data.table(
        learner_id = current_learner_id,
        fold = iter,
        n_features = if (!is.null(features)) length(features) else 0L,
        features = list(if (!is.null(features)) features else character(0))
      )
    }
  }

  if (length(results) == 0) {
    return(data.table::data.table(
      learner_id = character(),
      fold = integer(),
      n_features = integer(),
      features = list()
    ))
  }

  data.table::rbindlist(results)
}


#' @title Print Method for Signature Selection Result
#' @param x Result from select_best_signature
#' @param ... Additional arguments (ignored)
#' @export
print.SignatureSelectionResult <- function(x, ...) {
  mode <- attr(x, "mode")
  metric <- attr(x, "metric")
  constraints <- attr(x, "constraints")

  cat("=== Signature Selection Result ===\n\n")
  cat(sprintf("Mode: %s\n", mode))
  cat(sprintf("Metric: %s\n", metric))

  if (!is.null(constraints$auc_min)) cat(sprintf("AUC constraint: >= %.3f\n", constraints$auc_min))
  if (!is.null(constraints$stability_min)) cat(sprintf("Stability constraint: >= %.3f\n", constraints$stability_min))
  if (!is.null(constraints$k_max)) cat(sprintf("Feature count constraint: <= %d\n", constraints$k_max))

  cat(sprintf("\nCandidates evaluated: %d\n", nrow(x)))

  if ("selected" %in% names(x)) {
    selected <- x[selected == TRUE]
    if (nrow(selected) > 0) {
      cat("\n--- Selected Signature ---\n")
      cat(sprintf("Learner: %s\n", selected$learner_id[1]))
      cat(sprintf("Mean %s: %.4f (SE: %.4f)\n", metric, selected$mean_metric[1], selected$se_metric[1]))
      if (!is.na(selected$stability[1])) cat(sprintf("Stability: %.4f\n", selected$stability[1]))
      if (!is.na(selected$mean_k[1])) cat(sprintf("Features: %.1f\n", selected$mean_k[1]))
      if ("score" %in% names(selected)) cat(sprintf("Score: %.4f\n", selected$score[1]))
    }
  }

  invisible(x)
}


#' @title Plot Signature Selection Trade-offs
#'
#' @description
#' Visualize the trade-offs between performance, stability, and parsimony.
#' Shows all candidates with the selected/Pareto-optimal highlighted.
#'
#' @param selection_result Result from select_best_signature()
#' @param x_var Variable for x-axis: "mean_metric", "stability", or "mean_k"
#' @param y_var Variable for y-axis: "mean_metric", "stability", or "mean_k"
#' @param color_var Variable for color: "mean_metric", "stability", "mean_k", or "selected"
#'
#' @return A ggplot2 object
#'
#' @export
plot_signature_tradeoffs <- function(
    selection_result,
    x_var = "stability",
    y_var = "mean_metric",
    color_var = "mean_k"
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_signature_tradeoffs()")
  }

  dt <- data.table::as.data.table(selection_result)

  # Create plot
  p <- ggplot2::ggplot(dt, ggplot2::aes(
    x = .data[[x_var]],
    y = .data[[y_var]]
  ))

  # Color mapping
  if (color_var == "selected") {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(color = selected, size = selected),
      alpha = 0.8
    ) +
    ggplot2::scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) +
    ggplot2::scale_size_manual(values = c("FALSE" = 2, "TRUE" = 5))
  } else {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(color = .data[[color_var]]),
      size = 3, alpha = 0.8
    ) +
    ggplot2::scale_color_viridis_c()

    # Highlight selected
    if ("selected" %in% names(dt)) {
      selected_dt <- dt[selected == TRUE]
      if (nrow(selected_dt) > 0) {
        p <- p + ggplot2::geom_point(
          data = selected_dt,
          color = "red", size = 5, shape = 1, stroke = 2
        )
      }
    }
  }

  # Labels
  labels <- list(
    mean_metric = attr(selection_result, "metric") %||% "Performance",
    stability = "Nogueira Stability Index",
    mean_k = "Number of Features"
  )

  p <- p +
    ggplot2::labs(
      x = labels[[x_var]] %||% x_var,
      y = labels[[y_var]] %||% y_var,
      color = labels[[color_var]] %||% color_var,
      title = "Signature Selection Trade-offs",
      subtitle = sprintf("Mode: %s", attr(selection_result, "mode") %||% "unknown")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  p
}
