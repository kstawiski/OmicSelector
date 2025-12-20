#' @title Stability-Based Ensemble Selection
#'
#' @description
#' Builds an ensemble of models using stability-based feature selection.
#' The approach:
#' 1. Run N bootstrap resamples with feature selection
#' 2. Compute feature selection frequencies
#' 3. Create nested overlapping subsets based on stability tiers
#' 4. Train K base models on different stability tiers
#' 5. Combine predictions via stability-weighted voting
#'
#' This approach addresses the "reproducibility crisis" in biomarker discovery
#' by prioritizing stable, reproducible feature sets over single-run selections.
#'
#' @references
#' Meinshausen & Buhlmann (2010). Stability Selection. JRSS-B.
#' Nogueira et al. (2018). On the Stability of Feature Selection Algorithms.
#'
#' @importFrom R6 R6Class
#' @importFrom mlr3pipelines %>>%
#' @export
StabilityEnsemble <- R6::R6Class(

  "StabilityEnsemble",

  public = list(

    #' @field n_bootstrap Number of bootstrap resamples
    n_bootstrap = 100,

    #' @field filter Feature selection filter method
    filter = "anova",

    #' @field n_features Number of features to select per bootstrap
    n_features = 50,

    #' @field base_learner Base learner for ensemble members
    base_learner = "glmnet",

    #' @field threshold_mode "adaptive" (target set size) or "fixed" (percentage)
    threshold_mode = "adaptive",

    #' @field target_set_size Target number of stable features (for adaptive mode)
    target_set_size = c(10, 50),

    #' @field fixed_thresholds Thresholds for fixed mode (creates nested subsets)
    fixed_thresholds = c(0.9, 0.7, 0.5),

    #' @field aggregation "weighted" (stability-weighted) or "average" (simple)
    aggregation = "weighted",

    #' @field verbose Print progress messages
    verbose = TRUE,

    #' @field frequencies Feature selection frequencies from last fit
    frequencies = NULL,

    #' @field stable_features List of stable feature sets by tier
    stable_features = NULL,

    #' @field tier_weights Weights for each tier in ensemble
    tier_weights = NULL,

    #' @field fitted_models List of fitted models per tier
    fitted_models = NULL,

    #' @description

    #' Create a new StabilityEnsemble
    #'
    #' @param n_bootstrap Number of bootstrap resamples (default: 100)
    #' @param filter Feature selection filter: "anova", "mrmr", "importance"
    #' @param n_features Features to select per bootstrap (default: 50)
    #' @param base_learner Learner for ensemble: "glmnet", "ranger", "xgboost"
    #' @param threshold_mode "adaptive" or "fixed"
    #' @param target_set_size Target stable set size range for adaptive mode
    #' @param fixed_thresholds Thresholds for fixed mode
    #' @param aggregation "weighted" or "average"
    #' @param verbose Print progress
    #'
    #' @return A new StabilityEnsemble object
    initialize = function(n_bootstrap = 100,
                          filter = "anova",
                          n_features = 50,
                          base_learner = "glmnet",
                          threshold_mode = "adaptive",
                          target_set_size = c(10, 50),
                          fixed_thresholds = c(0.9, 0.7, 0.5),
                          aggregation = "weighted",
                          verbose = TRUE) {

      # Input validation
      stopifnot(
        "n_bootstrap must be >= 10" = n_bootstrap >= 10,
        "filter must be a string" = is.character(filter) && length(filter) == 1,
        "n_features must be > 0" = n_features > 0,
        "base_learner must be a string" = is.character(base_learner) && length(base_learner) == 1,
        "threshold_mode must be 'adaptive' or 'fixed'" =
          threshold_mode %in% c("adaptive", "fixed"),
        "target_set_size must be length 2" = length(target_set_size) == 2,
        "aggregation must be 'weighted' or 'average'" =
          aggregation %in% c("weighted", "average")
      )

      self$n_bootstrap <- n_bootstrap
      self$filter <- filter
      self$n_features <- n_features
      self$base_learner <- base_learner
      self$threshold_mode <- threshold_mode
      self$target_set_size <- target_set_size
      self$fixed_thresholds <- sort(fixed_thresholds, decreasing = TRUE)
      self$aggregation <- aggregation
      self$verbose <- verbose

      invisible(self)
    },

    #' @description
    #' Compute feature selection frequencies from bootstrap resamples
    #'
    #' @param task An mlr3 Task
    #' @param seed Random seed for reproducibility
    #'
    #' @return Named numeric vector of selection frequencies
    compute_frequencies = function(task, seed = NULL) {
      if (!requireNamespace("mlr3", quietly = TRUE)) {
        stop("Package 'mlr3' is required but not installed.", call. = FALSE)
      }
      if (!requireNamespace("mlr3filters", quietly = TRUE)) {
        stop("Package 'mlr3filters' is required but not installed.", call. = FALSE)
      }

      if (!is.null(seed)) set.seed(seed)

      feature_names <- task$feature_names
      n_features_total <- length(feature_names)

      # Track selection counts
      selection_counts <- setNames(rep(0, n_features_total), feature_names)

      # Get filter
      flt <- mlr3filters::flt(self$filter)

      if (self$verbose) {
        message(sprintf("Running %d bootstrap resamples...", self$n_bootstrap))
      }

      for (b in seq_len(self$n_bootstrap)) {
        # Bootstrap sample (with replacement)
        n <- task$nrow
        boot_idx <- sample(n, n, replace = TRUE)

        # Create subset task
        boot_task <- task$clone()
        boot_task$filter(boot_idx)

        # Calculate filter scores
        flt$calculate(boot_task)
        scores <- flt$scores

        # Select top n_features
        n_select <- min(self$n_features, length(scores))
        selected <- names(sort(scores, decreasing = TRUE))[1:n_select]

        # Update counts
        selection_counts[selected] <- selection_counts[selected] + 1

        if (self$verbose && b %% 20 == 0) {
          message(sprintf("  Bootstrap %d/%d", b, self$n_bootstrap))
        }
      }

      # Convert to frequencies
      self$frequencies <- selection_counts / self$n_bootstrap

      if (self$verbose) {
        n_stable <- sum(self$frequencies >= 0.5)
        message(sprintf("  Features with frequency >= 50%%: %d", n_stable))
      }

      self$frequencies
    },

    #' @description
    #' Identify stable feature subsets based on frequency thresholds
    #'
    #' @return List of feature vectors for each stability tier
    identify_stable_subsets = function() {
      if (is.null(self$frequencies)) {
        stop("Must call compute_frequencies() first", call. = FALSE)
      }

      freq <- self$frequencies
      freq_sorted <- sort(freq, decreasing = TRUE)

      if (self$threshold_mode == "adaptive") {
        # Target set size range
        min_size <- self$target_set_size[1]
        max_size <- self$target_set_size[2]

        # Find thresholds that give appropriate set sizes
        # Core tier: smallest set (highest threshold)
        # Balanced tier: medium set
        # Full tier: largest set

        # Create 3 nested tiers
        thresholds <- c()
        for (target in c(min_size, (min_size + max_size) / 2, max_size)) {
          if (target <= length(freq_sorted)) {
            # Threshold is frequency of target-th ranked feature
            thresh <- freq_sorted[min(target, length(freq_sorted))]
            thresholds <- c(thresholds, thresh)
          }
        }
        thresholds <- unique(sort(thresholds, decreasing = TRUE))

        if (self$verbose) {
          message(sprintf("Adaptive thresholds: %s",
                          paste(round(thresholds, 3), collapse = ", ")))
        }

      } else {
        thresholds <- self$fixed_thresholds
      }

      # Create nested overlapping subsets
      # Each tier includes all features from higher tiers (nested structure)
      self$stable_features <- list()
      self$tier_weights <- c()

      cumulative_features <- character(0)

      for (i in seq_along(thresholds)) {
        thresh <- thresholds[i]
        tier_features <- names(freq[freq >= thresh])

        # Nested: include all features from previous tiers
        tier_features <- unique(c(cumulative_features, tier_features))
        cumulative_features <- tier_features

        tier_name <- sprintf("tier_%d", i)
        self$stable_features[[tier_name]] <- tier_features

        # Weight based on average stability of features in this tier
        tier_freq <- freq[tier_features]
        tier_weight <- mean(tier_freq)
        self$tier_weights <- c(self$tier_weights, tier_weight)

        if (self$verbose) {
          message(sprintf("  %s: %d features (avg stability: %.3f)",
                          tier_name, length(tier_features), tier_weight))
        }
      }

      names(self$tier_weights) <- names(self$stable_features)

      # Normalize weights
      self$tier_weights <- self$tier_weights / sum(self$tier_weights)

      self$stable_features
    },

    #' @description
    #' Create mlr3 learner for a specific feature subset
    #'
    #' @param features Character vector of feature names
    #'
    #' @return An mlr3 Learner
    create_tier_learner = function(features) {
      if (!requireNamespace("mlr3", quietly = TRUE)) {
        stop("Package 'mlr3' is required but not installed.", call. = FALSE)
      }
      if (!requireNamespace("mlr3pipelines", quietly = TRUE)) {
        stop("Package 'mlr3pipelines' is required but not installed.", call. = FALSE)
      }

      learner_map <- list(
        "glmnet" = "classif.glmnet",
        "ranger" = "classif.ranger",
        "xgboost" = "classif.xgboost",
        "svm" = "classif.svm",
        "log_reg" = "classif.log_reg"
      )

      if (!self$base_learner %in% names(learner_map)) {
        stop(sprintf("Unknown base_learner: %s", self$base_learner), call. = FALSE)
      }

      # Create feature selector PipeOp
      po_select <- mlr3pipelines::po("select",
        selector = mlr3pipelines::selector_name(features)
      )

      # Create base learner
      base_lrn <- mlr3::lrn(learner_map[[self$base_learner]], predict_type = "prob")

      # Chain
      graph <- po_select %>>% mlr3pipelines::po("learner", base_lrn)
      graph_learner <- mlr3::as_learner(graph)

      graph_learner
    },

    #' @description
    #' Fit the stability ensemble on a task
    #'
    #' @param task An mlr3 Task
    #' @param seed Random seed
    #'
    #' @return The fitted StabilityEnsemble (invisibly)
    fit = function(task, seed = NULL) {
      if (!requireNamespace("mlr3", quietly = TRUE)) {
        stop("Package 'mlr3' is required but not installed.", call. = FALSE)
      }

      if (self$verbose) {
        message("Fitting Stability Ensemble...")
      }

      # Step 1: Compute frequencies
      self$compute_frequencies(task, seed = seed)

      # Step 2: Identify stable subsets
      self$identify_stable_subsets()

      # Step 3: Fit a model for each tier
      self$fitted_models <- list()

      for (tier_name in names(self$stable_features)) {
        features <- self$stable_features[[tier_name]]

        if (length(features) == 0) {
          warning(sprintf("Tier %s has no features, skipping", tier_name))
          next
        }

        if (self$verbose) {
          message(sprintf("  Fitting %s model on %d features...",
                          tier_name, length(features)))
        }

        learner <- self$create_tier_learner(features)
        learner$train(task)

        self$fitted_models[[tier_name]] <- learner
      }

      if (self$verbose) {
        message("Stability Ensemble fitting complete.")
      }

      invisible(self)
    },

    #' @description
    #' Predict on new data using the ensemble
    #'
    #' @param task An mlr3 Task (or data.frame with same features)
    #'
    #' @return A data.table with predicted probabilities
    predict = function(task) {
      if (is.null(self$fitted_models) || length(self$fitted_models) == 0) {
        stop("Must call fit() first", call. = FALSE)
      }

      # Get predictions from each tier
      tier_probs <- list()

      for (tier_name in names(self$fitted_models)) {
        model <- self$fitted_models[[tier_name]]
        pred <- model$predict(task)

        # Extract probability for positive class
        prob_matrix <- pred$prob
        tier_probs[[tier_name]] <- prob_matrix
      }

      # Aggregate predictions
      if (self$aggregation == "weighted") {
        # Stability-weighted average
        weights <- self$tier_weights[names(tier_probs)]
        weights <- weights / sum(weights)  # Normalize

        # Weighted sum
        combined <- Reduce(`+`, mapply(function(p, w) p * w,
                                        tier_probs, weights,
                                        SIMPLIFY = FALSE))
      } else {
        # Simple average
        combined <- Reduce(`+`, tier_probs) / length(tier_probs)
      }

      combined
    },

    #' @description
    #' Get feature importance based on selection frequencies
    #'
    #' @param top_n Number of top features to return (default: all)
    #'
    #' @return Named vector of feature frequencies, sorted descending
    get_feature_importance = function(top_n = NULL) {
      if (is.null(self$frequencies)) {
        stop("Must call compute_frequencies() or fit() first", call. = FALSE)
      }

      freq_sorted <- sort(self$frequencies, decreasing = TRUE)

      if (!is.null(top_n) && top_n < length(freq_sorted)) {
        freq_sorted <- freq_sorted[1:top_n]
      }

      freq_sorted
    },

    #' @description
    #' Get summary of the ensemble configuration and results
    #'
    #' @return A list with ensemble summary
    get_summary = function() {
      list(
        n_bootstrap = self$n_bootstrap,
        filter = self$filter,
        n_features_per_bootstrap = self$n_features,
        base_learner = self$base_learner,
        threshold_mode = self$threshold_mode,
        aggregation = self$aggregation,
        n_tiers = length(self$stable_features),
        tier_sizes = sapply(self$stable_features, length),
        tier_weights = self$tier_weights,
        top_10_features = self$get_feature_importance(10)
      )
    },

    #' @description
    #' Print summary
    print = function() {
      cat("Stability-Based Ensemble\n")
      cat("========================\n\n")

      cat(sprintf("Bootstrap resamples: %d\n", self$n_bootstrap))
      cat(sprintf("Filter method: %s\n", self$filter))
      cat(sprintf("Features per resample: %d\n", self$n_features))
      cat(sprintf("Base learner: %s\n", self$base_learner))
      cat(sprintf("Threshold mode: %s\n", self$threshold_mode))
      cat(sprintf("Aggregation: %s\n\n", self$aggregation))

      if (!is.null(self$stable_features)) {
        cat("Stability Tiers:\n")
        for (tier_name in names(self$stable_features)) {
          n_feat <- length(self$stable_features[[tier_name]])
          weight <- self$tier_weights[[tier_name]]
          cat(sprintf("  %s: %d features (weight: %.3f)\n",
                      tier_name, n_feat, weight))
        }
        cat("\n")
      }

      if (!is.null(self$frequencies)) {
        top5 <- self$get_feature_importance(5)
        cat("Top 5 Most Stable Features:\n")
        for (i in seq_along(top5)) {
          cat(sprintf("  %d. %s (%.1f%%)\n",
                      i, names(top5)[i], top5[i] * 100))
        }
      }

      invisible(self)
    }
  )
)


#' @title Create a Stability Ensemble with Presets
#'
#' @description
#' Convenience function to create a StabilityEnsemble with common presets.
#'
#' @param preset One of "default", "fast", "thorough"
#' @param ... Additional arguments passed to StabilityEnsemble$new()
#'
#' @return A StabilityEnsemble object
#'
#' @examples
#' \dontrun{
#' # Default: 100 bootstraps, glmnet base learner
#' ensemble <- create_stability_ensemble("default")
#'
#' # Fast: 50 bootstraps for quick exploration
#' ensemble <- create_stability_ensemble("fast")
#'
#' # Thorough: 500 bootstraps for publication-quality results
#' ensemble <- create_stability_ensemble("thorough")
#'
#' # Fit and predict
#' ensemble$fit(task)
#' probs <- ensemble$predict(task)
#' }
#'
#' @export
create_stability_ensemble <- function(preset = "default", ...) {
  presets <- list(
    default = list(
      n_bootstrap = 100,
      filter = "anova",
      n_features = 50,
      base_learner = "glmnet",
      threshold_mode = "adaptive",
      target_set_size = c(10, 50),
      aggregation = "weighted"
    ),
    fast = list(
      n_bootstrap = 50,
      filter = "anova",
      n_features = 30,
      base_learner = "glmnet",
      threshold_mode = "adaptive",
      target_set_size = c(10, 30),
      aggregation = "average"
    ),
    thorough = list(
      n_bootstrap = 500,
      filter = "mrmr",
      n_features = 100,
      base_learner = "glmnet",
      threshold_mode = "adaptive",
      target_set_size = c(20, 100),
      aggregation = "weighted"
    ),
    ranger = list(
      n_bootstrap = 100,
      filter = "importance",
      n_features = 50,
      base_learner = "ranger",
      threshold_mode = "adaptive",
      target_set_size = c(10, 50),
      aggregation = "weighted"
    )
  )

  if (!preset %in% names(presets)) {
    stop(sprintf("Unknown preset: %s. Available: %s",
                 preset, paste(names(presets), collapse = ", ")), call. = FALSE)
  }

  args <- modifyList(presets[[preset]], list(...))
  do.call(StabilityEnsemble$new, args)
}
