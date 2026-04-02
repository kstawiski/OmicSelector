#' @title Cross-Platform Transfer Learning Utilities
#'
#' @description
#' Utilities for adapting feature distributions between measurement platforms
#' before model transfer. This is designed for settings where standard frozen
#' batch correction is not applicable because the source training set contains
#' only one platform level.
#'
#' @details
#' `CrossPlatformAdapter` implements lightweight platform adaptation strategies
#' that can be applied before training on a source platform and predicting on a
#' target platform:
#'
#' - `"rank"`: within-sample percentile ranks across features
#' - `"quantile"`: quantile normalization to the source reference distribution
#' - `"zscore"`: within-sample z-scoring
#' - `"reference"`: subtraction of per-sample reference miRNA summary
#' - `"combat_pooled"`: pooled ComBat correction across source and target
#'
#' The pooled ComBat strategy is semi-supervised because it needs some target
#' labels to preserve biology while correcting platform effects.
#'
#' @name cross-platform
NULL


.as_transfer_matrix <- function(data, arg_name = "data") {
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  checkmate::assert_matrix(
    data,
    mode = "numeric",
    min.rows = 1L,
    min.cols = 1L,
    any.missing = FALSE,
    .var.name = arg_name
  )

  out <- unclass(data)
  storage.mode(out) <- "numeric"

  if (is.null(colnames(out))) {
    colnames(out) <- sprintf("feature_%d", seq_len(ncol(out)))
  }

  out
}


.align_transfer_matrices <- function(source_data, target_data) {
  source_mat <- .as_transfer_matrix(source_data, "source_data")
  target_mat <- .as_transfer_matrix(target_data, "target_data")

  if (!anyDuplicated(colnames(source_mat)) && !anyDuplicated(colnames(target_mat))) {
    if (all(colnames(source_mat) %in% colnames(target_mat))) {
      target_mat <- target_mat[, colnames(source_mat), drop = FALSE]
      return(list(
        source = source_mat,
        target = target_mat,
        features = colnames(source_mat)
      ))
    }
  }

  checkmate::assert_true(
    ncol(source_mat) == ncol(target_mat),
    .var.name = "source_data/target_data"
  )

  if (is.null(colnames(target_mat))) {
    colnames(target_mat) <- colnames(source_mat)
  }

  if (is.null(colnames(source_mat))) {
    colnames(source_mat) <- colnames(target_mat)
  }

  list(
    source = source_mat,
    target = target_mat,
    features = colnames(source_mat)
  )
}


.row_rank_normalize <- function(data) {
  out <- matrix(NA_real_, nrow = nrow(data), ncol = ncol(data), dimnames = dimnames(data))

  for (i in seq_len(nrow(data))) {
    if (ncol(data) == 1L) {
      out[i, 1L] <- 0.5
      next
    }

    ranks <- rank(data[i, ], ties.method = "average")
    out[i, ] <- (ranks - 1) / (ncol(data) - 1)
  }

  out
}


.row_zscore <- function(data) {
  out <- matrix(0, nrow = nrow(data), ncol = ncol(data), dimnames = dimnames(data))

  for (i in seq_len(nrow(data))) {
    row_mean <- mean(data[i, ])
    row_sd <- stats::sd(data[i, ])

    if (!is.finite(row_sd) || row_sd == 0) {
      next
    }

    out[i, ] <- (data[i, ] - row_mean) / row_sd
  }

  out
}


.reference_normalize <- function(data, reference_idx, summary_fun = c("mean", "median")) {
  summary_fun <- match.arg(summary_fun)

  if (length(reference_idx) == 0L) {
    stop("reference normalization requires at least one reference feature")
  }

  offsets <- if (summary_fun == "mean") {
    rowMeans(data[, reference_idx, drop = FALSE])
  } else {
    apply(data[, reference_idx, drop = FALSE], 1L, stats::median)
  }

  sweep(data, 1L, offsets, FUN = "-")
}


.build_quantile_reference <- function(data) {
  sorted <- t(apply(data, 1L, sort))
  colMeans(sorted)
}


.apply_quantile_reference <- function(data, reference_distribution) {
  out <- matrix(NA_real_, nrow = nrow(data), ncol = ncol(data), dimnames = dimnames(data))

  for (i in seq_len(nrow(data))) {
    ord <- order(data[i, ], decreasing = FALSE)
    out[i, ord] <- reference_distribution
  }

  out
}


.coerce_transfer_labels <- function(labels, n, arg_name = "labels") {
  if (is.null(labels)) {
    return(NULL)
  }

  checkmate::assert_atomic_vector(labels, len = n, any.missing = FALSE, .var.name = arg_name)

  if (is.factor(labels)) {
    return(labels)
  }

  if (is.logical(labels)) {
    return(factor(as.integer(labels), levels = c(0L, 1L)))
  }

  if (is.numeric(labels) || is.integer(labels)) {
    uniq <- sort(unique(as.numeric(labels)))
    if (length(uniq) == 2L && all(uniq %in% c(0, 1))) {
      return(factor(as.integer(labels), levels = c(0L, 1L)))
    }
  }

  factor(labels)
}


.compute_cohens_d <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  nx <- length(x)
  ny <- length(y)

  if (nx < 2L || ny < 2L) {
    return(NA_real_)
  }

  sx2 <- stats::var(x)
  sy2 <- stats::var(y)
  pooled_sd <- sqrt(((nx - 1) * sx2 + (ny - 1) * sy2) / (nx + ny - 2))

  if (!is.finite(pooled_sd) || pooled_sd == 0) {
    mean_diff <- mean(x) - mean(y)
    if (isTRUE(all.equal(mean_diff, 0))) {
      return(0)
    }
    return(sign(mean_diff) * Inf)
  }

  (mean(x) - mean(y)) / pooled_sd
}


.compute_binary_auc <- function(scores, labels) {
  n_pos <- sum(labels == 1L)
  n_neg <- sum(labels == 0L)

  if (n_pos == 0L || n_neg == 0L) {
    return(NA_real_)
  }

  ranks <- rank(scores, ties.method = "average")
  u_stat <- sum(ranks[labels == 1L]) - n_pos * (n_pos + 1) / 2

  u_stat / (n_pos * n_neg)
}


.extract_prediction_vector <- function(pred) {
  if (is.list(pred) && !is.data.frame(pred)) {
    for (field in c("prob", "response", "predictions", "prediction", "scores")) {
      if (!is.null(pred[[field]])) {
        return(.extract_prediction_vector(pred[[field]]))
      }
    }
  }

  if (is.data.frame(pred)) {
    pred <- as.matrix(pred)
  }

  if (is.matrix(pred)) {
    if (ncol(pred) == 1L) {
      pred <- pred[, 1L]
    } else {
      pred <- pred[, ncol(pred)]
    }
  }

  if (is.logical(pred)) {
    return(as.integer(pred))
  }

  if (is.factor(pred)) {
    return(as.character(pred))
  }

  pred
}


.train_cross_platform_model <- function(model, x, y) {
  if (inherits(model, "Learner")) {
    trained_model <- model$clone(deep = TRUE)
    y_factor <- .coerce_transfer_labels(y, nrow(x), "source_labels")
    backend <- data.frame(.target = y_factor, as.data.frame(x), check.names = FALSE)
    positive <- if (nlevels(y_factor) == 2L) levels(y_factor)[2L] else NULL
    task <- mlr3::TaskClassif$new(
      id = "cross_platform_source",
      backend = backend,
      target = ".target",
      positive = positive
    )
    trained_model$train(task)
    return(trained_model)
  }

  if (is.function(model)) {
    return(model(x, y))
  }

  model
}


.predict_cross_platform_model <- function(model, newdata) {
  if (inherits(model, "Learner")) {
    train_task <- model$state$train_task
    default_label <- if (!is.null(train_task$positive)) {
      train_task$positive
    } else {
      train_task$class_names[1L]
    }

    backend <- data.frame(
      .target = factor(rep(default_label, nrow(newdata)), levels = train_task$class_names),
      as.data.frame(newdata),
      check.names = FALSE
    )

    task <- mlr3::TaskClassif$new(
      id = "cross_platform_target",
      backend = backend,
      target = ".target",
      positive = train_task$positive
    )

    pred <- model$predict(task)
    if (!is.null(pred$prob)) {
      return(as.numeric(pred$prob[, ncol(pred$prob)]))
    }
    return(as.character(pred$response))
  }

  if (!is.null(model$predict_proba) && is.function(model$predict_proba)) {
    return(.extract_prediction_vector(model$predict_proba(newdata)))
  }

  if (!is.null(model$predict) &&
      is.function(model$predict) &&
      !inherits(model, c("glm", "lm", "train"))) {
    pred <- tryCatch(
      model$predict(newdata),
      error = function(e) NULL
    )
    if (!is.null(pred)) {
      return(.extract_prediction_vector(pred))
    }
  }

  pred <- tryCatch(
    stats::predict(model, newdata = as.data.frame(newdata), type = "response"),
    error = function(e) NULL
  )

  if (is.null(pred)) {
    pred <- stats::predict(model, newdata = as.data.frame(newdata))
  }

  .extract_prediction_vector(pred)
}


.evaluate_transfer_predictions <- function(predictions, labels) {
  truth_factor <- .coerce_transfer_labels(labels, length(labels), "target_labels")
  levels_truth <- levels(truth_factor)

  result <- list(
    n_samples = length(labels),
    levels = levels_truth
  )

  if (length(levels_truth) == 2L) {
    truth_binary <- as.integer(truth_factor == levels_truth[2L])
  } else {
    truth_binary <- NULL
  }

  if (is.numeric(predictions)) {
    result$scores <- as.numeric(predictions)

    if (!is.null(truth_binary)) {
      result$auc <- .compute_binary_auc(result$scores, truth_binary)
    }

    if (all(result$scores >= 0 & result$scores <= 1) && length(levels_truth) == 2L) {
      predicted_labels <- factor(
        ifelse(result$scores >= 0.5, levels_truth[2L], levels_truth[1L]),
        levels = levels_truth
      )
    } else {
      predicted_labels <- NULL
    }
  } else {
    predicted_labels <- factor(as.character(predictions), levels = levels_truth)
  }

  if (!is.null(predicted_labels)) {
    result$predicted_labels <- predicted_labels
    result$accuracy <- mean(predicted_labels == truth_factor)

    if (length(levels_truth) == 2L) {
      pred_binary <- as.integer(predicted_labels == levels_truth[2L])
      tp <- sum(pred_binary == 1L & truth_binary == 1L)
      tn <- sum(pred_binary == 0L & truth_binary == 0L)
      fp <- sum(pred_binary == 1L & truth_binary == 0L)
      fn <- sum(pred_binary == 0L & truth_binary == 1L)

      result$sensitivity <- if ((tp + fn) == 0L) NA_real_ else tp / (tp + fn)
      result$specificity <- if ((tn + fp) == 0L) NA_real_ else tn / (tn + fp)
      result$confusion <- data.table::data.table(
        metric = c("tp", "tn", "fp", "fn"),
        value = c(tp, tn, fp, fn)
      )
    }
  }

  result
}


#' @title CrossPlatformAdapter R6 Class
#'
#' @description
#' An adapter for learning a source-platform representation and applying it to
#' target-platform data before model transfer.
#'
#' @param strategy Adaptation strategy. One of `"rank"`, `"quantile"`,
#'   `"zscore"`, `"reference"`, or `"combat_pooled"`.
#' @param reference_features Optional character vector of reference miRNAs used
#'   when `strategy = "reference"`.
#' @param reference_summary Summary statistic for reference normalization:
#'   `"mean"` or `"median"`.
#' @param source_platform Label used for source samples in pooled ComBat.
#' @param target_platform Label used for target samples in pooled ComBat.
#' @param combat_parametric Logical, use parametric priors in pooled ComBat.
#' @param combat_mean_only Logical, only correct means in pooled ComBat.
#'
#' @return A `CrossPlatformAdapter` object.
#'
#' @examples
#' set.seed(42)
#' src <- matrix(rnorm(100 * 5), 100, 5)
#' tgt <- matrix(rnorm(50 * 5, mean = 2), 50, 5)
#'
#' adapter <- CrossPlatformAdapter$new(strategy = "rank")
#' adapter$fit(src)
#' adapted_tgt <- adapter$adapt(tgt)
#' head(adapted_tgt)
#'
#' @export
CrossPlatformAdapter <- R6::R6Class(
  "CrossPlatformAdapter",
  public = list(
    initialize = function(strategy = c("rank", "quantile", "zscore", "reference", "combat_pooled"),
                          reference_features = NULL,
                          reference_summary = c("mean", "median"),
                          source_platform = "source",
                          target_platform = "target",
                          combat_parametric = TRUE,
                          combat_mean_only = FALSE) {
      private$.strategy <- match.arg(strategy)
      private$.reference_features <- reference_features
      private$.reference_summary <- match.arg(reference_summary)
      private$.source_platform <- source_platform
      private$.target_platform <- target_platform
      private$.combat_parametric <- combat_parametric
      private$.combat_mean_only <- combat_mean_only
      private$.fitted <- FALSE
    },

    #' @description
    #' Learn a source-platform reference representation.
    #'
    #' @param source_data Numeric matrix with samples in rows and features in columns.
    #' @param source_labels Optional labels for the source data. Required for
    #'   `strategy = "combat_pooled"`.
    #'
    #' @return The adapter itself, invisibly.
    fit = function(source_data, source_labels = NULL) {
      source_data <- .as_transfer_matrix(source_data, "source_data")

      private$.feature_names <- colnames(source_data)
      private$.source_raw <- source_data
      private$.source_labels <- .coerce_transfer_labels(
        source_labels,
        nrow(source_data),
        "source_labels"
      )

      if (private$.strategy == "reference") {
        checkmate::assert_character(
          private$.reference_features,
          min.len = 1L,
          unique = TRUE,
          any.missing = FALSE,
          .var.name = "reference_features"
        )

        missing_ref <- setdiff(private$.reference_features, private$.feature_names)
        if (length(missing_ref) > 0L) {
          stop(
            sprintf(
              "reference_features missing from source_data: %s",
              paste(missing_ref, collapse = ", ")
            )
          )
        }

        private$.reference_idx <- match(private$.reference_features, private$.feature_names)
      }

      if (private$.strategy == "quantile") {
        private$.quantile_reference <- .build_quantile_reference(source_data)
      }

      if (private$.strategy != "combat_pooled") {
        private$.source_adapted <- private$.apply_strategy(source_data)
      } else {
        private$.source_adapted <- source_data
      }

      private$.fitted <- TRUE
      invisible(self)
    },

    #' @description
    #' Adapt target-platform data to the learned source representation.
    #'
    #' @param target_data Numeric matrix with samples in rows and features in columns.
    #' @param target_labels Optional target labels. Required when
    #'   `strategy = "combat_pooled"`.
    #'
    #' @return Adapted target matrix with samples in rows and features in columns.
    adapt = function(target_data, target_labels = NULL) {
      if (!private$.fitted) {
        stop("CrossPlatformAdapter must be fitted before adapt().")
      }

      aligned <- .align_transfer_matrices(private$.source_raw, target_data)
      target_mat <- aligned$target

      if (private$.strategy == "combat_pooled") {
        if (!exists("FrozenComBat", mode = "function")) {
          stop(
            "combat_pooled requires FrozenComBat to be available. ",
            "Source 'R/frozen-combat.R' or load the OmicSelector package first."
          )
        }

        if (is.null(private$.source_labels)) {
          stop("combat_pooled requires source_labels during fit().")
        }

        target_labels <- .coerce_transfer_labels(target_labels, nrow(target_mat), "target_labels")
        if (is.null(target_labels)) {
          stop("combat_pooled requires target_labels during adapt().")
        }

        combined <- rbind(private$.source_raw, target_mat)
        batch <- c(
          rep(private$.source_platform, nrow(private$.source_raw)),
          rep(private$.target_platform, nrow(target_mat))
        )
        covariates <- data.frame(
          label = factor(c(as.character(private$.source_labels), as.character(target_labels)))
        )

        frozen <- FrozenComBat$new(
          parametric = private$.combat_parametric,
          mean_only = private$.combat_mean_only
        )

        corrected <- frozen$fit_transform(
          data = combined,
          batch = batch,
          covariates = covariates
        )

        private$.source_adapted <- corrected[seq_len(nrow(private$.source_raw)), , drop = FALSE]
        private$.target_adapted <- corrected[nrow(private$.source_raw) + seq_len(nrow(target_mat)), , drop = FALSE]
      } else {
        private$.target_adapted <- private$.apply_strategy(target_mat)
      }

      private$.target_adapted
    },

    #' @description
    #' Return the stored source matrix.
    #'
    #' @param adapted Logical, return the adapted source matrix if `TRUE`,
    #'   otherwise return the raw source matrix.
    #'
    #' @return A numeric matrix.
    get_source_data = function(adapted = TRUE) {
      if (!private$.fitted) {
        stop("CrossPlatformAdapter not fitted yet.")
      }

      if (adapted) {
        private$.source_adapted
      } else {
        private$.source_raw
      }
    },

    #' @description
    #' Check whether the adapter has been fitted.
    #'
    #' @return Logical scalar.
    is_fitted = function() {
      private$.fitted
    },

    #' @description
    #' Print a compact adapter summary.
    print = function() {
      cat("=== CrossPlatformAdapter ===\n")
      cat(sprintf("  Fitted: %s\n", private$.fitted))
      cat(sprintf("  Strategy: %s\n", private$.strategy))
      if (private$.fitted) {
        cat(sprintf("  Features: %d\n", ncol(private$.source_raw)))
        cat(sprintf("  Source samples: %d\n", nrow(private$.source_raw)))
      }
      invisible(self)
    }
  ),
  private = list(
    .strategy = NULL,
    .reference_features = NULL,
    .reference_summary = "mean",
    .source_platform = "source",
    .target_platform = "target",
    .combat_parametric = TRUE,
    .combat_mean_only = FALSE,
    .fitted = FALSE,
    .feature_names = NULL,
    .reference_idx = integer(),
    .quantile_reference = NULL,
    .source_raw = NULL,
    .source_adapted = NULL,
    .source_labels = NULL,
    .target_adapted = NULL,

    .apply_strategy = function(data) {
      switch(
        private$.strategy,
        rank = .row_rank_normalize(data),
        quantile = .apply_quantile_reference(data, private$.quantile_reference),
        zscore = .row_zscore(data),
        reference = .reference_normalize(
          data,
          reference_idx = private$.reference_idx,
          summary_fun = private$.reference_summary
        ),
        combat_pooled = data,
        stop(sprintf("Unknown strategy: %s", private$.strategy))
      )
    }
  )
)


#' @title Quantify Domain Shift Between Source and Target Platforms
#'
#' @description
#' Computes a per-feature domain-shift summary using the Kolmogorov-Smirnov test
#' and Cohen's d effect size.
#'
#' @param source_data Numeric matrix with source samples in rows and features in columns.
#' @param target_data Numeric matrix with target samples in rows and features in columns.
#'
#' @return A `data.table` with columns `feature`, `ks_stat`, `ks_p`, and
#'   `effect_size_d`.
#'
#' @examples
#' set.seed(42)
#' src <- matrix(rnorm(100 * 5), 100, 5)
#' tgt <- matrix(rnorm(50 * 5, mean = 2), 50, 5)
#'
#' compute_domain_shift(src, tgt)
#'
#' @export
compute_domain_shift <- function(source_data, target_data) {
  aligned <- .align_transfer_matrices(source_data, target_data)
  source_mat <- aligned$source
  target_mat <- aligned$target
  features <- aligned$features

  stats_list <- vector("list", length = ncol(source_mat))

  for (j in seq_len(ncol(source_mat))) {
    ks <- suppressWarnings(
      stats::ks.test(source_mat[, j], target_mat[, j], exact = FALSE)
    )

    stats_list[[j]] <- data.table::data.table(
      feature = features[j],
      ks_stat = unname(ks$statistic),
      ks_p = ks$p.value,
      effect_size_d = .compute_cohens_d(source_mat[, j], target_mat[, j])
    )
  }

  result <- data.table::rbindlist(stats_list)
  result[order(-ks_stat, ks_p)]
}


#' @title Convenience Wrapper for Cross-Platform Transfer
#'
#' @description
#' Fits a cross-platform adapter on the source data, adapts source and target
#' matrices, fits or reuses a prediction model, generates predictions on the
#' target platform, and evaluates performance when target labels are available.
#'
#' @param model A training function, an `mlr3` learner, or an already fitted
#'   model object accepted by [stats::predict()].
#' @param source_data Numeric source matrix with samples in rows and features in columns.
#' @param source_labels Source labels used for model fitting.
#' @param target_data Numeric target matrix with samples in rows and features in columns.
#' @param target_labels Optional target labels for evaluation.
#' @param strategy Adaptation strategy passed to [CrossPlatformAdapter].
#' @param reference_features Optional reference miRNA names for
#'   `strategy = "reference"`.
#' @param reference_summary Summary statistic for reference normalization.
#' @param ... Additional arguments forwarded to [CrossPlatformAdapter].
#'
#' @return A list with the fitted adapter, adapted matrices, fitted model,
#'   predictions, optional evaluation metrics, and domain-shift summaries before
#'   and after adaptation.
#'
#' @examples
#' set.seed(42)
#' src <- matrix(rnorm(100 * 5), 100, 5)
#' tgt <- matrix(rnorm(50 * 5, mean = 1.5), 50, 5)
#' y_src <- rbinom(100, 1, 0.4)
#' y_tgt <- rbinom(50, 1, 0.4)
#'
#' fit_glm <- function(x, y) {
#'   stats::glm(
#'     y ~ .,
#'     data = data.frame(y = as.integer(y), x),
#'     family = stats::binomial()
#'   )
#' }
#'
#' transfer <- cross_platform_transfer(
#'   model = fit_glm,
#'   source_data = src,
#'   source_labels = y_src,
#'   target_data = tgt,
#'   target_labels = y_tgt,
#'   strategy = "rank"
#' )
#'
#' transfer$evaluation$auc
#'
#' @export
cross_platform_transfer <- function(model,
                                    source_data,
                                    source_labels,
                                    target_data,
                                    target_labels = NULL,
                                    strategy = "rank",
                                    reference_features = NULL,
                                    reference_summary = c("mean", "median"),
                                    ...) {
  checkmate::assert_true(length(source_labels) == nrow(source_data), .var.name = "source_labels")
  if (!is.null(target_labels)) {
    checkmate::assert_true(length(target_labels) == nrow(target_data), .var.name = "target_labels")
  }

  adapter <- CrossPlatformAdapter$new(
    strategy = strategy,
    reference_features = reference_features,
    reference_summary = match.arg(reference_summary),
    ...
  )

  adapter$fit(source_data = source_data, source_labels = source_labels)
  adapted_target <- adapter$adapt(target_data = target_data, target_labels = target_labels)
  adapted_source <- adapter$get_source_data(adapted = TRUE)

  fitted_model <- .train_cross_platform_model(model, adapted_source, source_labels)
  predictions <- .predict_cross_platform_model(fitted_model, adapted_target)

  list(
    strategy = strategy,
    adapter = adapter,
    fitted_model = fitted_model,
    source_adapted = adapted_source,
    target_adapted = adapted_target,
    predictions = predictions,
    evaluation = if (!is.null(target_labels)) {
      .evaluate_transfer_predictions(predictions, target_labels)
    } else {
      NULL
    },
    shift_before = compute_domain_shift(source_data, target_data),
    shift_after = compute_domain_shift(adapted_source, adapted_target)
  )
}
