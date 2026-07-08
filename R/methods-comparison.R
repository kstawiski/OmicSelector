#' @title Methods Comparison Utilities
#'
#' @description
#' Helpers for comparing biomarker preprocessing pipelines on the same cohort
#' under leakage-free cross-validation. Includes:
#' \itemize{
#'   \item \code{os_log_transform_adaptive()} — per-feature pseudocount log2
#'         transform that survives mixtures of raw counts and pre-normalized
#'         ratios.
#'   \item \code{os_plate_median_frozen} — R6 class for fold-frozen plate-median
#'         (or batch-median) centering, simpler and more deployable than ComBat
#'         when only a location shift is suspected.
#'   \item \code{os_paired_delong()} — paired DeLong test wrapper for AUC
#'         comparison of two scoring functions on the same cohort.
#'   \item \code{os_oof_pipeline_compare()} — grouped-CV driver that runs each
#'         pipeline train-fold-frozen, pools out-of-fold predictions per seed,
#'         and reports paired DeLong tests vs a reference pipeline.
#' }
#'
#' @details
#' These functions exist to make the "is method X better than the current
#' pipeline?" question answerable in a single call instead of an ad-hoc script.
#' They are meant to compose with existing OmicSelector primitives such as
#' \code{\link{clr_transform}}, \code{\link{FrozenComBat}}, \code{\link{ws_rank}},
#' \code{\link{os_panel_null_benchmark}}, and \code{\link{os_clustered_bootstrap_auc}}.
#'
#' @references
#' DeLong ER, DeLong DM, Clarke-Pearson DL. (1988). Comparing the areas under
#' two or more correlated receiver operating characteristic curves: a
#' nonparametric approach. \emph{Biometrics}, 44(3), 837-845.
#'
#' Robin X, Turck N, Hainard A, et al. (2011). pROC: an open-source package
#' for R and S+ to analyze and compare ROC curves. \emph{BMC Bioinformatics},
#' 12, 77.
#'
#' Stouffer SA, Suchman EA, DeVinney LC, Star SA, Williams RM. (1949).
#' \emph{The American Soldier: Adjustment During Army Life}. Princeton
#' University Press.
#'
#' @name methods-comparison
NULL


# =============================================================================
# 1. Adaptive pseudocount log2 transform
# =============================================================================

#' @title Adaptive Pseudocount log2 Transform
#'
#' @description
#' Compute \code{log2(x + pseudocount)} where the pseudocount is chosen
#' per-feature to avoid distorting the dynamic range. Critical when the same
#' analysis matrix mixes raw counts (typical scale 1-10000) with pre-normalized
#' ratios (typical scale 0.001-0.1) — adding a fixed +1 to a 0.003 ratio
#' compresses its log2 to near zero and erases the signal.
#'
#' @param x Numeric matrix or data.frame with samples in rows and features in
#'   columns. Must contain non-negative values; \code{NA} is propagated.
#' @param method One of \code{"half_min"}, \code{"fixed"}, or \code{"epsilon"}.
#'   \itemize{
#'     \item \code{"half_min"}: per-feature pseudocount = 0.5 * smallest
#'           positive value in that column. Falls back to 1e-6 if a column has
#'           no positive values. Default.
#'     \item \code{"fixed"}: same pseudocount across features, given by
#'           \code{value} (default 1.0). Backward-compatible with the
#'           \code{log2(x + 1)} convention.
#'     \item \code{"epsilon"}: pseudocount = 1e-6 for all features.
#'   }
#' @param value Numeric scalar; the pseudocount used when \code{method = "fixed"}.
#'   Ignored otherwise. Default 1.0.
#'
#' @return A numeric matrix with the same dimensions and dimnames as \code{x},
#'   carrying an attribute \code{"pseudocount"} (length-\code{ncol(x)} vector)
#'   recording the per-feature pseudocount used. The same pseudocount vector
#'   should be reused on held-out data to avoid leakage.
#'
#' @examples
#' x <- matrix(c(0.001, 0.05, 0.13,
#'               2,     30,   1500), nrow = 2, byrow = TRUE)
#' colnames(x) <- c("ratio_a", "ratio_b", "count_c")
#' os_log_transform_adaptive(x, method = "half_min")
#' os_log_transform_adaptive(x, method = "fixed", value = 1)
#'
#' @export
os_log_transform_adaptive <- function(x,
                                      method = c("half_min", "fixed", "epsilon"),
                                      value = 1.0) {
  method <- match.arg(method)
  checkmate::assert_number(value, lower = 0, finite = TRUE)
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  checkmate::assert_matrix(x, mode = "numeric")
  if (any(x < 0, na.rm = TRUE)) {
    stop("os_log_transform_adaptive(): negative values are not allowed.",
         call. = FALSE)
  }

  n_feat <- ncol(x)
  pseudo <- switch(
    method,
    half_min = vapply(seq_len(n_feat), function(j) {
      col <- x[, j]
      pos <- col[!is.na(col) & col > 0]
      if (length(pos) == 0L) return(1e-6)
      0.5 * min(pos)
    }, numeric(1)),
    fixed   = rep(value, n_feat),
    epsilon = rep(1e-6, n_feat)
  )
  names(pseudo) <- colnames(x)

  pseudo_mat <- matrix(pseudo, nrow = nrow(x), ncol = n_feat, byrow = TRUE)
  out <- log2(x + pseudo_mat)
  attr(out, "pseudocount") <- pseudo
  out
}


# =============================================================================
# 2. Frozen plate-median centering (R6, mirrors FrozenComBat API)
# =============================================================================

#' @title Frozen Plate-Median Correction
#'
#' @description
#' R6 class for fold-frozen per-plate (or per-batch) location-shift correction.
#' Subtracts the per-plate median per feature, then adds back the grand median
#' of training-plate medians, so each plate is centered without changing the
#' overall mean. Optionally divides by a per-plate scale (IQR or MAD) before
#' re-multiplying by the grand scale.
#'
#' @details
#' This is intentionally simpler than \code{\link{FrozenComBat}}: it assumes
#' a pure location-shift (or location-and-scale) batch model and is
#' single-plate deployable. For ddPCR-style data where plate effects are
#' largely additive in log-space this is often the right tool.
#'
#' Use \code{\link{FrozenComBat}} when batch effects are believed to interact
#' with biological covariates and an empirical-Bayes shrinkage is desirable.
#'
#' @section Public methods:
#' \describe{
#'   \item{\code{$new(scale)}}{Constructor; \code{scale} ∈ \{"none", "iqr", "mad"\}.}
#'   \item{\code{$fit(data, plate)}}{Learn per-plate medians (and optional scales)
#'         from training data only.}
#'   \item{\code{$transform(data, plate)}}{Apply frozen parameters to new data.
#'         Plates not seen during fit fall back to grand-median centering and
#'         emit a warning.}
#'   \item{\code{$fit_transform(data, plate)}}{Convenience for fit + transform.}
#'   \item{\code{$is_fitted()}}{Logical.}
#'   \item{\code{$get_plate_levels()}}{Character vector of training plate labels.}
#' }
#'
#' @examples
#' set.seed(7)
#' plate <- rep(c("P1", "P2", "P3"), each = 10)
#' x <- matrix(rnorm(60, mean = 0, sd = 0.5), nrow = 30)
#' x[plate == "P2", ] <- x[plate == "P2", ] + 3
#' x[plate == "P3", ] <- x[plate == "P3", ] - 2
#' fc <- os_plate_median_frozen$new()
#' fc$fit(x[1:24, ], plate[1:24])
#' corrected <- fc$transform(x[25:30, ], plate[25:30])
#'
#' @export
os_plate_median_frozen <- R6::R6Class(
  "os_plate_median_frozen",
  public = list(

    #' @description Construct a new os_plate_median_frozen object.
    #' @param scale One of \code{"none"} (default), \code{"iqr"}, or \code{"mad"}.
    initialize = function(scale = c("none", "iqr", "mad")) {
      private$.scale <- match.arg(scale)
      private$.fitted <- FALSE
    },

    #' @description Fit per-plate centering parameters on training data.
    #' @param data Numeric matrix; samples in rows, features in columns.
    #' @param plate Character or factor vector; length \code{nrow(data)}.
    fit = function(data, plate) {
      checkmate::assert_matrix(data, mode = "numeric", min.rows = 2L, min.cols = 1L)
      checkmate::assert_true(length(plate) == nrow(data))
      plate <- as.character(plate)
      if (anyNA(plate)) stop("os_plate_median_frozen: NA in plate vector.", call. = FALSE)
      levels_plate <- sort(unique(plate))
      if (length(levels_plate) < 1L) {
        stop("os_plate_median_frozen: at least one plate level required.", call. = FALSE)
      }

      n_feat <- ncol(data)
      med_mat <- matrix(NA_real_, nrow = length(levels_plate), ncol = n_feat,
                        dimnames = list(levels_plate, colnames(data)))
      scl_mat <- matrix(1.0, nrow = length(levels_plate), ncol = n_feat,
                        dimnames = list(levels_plate, colnames(data)))
      for (i in seq_along(levels_plate)) {
        rows <- plate == levels_plate[i]
        sub <- data[rows, , drop = FALSE]
        med_mat[i, ] <- apply(sub, 2L, stats::median, na.rm = TRUE)
        if (private$.scale == "iqr") {
          scl_mat[i, ] <- pmax(apply(sub, 2L, stats::IQR, na.rm = TRUE), 1e-8)
        } else if (private$.scale == "mad") {
          scl_mat[i, ] <- pmax(apply(sub, 2L, stats::mad, na.rm = TRUE), 1e-8)
        }
      }

      private$.plate_levels <- levels_plate
      private$.plate_medians <- med_mat
      private$.plate_scales  <- scl_mat
      private$.grand_median  <- apply(med_mat, 2L, stats::median, na.rm = TRUE)
      private$.grand_scale   <- if (private$.scale == "none") {
        rep(1.0, n_feat)
      } else {
        apply(scl_mat, 2L, stats::median, na.rm = TRUE)
      }
      private$.feature_names <- colnames(data)
      private$.fitted <- TRUE
      invisible(self)
    },

    #' @description Apply frozen correction to new data.
    #' @param data Numeric matrix.
    #' @param plate Character or factor vector.
    transform = function(data, plate) {
      if (!private$.fitted) {
        stop("os_plate_median_frozen$transform() called before $fit().", call. = FALSE)
      }
      checkmate::assert_matrix(data, mode = "numeric", min.rows = 1L)
      checkmate::assert_true(length(plate) == nrow(data))
      plate <- as.character(plate)

      out <- data
      for (i in seq_len(nrow(data))) {
        p <- plate[i]
        if (p %in% private$.plate_levels) {
          med <- private$.plate_medians[p, ]
          scl <- private$.plate_scales[p, ]
        } else {
          warning(sprintf(
            "os_plate_median_frozen: plate '%s' not seen during fit; falling back to grand-median centering.",
            p), call. = FALSE)
          med <- private$.grand_median
          scl <- private$.grand_scale
        }
        out[i, ] <- (data[i, ] - med) / scl * private$.grand_scale + private$.grand_median
      }
      out
    },

    #' @description Fit and transform in one call.
    #' @param data Numeric matrix.
    #' @param plate Plate vector.
    fit_transform = function(data, plate) {
      self$fit(data, plate)
      self$transform(data, plate)
    },

    #' @description Logical: has the object been fitted?
    is_fitted = function() private$.fitted,

    #' @description Character vector of plate labels seen during fit.
    get_plate_levels = function() {
      if (!private$.fitted) stop("not fitted", call. = FALSE)
      private$.plate_levels
    }
  ),
  private = list(
    .scale = "none",
    .fitted = FALSE,
    .plate_levels = NULL,
    .plate_medians = NULL,
    .plate_scales = NULL,
    .grand_median = NULL,
    .grand_scale = NULL,
    .feature_names = NULL
  )
)


# =============================================================================
# 3. Paired DeLong test wrapper
# =============================================================================

#' @title Paired DeLong Test for AUC Comparison
#'
#' @description
#' Compare two prediction-score vectors on the same labels using the paired
#' DeLong test (Robin et al. 2011 implementation in \code{\link[pROC]{roc.test}}).
#'
#' @param y_true Integer or logical vector of true labels (0/1).
#' @param scores_a Numeric vector of scores from pipeline A; same length as
#'   \code{y_true}.
#' @param scores_b Numeric vector of scores from pipeline B; same length as
#'   \code{y_true}.
#' @param alternative One of \code{"two.sided"} (default), \code{"greater"}
#'   (A's AUC > B's), or \code{"less"}.
#'
#' @return A named list with elements \code{auc_a}, \code{auc_b}, \code{delta}
#'   (= auc_a - auc_b), \code{z}, \code{p_value}, \code{ci_delta} (95% CI of
#'   the AUC difference; \code{NA} if pROC does not return one), \code{n_pos},
#'   \code{n_neg}, \code{alternative}, and \code{method} (= "DeLong").
#'
#' @examples
#' set.seed(7)
#' y <- rep(c(0, 1), each = 50)
#' a <- y + rnorm(100, sd = 0.5)            # informative
#' b <- rnorm(100)                          # random
#' os_paired_delong(y, a, b, alternative = "greater")
#'
#' @export
os_paired_delong <- function(y_true, scores_a, scores_b,
                             alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  checkmate::assert_atomic_vector(y_true)
  checkmate::assert_numeric(scores_a, len = length(y_true))
  checkmate::assert_numeric(scores_b, len = length(y_true))
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("os_paired_delong() requires the 'pROC' package. Install it first.", call. = FALSE)
  }

  y_true <- as.integer(as.logical(y_true))
  if (length(unique(stats::na.omit(y_true))) < 2L) {
    stop("os_paired_delong(): y_true must contain both classes.", call. = FALSE)
  }

  roc_a <- pROC::roc(response = y_true, predictor = scores_a, quiet = TRUE,
                     direction = "<")
  roc_b <- pROC::roc(response = y_true, predictor = scores_b, quiet = TRUE,
                     direction = "<")

  test_result <- tryCatch(
    pROC::roc.test(roc_a, roc_b, method = "delong", paired = TRUE,
                   alternative = alternative),
    error = function(e) e
  )

  auc_a <- as.numeric(pROC::auc(roc_a))
  auc_b <- as.numeric(pROC::auc(roc_b))
  delta <- auc_a - auc_b

  if (inherits(test_result, "error")) {
    return(list(
      auc_a = auc_a, auc_b = auc_b, delta = delta,
      z = NA_real_, p_value = NA_real_,
      ci_delta = c(NA_real_, NA_real_),
      n_pos = sum(y_true == 1L), n_neg = sum(y_true == 0L),
      alternative = alternative, method = "DeLong",
      error = conditionMessage(test_result)
    ))
  }

  ci_delta <- if (!is.null(test_result$conf.int)) {
    as.numeric(test_result$conf.int)
  } else {
    c(NA_real_, NA_real_)
  }

  list(
    auc_a       = auc_a,
    auc_b       = auc_b,
    delta       = delta,
    z           = unname(test_result$statistic),
    p_value     = test_result$p.value,
    ci_delta    = ci_delta,
    n_pos       = sum(y_true == 1L),
    n_neg       = sum(y_true == 0L),
    alternative = alternative,
    method      = "DeLong"
  )
}


# =============================================================================
# 4. OOF pipeline comparison driver
# =============================================================================

#' @title Out-of-Fold Pipeline Comparison
#'
#' @description
#' Run multiple preprocessing-and-scoring pipelines on the same cohort under
#' grouped cross-validation, pool out-of-fold predictions per seed, and report
#' AUCs with clustered bootstrap CIs plus paired DeLong tests vs a reference
#' pipeline.
#'
#' @param y Integer or logical vector of labels (0/1), length \code{N}.
#' @param groups Vector of grouping IDs, length \code{N}; samples within the
#'   same group never appear on opposite sides of a fold split.
#' @param pipelines Named list of pipeline functions. Each function must have
#'   signature
#'   \code{function(X_train, y_train, X_test, batch_train, batch_test,
#'                  plate_train, plate_test) -> list(scores_train, scores_test)}.
#'   The function is responsible for any preprocessing, fold-frozen parameter
#'   estimation, and learner fitting; it returns numeric prediction scores
#'   for the train rows (informational) and test rows (used).
#' @param X Either a single feature matrix shared across all pipelines, or a
#'   named list of matrices keyed by pipeline name. When a list is supplied
#'   its names must match \code{names(pipelines)} exactly.
#' @param batch Optional character/factor vector, length \code{N}, passed to
#'   each pipeline. Use for ComBat-style corrections.
#' @param plate Optional character/factor vector, length \code{N}, passed to
#'   each pipeline. Use for plate-median corrections.
#' @param n_folds Integer; number of CV folds (default 5).
#' @param seeds Integer vector of seeds (default \code{c(7L, 42L, 2026L)}).
#' @param reference Character; name of the pipeline used as comparator in
#'   pairwise DeLong tests. Set to \code{NULL} to skip pairwise tests.
#' @param verbose Logical; print fold-level progress.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{per_seed_auc}}{data.frame: pipeline, seed, auc, ci_lo, ci_hi.}
#'     \item{\code{summary}}{data.frame: pipeline, mean_auc, mean_ci_lo,
#'           mean_ci_hi (across seeds).}
#'     \item{\code{oof_predictions}}{Named list of length \code{length(seeds)}
#'           — each element is a named list of OOF score vectors keyed by
#'           pipeline.}
#'     \item{\code{pairwise_delong}}{data.frame of paired DeLong results vs
#'           \code{reference} per seed (NULL if \code{reference} is NULL).}
#'     \item{\code{pairwise_summary}}{data.frame of meta-combined deltas and
#'           Stouffer-combined p-values across seeds (NULL if reference NULL).}
#'   }
#'
#' @export
os_oof_pipeline_compare <- function(y,
                                    groups,
                                    pipelines,
                                    X = NULL,
                                    batch = NULL,
                                    plate = NULL,
                                    n_folds = 5L,
                                    seeds = c(7L, 42L, 2026L),
                                    reference = NULL,
                                    verbose = FALSE) {
  checkmate::assert_atomic_vector(y)
  checkmate::assert_atomic_vector(groups, len = length(y))
  checkmate::assert_list(pipelines, min.len = 1L, names = "named")
  pipe_names <- names(pipelines)
  for (nm in pipe_names) {
    if (!is.function(pipelines[[nm]])) {
      stop("Pipeline '", nm, "' is not a function.", call. = FALSE)
    }
  }
  if (!is.null(reference)) {
    if (!reference %in% pipe_names) {
      stop("reference '", reference, "' is not in names(pipelines).", call. = FALSE)
    }
  }
  N <- length(y)
  if (is.null(X)) X <- matrix(0, nrow = N, ncol = 1)  # dummy; pipelines may ignore
  if (is.matrix(X) || is.data.frame(X)) {
    X_per_pipe <- setNames(replicate(length(pipe_names), X, simplify = FALSE), pipe_names)
  } else if (is.list(X)) {
    if (!setequal(names(X), pipe_names)) {
      stop("When X is a list, names(X) must equal names(pipelines).", call. = FALSE)
    }
    X_per_pipe <- X
  } else {
    stop("X must be a matrix, data.frame, or named list of either.", call. = FALSE)
  }
  if (!is.null(batch)) checkmate::assert_atomic_vector(batch, len = N)
  if (!is.null(plate)) checkmate::assert_atomic_vector(plate, len = N)
  y_int <- as.integer(as.logical(y))

  per_seed_rows <- list()
  oof_predictions <- list()
  pairwise_rows <- list()

  for (seed in seeds) {
    set.seed(seed)
    folds <- .grouped_stratified_folds(y_int, groups, n_folds = n_folds, seed = seed)
    if (verbose) message(sprintf("[seed=%d] folds built (sizes=%s)", seed,
                                 paste(vapply(folds, length, integer(1)), collapse = ",")))

    oof_for_seed <- setNames(replicate(length(pipe_names),
                                       rep(NA_real_, N), simplify = FALSE),
                             pipe_names)
    for (k in seq_along(folds)) {
      te <- folds[[k]]
      tr <- setdiff(seq_len(N), te)
      bt <- if (is.null(batch)) NULL else batch[tr]
      be <- if (is.null(batch)) NULL else batch[te]
      pt <- if (is.null(plate)) NULL else plate[tr]
      pe <- if (is.null(plate)) NULL else plate[te]
      for (nm in pipe_names) {
        Xp <- X_per_pipe[[nm]]
        Xtr <- if (is.matrix(Xp)) Xp[tr, , drop = FALSE] else as.matrix(Xp)[tr, , drop = FALSE]
        Xte <- if (is.matrix(Xp)) Xp[te, , drop = FALSE] else as.matrix(Xp)[te, , drop = FALSE]
        out <- tryCatch(
          pipelines[[nm]](Xtr, y_int[tr], Xte, bt, be, pt, pe),
          error = function(e) {
            warning(sprintf("Pipeline '%s' failed at seed=%d fold=%d: %s",
                            nm, seed, k, conditionMessage(e)), call. = FALSE)
            list(train = rep(NA_real_, length(tr)), test = rep(NA_real_, length(te)))
          }
        )
        if (!is.list(out) || !all(c("test") %in% names(out))) {
          stop("Pipeline '", nm, "' must return list(train=..., test=...).", call. = FALSE)
        }
        scores_te <- as.numeric(out$test)
        if (length(scores_te) != length(te)) {
          stop(sprintf("Pipeline '%s' returned %d test scores; expected %d.",
                       nm, length(scores_te), length(te)), call. = FALSE)
        }
        oof_for_seed[[nm]][te] <- scores_te
      }
      if (verbose) message(sprintf("  fold %d/%d done (n_test=%d)", k, length(folds), length(te)))
    }
    oof_predictions[[as.character(seed)]] <- oof_for_seed

    # Per-pipeline AUC + clustered bootstrap CI on OOF predictions.
    for (nm in pipe_names) {
      sc <- oof_for_seed[[nm]]
      keep <- !is.na(sc) & !is.na(y_int)
      if (sum(keep) < 10L || length(unique(y_int[keep])) < 2L) {
        per_seed_rows[[length(per_seed_rows) + 1L]] <- data.frame(
          pipeline = nm, seed = seed, auc = NA_real_,
          ci_lo = NA_real_, ci_hi = NA_real_, n = sum(keep),
          stringsAsFactors = FALSE
        )
        next
      }
      auc_val <- .oms_auc_simple(y_int[keep], sc[keep])
      ci <- tryCatch(
        os_clustered_bootstrap_auc(y_int[keep], sc[keep], cluster_id = groups[keep],
                                   n_boot = 2000L, alpha = 0.05, seed = 42L),
        error = function(e) list(ci_lower = NA_real_, ci_upper = NA_real_)
      )
      ci_v <- c(if (is.null(ci$ci_lower)) NA_real_ else ci$ci_lower,
                if (is.null(ci$ci_upper)) NA_real_ else ci$ci_upper)
      per_seed_rows[[length(per_seed_rows) + 1L]] <- data.frame(
        pipeline = nm, seed = seed, auc = auc_val,
        ci_lo = ci_v[1], ci_hi = ci_v[2], n = sum(keep),
        stringsAsFactors = FALSE
      )
    }

    # Paired DeLong vs reference, per seed.
    if (!is.null(reference)) {
      ref_scores <- oof_for_seed[[reference]]
      for (nm in setdiff(pipe_names, reference)) {
        sc <- oof_for_seed[[nm]]
        keep <- !is.na(sc) & !is.na(ref_scores) & !is.na(y_int)
        if (sum(keep) < 10L || length(unique(y_int[keep])) < 2L) {
          next
        }
        dl <- tryCatch(
          os_paired_delong(y_int[keep], sc[keep], ref_scores[keep],
                           alternative = "two.sided"),
          error = function(e) list(auc_a = NA, auc_b = NA, delta = NA,
                                   z = NA, p_value = NA,
                                   ci_delta = c(NA, NA),
                                   n_pos = NA, n_neg = NA,
                                   alternative = "two.sided", method = "DeLong")
        )
        pairwise_rows[[length(pairwise_rows) + 1L]] <- data.frame(
          pipeline = nm, reference = reference, seed = seed,
          auc_pipeline = dl$auc_a, auc_reference = dl$auc_b,
          delta = dl$delta, z = dl$z, p_value = dl$p_value,
          ci_delta_lo = dl$ci_delta[1], ci_delta_hi = dl$ci_delta[2],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  per_seed_df <- do.call(rbind, per_seed_rows)
  summary_df <- do.call(rbind, lapply(split(per_seed_df, per_seed_df$pipeline), function(df) {
    data.frame(pipeline = unique(df$pipeline),
               mean_auc = mean(df$auc, na.rm = TRUE),
               sd_auc   = stats::sd(df$auc, na.rm = TRUE),
               mean_ci_lo = mean(df$ci_lo, na.rm = TRUE),
               mean_ci_hi = mean(df$ci_hi, na.rm = TRUE),
               n_seeds  = sum(!is.na(df$auc)),
               stringsAsFactors = FALSE)
  }))
  rownames(summary_df) <- NULL
  summary_df <- summary_df[order(-summary_df$mean_auc), , drop = FALSE]

  pairwise_df <- if (length(pairwise_rows) > 0L) do.call(rbind, pairwise_rows) else NULL
  pairwise_summary <- if (!is.null(pairwise_df)) {
    do.call(rbind, lapply(split(pairwise_df, pairwise_df$pipeline), function(df) {
      # Stouffer-combined two-sided p across seeds: combine via signed z.
      # For two-sided we work in the |z| scale; preserve sign by majority.
      z <- df$z
      keep <- !is.na(z)
      if (sum(keep) == 0L) {
        return(data.frame(pipeline = unique(df$pipeline),
                          mean_delta = NA_real_,
                          stouffer_z = NA_real_,
                          stouffer_p = NA_real_,
                          n_seeds = 0L, stringsAsFactors = FALSE))
      }
      z <- z[keep]
      sgn <- sign(mean(df$delta[keep], na.rm = TRUE))
      stf_z <- sum(abs(z)) / sqrt(length(z))
      stf_z <- sgn * stf_z
      stf_p <- 2 * stats::pnorm(abs(stf_z), lower.tail = FALSE)
      data.frame(pipeline = unique(df$pipeline),
                 mean_delta = mean(df$delta[keep], na.rm = TRUE),
                 stouffer_z = stf_z,
                 stouffer_p = stf_p,
                 n_seeds = length(z),
                 stringsAsFactors = FALSE)
    }))
  } else NULL

  list(
    per_seed_auc      = per_seed_df,
    summary           = summary_df,
    oof_predictions   = oof_predictions,
    pairwise_delong   = pairwise_df,
    pairwise_summary  = pairwise_summary
  )
}


# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.oms_auc_simple <- function(y, scores) {
  # Mann-Whitney U / nrow_pos / nrow_neg
  pos <- scores[y == 1L]
  neg <- scores[y == 0L]
  if (length(pos) == 0L || length(neg) == 0L) return(NA_real_)
  r <- rank(c(pos, neg))
  n_pos <- length(pos)
  n_neg <- length(neg)
  (sum(r[seq_len(n_pos)]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}


#' @keywords internal
#' @noRd
.grouped_stratified_folds <- function(y, groups, n_folds = 5L, seed = 7L) {
  # Build folds that (a) keep all rows of one group on the same side and
  # (b) approximately balance class counts per fold.
  set.seed(seed)
  groups <- as.character(groups)
  uniq <- unique(groups)
  # Class label per group: any() of the y values within that group (assumes
  # groups have homogeneous labels — typical for sample_id-grouped omics CV).
  grp_labels <- vapply(uniq, function(g) {
    yy <- y[groups == g]
    as.integer(any(yy == 1L, na.rm = TRUE))
  }, integer(1))
  pos_grps <- uniq[grp_labels == 1L]
  neg_grps <- uniq[grp_labels == 0L]
  pos_grps <- sample(pos_grps)
  neg_grps <- sample(neg_grps)
  pos_assign <- rep_len(seq_len(n_folds), length(pos_grps))
  neg_assign <- rep_len(seq_len(n_folds), length(neg_grps))
  fold_of_group <- setNames(integer(0), character(0))
  fold_of_group[pos_grps] <- pos_assign
  fold_of_group[neg_grps] <- neg_assign

  fold_of_row <- fold_of_group[groups]
  folds <- lapply(seq_len(n_folds), function(k) which(fold_of_row == k))
  folds
}
