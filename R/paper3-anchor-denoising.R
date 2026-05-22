# anchor_denoising_helpers.R
#
# Reference/anchor denoising primitives for Worker E3.
# All training-dependent objects are explicit fit objects. Scoring helpers do
# not silently fit on evaluation data because that would leak test information.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

.ad_clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}

.ad_normalise_mirna_names <- function(x) {
  out <- as.character(x)
  out <- sub("^hsa-mir-", "hsa-miR-", out, ignore.case = TRUE)
  out <- sub("^hsa-let-", "hsa-let-", out, ignore.case = TRUE)
  make.unique(out)
}

.ad_log_expr <- function(expr, pseudocount = NULL) {
  X <- as.matrix(expr)
  storage.mode(X) <- "double"
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("feature_", seq_len(ncol(X)))
  }
  colnames(X) <- .ad_normalise_mirna_names(colnames(X))
  if (any(X < 0, na.rm = TRUE)) {
    return(X)
  }
  if (is.null(pseudocount)) {
    pc <- 1e-6 * rowSums(X, na.rm = TRUE)
    pc[!is.finite(pc) | pc <= 0] <- .Machine$double.eps
  } else {
    pc <- rep(as.numeric(pseudocount), nrow(X))
    pc[!is.finite(pc) | pc <= 0] <- .Machine$double.eps
  }
  log(sweep(pmax(X, 0), 1L, pc, `+`))
}

.ad_infer_kit_col <- function(meta, kit_col = NULL) {
  if (!is.null(kit_col) && kit_col %in% names(meta)) return(kit_col)
  candidates <- c("kit_label", "library_kit", "extraction_kit_family",
                  "extraction_kit", "kit", "batch", "platform_id",
                  "modality")
  hit <- candidates[candidates %in% names(meta)][1L]
  if (length(hit) == 0L || is.na(hit)) NA_character_ else hit
}

.ad_kit_labels <- function(meta, kit_col = NULL) {
  meta <- as.data.frame(meta)
  kc <- .ad_infer_kit_col(meta, kit_col)
  if (is.na(kc)) {
    kit <- rep("unknown", nrow(meta))
  } else {
    kit <- .ad_clean_chr(meta[[kc]])
    kit[is.na(kit)] <- "unknown"
  }
  kit
}

.ad_weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

.ad_auc_fast <- function(y, score, sign_agnostic = TRUE) {
  y <- as.integer(y)
  score <- as.numeric(score)
  ok <- is.finite(score) & !is.na(y)
  y <- y[ok]
  score <- score[ok]
  if (length(unique(y)) != 2L || length(unique(score)) < 2L) return(NA_real_)
  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  if (n_pos < 1L || n_neg < 1L) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  auc <- (sum(ranks[y == 1L]) - n_pos * (n_pos + 1L) / 2) / (n_pos * n_neg)
  if (!is.finite(auc)) return(NA_real_)
  if (isTRUE(sign_agnostic)) max(auc, 1 - auc) else auc
}

.ad_hanley_var <- function(auc, n_pos, n_neg) {
  auc <- as.numeric(auc)
  if (!is.finite(auc) || n_pos < 1L || n_neg < 1L) return(NA_real_)
  q1 <- auc / (2 - auc)
  q2 <- 2 * auc^2 / (1 + auc)
  v <- (auc * (1 - auc) +
          (n_pos - 1) * (q1 - auc^2) +
          (n_neg - 1) * (q2 - auc^2)) / (n_pos * n_neg)
  if (!is.finite(v) || v < 0) NA_real_ else max(v, 1e-10)
}

.ad_hanley_ci <- function(auc, n_pos, n_neg) {
  v <- .ad_hanley_var(auc, n_pos, n_neg)
  if (!is.finite(v)) return(c(se = NA_real_, lo = NA_real_, hi = NA_real_))
  se <- sqrt(v)
  c(se = se, lo = max(0, auc - 1.96 * se), hi = min(1, auc + 1.96 * se))
}

.ad_find_feature <- function(feature_names, aliases) {
  idx <- which(tolower(feature_names) %in% tolower(aliases))
  if (length(idx) == 0L) NA_character_ else feature_names[idx[1L]]
}

.ad_hemo_index <- function(X_log) {
  hemo <- .ad_find_feature(colnames(X_log), c("hsa-miR-451a", "MIMAT0001631"))
  ref <- .ad_find_feature(colnames(X_log), c("hsa-miR-23a-3p", "MIMAT0000078"))
  if (is.na(hemo) || is.na(ref)) return(NULL)
  as.numeric(X_log[, hemo]) - as.numeric(X_log[, ref])
}

#' Fit kit-stable denominator anchors on training data only.
#'
#' Lower kit_variance_ratio means lower between-kit signal relative to
#' within-kit noise. The returned anchor table records the training cohorts and
#' denominator used by downstream scoring. Evaluation data must not be included.
#'
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param kit_col Optional kit column. When NULL, a kit-like column is inferred.
#' @param cohort_col Column identifying training cohorts.
#' @param anchor_fraction Fraction of eligible features selected as anchors.
#' @param min_anchors Minimum number of anchors to select.
#' @param min_samples_per_kit Minimum samples per kit stratum.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @export
fit_kit_stable_anchors <- function(expr, meta, kit_col = NULL,
                                   cohort_col = "cohort",
                                   anchor_fraction = 0.05,
                                   min_anchors = 8L,
                                   min_samples_per_kit = 5L,
                                   pseudocount = NULL) {
  X_log <- .ad_log_expr(expr, pseudocount = pseudocount)
  meta <- as.data.frame(meta)
  if (nrow(meta) != nrow(X_log)) {
    stop("fit_kit_stable_anchors: nrow(meta) must match nrow(expr).",
         call. = FALSE)
  }
  kit <- .ad_kit_labels(meta, kit_col)
  usable_kit <- !is.na(kit) & nzchar(kit)
  kit <- kit[usable_kit]
  X_log <- X_log[usable_kit, , drop = FALSE]
  kit_tab <- table(kit)
  keep_kits <- names(kit_tab)[kit_tab >= min_samples_per_kit]
  keep_rows <- kit %in% keep_kits
  kit <- kit[keep_rows]
  X_log <- X_log[keep_rows, , drop = FALSE]

  train_cohorts <- if (cohort_col %in% names(meta)) {
    sort(unique(.ad_clean_chr(meta[[cohort_col]])[usable_kit][keep_rows]))
  } else {
    "unknown_train_pool"
  }
  train_cohorts <- train_cohorts[!is.na(train_cohorts)]

  if (length(unique(kit)) < 2L) {
    empty <- data.frame(
      feature = character(),
      kit_variance_ratio = numeric(),
      selection_quantile = numeric(),
      train_cohorts_csv = character(),
      n_train_samples = integer(),
      stringsAsFactors = FALSE)
    fit <- list(
      anchors = character(),
      anchor_table = empty,
      kit_col = .ad_infer_kit_col(meta, kit_col),
      status = "not_evaluable_insufficient_training_kit_diversity",
      n_train_samples = nrow(X_log),
      train_cohorts = train_cohorts)
    class(fit) <- "kit_stable_anchor_fit"
    return(fit)
  }

  ratios <- vapply(seq_len(ncol(X_log)), function(j) {
    x <- X_log[, j]
    ok <- is.finite(x) & !is.na(kit)
    if (sum(ok) < min_samples_per_kit * 2L ||
        length(unique(kit[ok])) < 2L) return(NA_real_)
    g <- kit[ok]
    xv <- x[ok]
    tab <- table(g)
    if (sum(tab >= min_samples_per_kit) < 2L) return(NA_real_)
    keep <- g %in% names(tab)[tab >= min_samples_per_kit]
    g <- g[keep]
    xv <- xv[keep]
    group_means <- tapply(xv, g, mean, na.rm = TRUE)
    group_vars <- tapply(xv, g, function(v) {
      if (length(v) < 2L) NA_real_ else var(v, na.rm = TRUE)
    })
    group_n <- tapply(xv, g, length)
    if (length(group_means) < 2L) return(NA_real_)
    global_mean <- mean(xv, na.rm = TRUE)
    between <- sum(group_n * (group_means - global_mean)^2) /
      max(1L, length(group_means) - 1L)
    within <- weighted.mean(group_vars, group_n, na.rm = TRUE)
    between / (within + .Machine$double.eps)
  }, numeric(1L))
  names(ratios) <- colnames(X_log)

  eligible <- names(ratios)[is.finite(ratios)]
  if (length(eligible) == 0L) {
    empty <- data.frame(
      feature = character(),
      kit_variance_ratio = numeric(),
      selection_quantile = numeric(),
      train_cohorts_csv = character(),
      n_train_samples = integer(),
      stringsAsFactors = FALSE)
    fit <- list(
      anchors = character(),
      anchor_table = empty,
      kit_col = .ad_infer_kit_col(meta, kit_col),
      status = "not_evaluable_no_finite_kit_variance_ratios",
      n_train_samples = nrow(X_log),
      train_cohorts = train_cohorts)
    class(fit) <- "kit_stable_anchor_fit"
    return(fit)
  }

  n_select <- max(as.integer(min_anchors),
                  ceiling(anchor_fraction * length(eligible)))
  n_select <- min(n_select, length(eligible))
  ord <- eligible[order(ratios[eligible], decreasing = FALSE)]
  anchors <- ord[seq_len(n_select)]
  rank_all <- rank(ratios[eligible], ties.method = "first")
  names(rank_all) <- eligible
  anchor_table <- data.frame(
    feature = anchors,
    kit_variance_ratio = as.numeric(ratios[anchors]),
    selection_quantile = as.numeric(rank_all[anchors]) / length(eligible),
    train_cohorts_csv = paste(train_cohorts, collapse = ","),
    n_train_samples = nrow(X_log),
    stringsAsFactors = FALSE)

  fit <- list(
    anchors = anchors,
    anchor_table = anchor_table,
    ratios = ratios,
    kit_col = .ad_infer_kit_col(meta, kit_col),
    status = "ok",
    anchor_fraction = anchor_fraction,
    min_anchors = min_anchors,
    n_eligible_features = length(eligible),
    n_train_samples = nrow(X_log),
    train_cohorts = train_cohorts)
  class(fit) <- "kit_stable_anchor_fit"
  fit
}

.ad_require_anchor_fit <- function(anchor_fit, train_expr, train_meta, ...) {
  if (!is.null(anchor_fit)) return(anchor_fit)
  if (!is.null(train_expr) && !is.null(train_meta)) {
    return(fit_kit_stable_anchors(train_expr, train_meta, ...))
  }
  stop("A train-only anchor_fit is required. Pass anchor_fit, or pass ",
       "train_expr/train_meta explicitly. Refusing to fit on evaluation data.",
       call. = FALSE)
}

#' Score a panel against train-only kit-stable anchors.
#'
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param panel Character vector of panel features.
#' @param anchor_fit Optional fit from `fit_kit_stable_anchors()`.
#' @param train_expr Optional training matrix used to fit anchors when
#'   `anchor_fit` is NULL.
#' @param train_meta Optional training metadata used with `train_expr`.
#' @param kit_col Optional kit column.
#' @param min_eval_anchors Minimum anchors required in evaluation data.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @param ... Additional arguments passed to `fit_kit_stable_anchors()`.
#' @export
score_kit_stable_anchor_rclr <- function(expr, meta, panel,
                                         anchor_fit = NULL,
                                         train_expr = NULL,
                                         train_meta = NULL,
                                         kit_col = NULL,
                                         min_eval_anchors = 3L,
                                         pseudocount = NULL,
                                         ...) {
  fit <- .ad_require_anchor_fit(anchor_fit, train_expr, train_meta,
                                kit_col = kit_col, pseudocount = pseudocount,
                                ...)
  X_log <- .ad_log_expr(expr, pseudocount = pseudocount)
  panel <- intersect(panel, colnames(X_log))
  anchors <- intersect(fit$anchors, colnames(X_log))
  if (length(panel) == 0L || length(anchors) < min_eval_anchors) {
    out <- rep(NA_real_, nrow(X_log))
    attr(out, "status") <- "not_evaluable_insufficient_panel_or_anchor_overlap"
    attr(out, "n_anchors_eval") <- length(anchors)
    return(out)
  }
  denom <- rowMeans(X_log[, anchors, drop = FALSE], na.rm = TRUE)
  score <- rowSums(sweep(X_log[, panel, drop = FALSE], 1L, denom, "-"),
                   na.rm = TRUE)
  attr(score, "status") <- fit$status
  attr(score, "n_anchors_eval") <- length(anchors)
  attr(score, "anchors_eval") <- anchors
  score
}

#' Fit the hemolysis plus kit-stable dual-anchor denominator.
#'
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param kit_col Optional kit column.
#' @param cohort_col Column identifying training cohorts.
#' @param anchor_fit Optional precomputed `kit_stable_anchor_fit`.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @param ... Additional arguments passed to `fit_kit_stable_anchors()`.
#' @export
fit_hemolysis_kit_dual_anchor <- function(expr, meta, kit_col = NULL,
                                          cohort_col = "cohort",
                                          anchor_fit = NULL,
                                          pseudocount = NULL,
                                          ...) {
  fit <- anchor_fit %||% fit_kit_stable_anchors(
    expr, meta, kit_col = kit_col, cohort_col = cohort_col,
    pseudocount = pseudocount, ...)
  X_log <- .ad_log_expr(expr, pseudocount = pseudocount)
  anchors <- intersect(fit$anchors, colnames(X_log))
  hemo <- .ad_hemo_index(X_log)
  hemo_ok <- !is.null(hemo) && length(anchors) >= 3L &&
    sum(is.finite(hemo)) >= 10L
  slope <- 0
  hemo_center <- NA_real_
  if (hemo_ok) {
    anchor_center <- rowMeans(X_log[, anchors, drop = FALSE], na.rm = TRUE)
    hemo_center <- median(hemo[is.finite(hemo)], na.rm = TRUE)
    df <- data.frame(anchor_center = anchor_center,
                     hemo_centered = hemo - hemo_center)
    ok <- is.finite(df$anchor_center) & is.finite(df$hemo_centered)
    if (sum(ok) >= 10L && requireNamespace("MASS", quietly = TRUE)) {
      mod <- tryCatch(MASS::rlm(anchor_center ~ hemo_centered,
                                data = df[ok, , drop = FALSE],
                                method = "M",
                                psi = MASS::psi.huber),
                      error = function(e) NULL,
                      warning = function(w) NULL)
      if (!is.null(mod) && "hemo_centered" %in% names(coef(mod))) {
        slope <- as.numeric(coef(mod)[["hemo_centered"]])
      }
    } else if (sum(ok) >= 10L) {
      mod <- tryCatch(stats::lm(anchor_center ~ hemo_centered,
                                data = df[ok, , drop = FALSE]),
                      error = function(e) NULL)
      if (!is.null(mod) && "hemo_centered" %in% names(coef(mod))) {
        slope <- as.numeric(coef(mod)[["hemo_centered"]])
      }
    }
  }
  out <- list(
    anchor_fit = fit,
    hemo_ok = isTRUE(hemo_ok),
    hemo_slope = ifelse(is.finite(slope), slope, 0),
    hemo_center = hemo_center,
    status = if (isTRUE(hemo_ok)) "ok" else "anchor_only_hemo_proxy_missing_or_weak",
    n_train_samples = nrow(X_log))
  class(out) <- "hemolysis_kit_dual_anchor_fit"
  out
}

#' Score a panel against the dual hemolysis plus kit-stable denominator.
#'
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param panel Character vector of panel features.
#' @param dual_fit Optional fit from `fit_hemolysis_kit_dual_anchor()`.
#' @param train_expr Optional training matrix used to fit `dual_fit`.
#' @param train_meta Optional training metadata used with `train_expr`.
#' @param kit_col Optional kit column.
#' @param min_eval_anchors Minimum anchors required in evaluation data.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @param ... Additional arguments passed to `fit_hemolysis_kit_dual_anchor()`.
#' @export
score_hemolysis_kit_dual_anchor <- function(expr, meta, panel,
                                            dual_fit = NULL,
                                            train_expr = NULL,
                                            train_meta = NULL,
                                            kit_col = NULL,
                                            min_eval_anchors = 3L,
                                            pseudocount = NULL,
                                            ...) {
  if (is.null(dual_fit)) {
    if (is.null(train_expr) || is.null(train_meta)) {
      stop("A train-only dual_fit is required. Pass dual_fit, or pass ",
           "train_expr/train_meta explicitly. Refusing to fit on evaluation data.",
           call. = FALSE)
    }
    dual_fit <- fit_hemolysis_kit_dual_anchor(
      train_expr, train_meta, kit_col = kit_col, pseudocount = pseudocount, ...)
  }
  X_log <- .ad_log_expr(expr, pseudocount = pseudocount)
  panel <- intersect(panel, colnames(X_log))
  anchors <- intersect(dual_fit$anchor_fit$anchors, colnames(X_log))
  if (length(panel) == 0L || length(anchors) < min_eval_anchors) {
    out <- rep(NA_real_, nrow(X_log))
    attr(out, "status") <- "not_evaluable_insufficient_panel_or_anchor_overlap"
    attr(out, "n_anchors_eval") <- length(anchors)
    return(out)
  }
  denom <- rowMeans(X_log[, anchors, drop = FALSE], na.rm = TRUE)
  hemo <- .ad_hemo_index(X_log)
  if (isTRUE(dual_fit$hemo_ok) && !is.null(hemo) &&
      is.finite(dual_fit$hemo_center)) {
    denom <- denom + dual_fit$hemo_slope * (hemo - dual_fit$hemo_center)
  }
  score <- rowSums(sweep(X_log[, panel, drop = FALSE], 1L, denom, "-"),
                   na.rm = TRUE)
  attr(score, "status") <- dual_fit$status
  attr(score, "n_anchors_eval") <- length(anchors)
  attr(score, "hemo_ok") <- isTRUE(dual_fit$hemo_ok)
  score
}

.ad_rin_values <- function(meta, rin_col = NULL) {
  meta <- as.data.frame(meta)
  if (!is.null(rin_col) && rin_col %in% names(meta)) {
    raw <- meta[[rin_col]]
  } else {
    candidates <- grep("rin|rna_integrity|rqi", names(meta),
                       ignore.case = TRUE, value = TRUE)
    raw <- if (length(candidates) > 0L) meta[[candidates[1L]]] else rep(NA, nrow(meta))
  }
  raw <- suppressWarnings(as.numeric(gsub("[^0-9.+-]", "", as.character(raw))))
  raw[!is.finite(raw) | raw <= 0] <- NA_real_
  raw
}

#' Fit a RIN-weighted frozen reference profile on training samples only.
#'
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param rin_col Optional RNA-integrity-number column.
#' @param cohort_col Column identifying training cohorts.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @export
fit_rin_weighted_reference <- function(expr, meta, rin_col = NULL,
                                       cohort_col = "cohort",
                                       pseudocount = NULL) {
  X_log <- .ad_log_expr(expr, pseudocount = pseudocount)
  meta <- as.data.frame(meta)
  if (nrow(meta) != nrow(X_log)) {
    stop("fit_rin_weighted_reference: nrow(meta) must match nrow(expr).",
         call. = FALSE)
  }
  rin <- .ad_rin_values(meta, rin_col = rin_col)
  fallback <- median(rin[is.finite(rin)], na.rm = TRUE)
  if (!is.finite(fallback)) fallback <- 1
  w <- rin
  w[!is.finite(w)] <- fallback
  w <- w / median(w[is.finite(w) & w > 0], na.rm = TRUE)
  w[!is.finite(w) | w <= 0] <- 1
  ref <- vapply(seq_len(ncol(X_log)), function(j) {
    .ad_weighted_mean(X_log[, j], w)
  }, numeric(1L))
  names(ref) <- colnames(X_log)
  train_cohorts <- if (cohort_col %in% names(meta)) {
    sort(unique(.ad_clean_chr(meta[[cohort_col]])))
  } else {
    "unknown_train_pool"
  }
  out <- list(
    reference = ref,
    rin_col = rin_col %||% grep("rin|rna_integrity|rqi", names(meta),
                                ignore.case = TRUE, value = TRUE)[1L] %||% NA_character_,
    median_weight = fallback,
    rin_missing_fraction = mean(!is.finite(rin)),
    n_train_samples = nrow(X_log),
    train_cohorts = train_cohorts,
    status = if (all(!is.finite(rin))) "ok_all_rin_missing_median_weight_fallback" else "ok")
  class(out) <- "rin_weighted_reference_fit"
  out
}

#' Score a panel against a RIN-weighted frozen training reference profile.
#'
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param panel Character vector of panel features.
#' @param rin_fit Optional fit from `fit_rin_weighted_reference()`.
#' @param train_expr Optional training matrix used to fit `rin_fit`.
#' @param train_meta Optional training metadata used with `train_expr`.
#' @param rin_col Optional RNA-integrity-number column.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @export
score_rin_weighted_score <- function(expr, meta, panel,
                                     rin_fit = NULL,
                                     train_expr = NULL,
                                     train_meta = NULL,
                                     rin_col = NULL,
                                     pseudocount = NULL) {
  if (is.null(rin_fit)) {
    if (is.null(train_expr) || is.null(train_meta)) {
      stop("A train-only rin_fit is required. Pass rin_fit, or pass ",
           "train_expr/train_meta explicitly. Refusing to fit on evaluation data.",
           call. = FALSE)
    }
    rin_fit <- fit_rin_weighted_reference(train_expr, train_meta,
                                          rin_col = rin_col,
                                          pseudocount = pseudocount)
  }
  X_log <- .ad_log_expr(expr, pseudocount = pseudocount)
  panel <- intersect(panel, intersect(colnames(X_log), names(rin_fit$reference)))
  if (length(panel) == 0L) {
    out <- rep(NA_real_, nrow(X_log))
    attr(out, "status") <- "not_evaluable_no_panel_reference_overlap"
    return(out)
  }
  score <- rowSums(sweep(X_log[, panel, drop = FALSE], 2L,
                         rin_fit$reference[panel], "-"),
                   na.rm = TRUE)
  attr(score, "status") <- rin_fit$status
  attr(score, "rin_missing_fraction") <- rin_fit$rin_missing_fraction
  score
}

#' Optional ComBat-seq followed by rCLR scoring.
#'
#' @param expr Numeric count matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param panel Character vector of panel features.
#' @param batch_col Optional batch column. Defaults to inferred kit column.
#' @param disease_col Optional disease-label column passed to ComBat-seq.
#' @param allow_transductive Logical. Must be TRUE because this comparator is
#'   transductive and should not be treated as a frozen deployable scorer.
#' @param pseudocount Additive pseudocount before rCLR scoring.
#' @export
score_combat_seq_then_rclr <- function(expr, meta, panel,
                                       batch_col = NULL,
                                       disease_col = NULL,
                                       allow_transductive = FALSE,
                                       pseudocount = 0.5) {
  if (!isTRUE(allow_transductive)) {
    stop("score_combat_seq_then_rclr requires allow_transductive=TRUE. ",
         "ComBat-seq has no frozen train/apply interface here, so using all ",
         "evaluation samples is transductive and must not be mixed with ",
         "train-only LOCO claims.", call. = FALSE)
  }
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("score_combat_seq_then_rclr requires the sva package.", call. = FALSE)
  }
  X <- as.matrix(expr)
  X[!is.finite(X) | X < 0] <- 0
  counts <- round(X)
  meta <- as.data.frame(meta)
  batch_col <- batch_col %||% .ad_infer_kit_col(meta)
  if (is.na(batch_col) || !batch_col %in% names(meta)) {
    stop("score_combat_seq_then_rclr: no usable batch/kit column in meta.",
         call. = FALSE)
  }
  batch <- .ad_clean_chr(meta[[batch_col]])
  batch[is.na(batch)] <- "unknown"
  group <- NULL
  if (!is.null(disease_col) && disease_col %in% names(meta)) {
    group <- as.factor(meta[[disease_col]])
  }
  corrected <- sva::ComBat_seq(counts = t(counts), batch = batch, group = group)
  corrected <- t(corrected)
  colnames(corrected) <- colnames(expr)
  rownames(corrected) <- rownames(expr)
  score_rclr_baseline(corrected, meta, panel, pseudocount = pseudocount)
}

#' Baseline rCLR panel sum using the standard trimmed denominator.
#'
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param panel Character vector of panel features.
#' @param trim_upper,trim_lower Feature trimming fractions for rCLR denominator
#'   selection.
#' @param exclude_features Features excluded from denominator selection.
#' @param min_centering_size Minimum denominator size before fallback.
#' @param pseudocount Optional additive pseudocount before log transform.
score_rclr_baseline <- function(expr, meta, panel,
                                trim_upper = 0.10,
                                trim_lower = 0.05,
                                exclude_features = c("hsa-miR-451a",
                                                     "hsa-miR-16-5p",
                                                     "hsa-miR-486-5p",
                                                     "hsa-miR-144-3p",
                                                     "hsa-miR-223-3p"),
                                min_centering_size = 8L,
                                pseudocount = NULL) {
  X_log <- .ad_log_expr(expr, pseudocount = pseudocount)
  panel <- intersect(panel, colnames(X_log))
  if (length(panel) == 0L) return(rep(NA_real_, nrow(X_log)))
  features <- colnames(X_log)
  is_excluded <- features %in% exclude_features
  out <- matrix(NA_real_, nrow = nrow(X_log), ncol = length(panel))
  colnames(out) <- panel
  for (i in seq_len(nrow(X_log))) {
    x <- X_log[i, ]
    finite <- is.finite(x)
    candidates <- finite & !is_excluded
    n_cand <- sum(candidates)
    if (n_cand < min_centering_size) {
      center <- mean(x[finite], na.rm = TRUE)
    } else {
      ranks <- rank(-x, ties.method = "first")
      upper_cut <- floor(trim_upper * n_cand)
      lower_cut <- floor(trim_lower * n_cand)
      cand_ranks <- sort(ranks[candidates])
      keep_lo <- cand_ranks[upper_cut + 1L]
      keep_hi <- cand_ranks[n_cand - lower_cut]
      keep <- candidates & ranks >= keep_lo & ranks <= keep_hi
      center <- if (sum(keep) >= min_centering_size) {
        mean(x[keep], na.rm = TRUE)
      } else {
        mean(x[finite], na.rm = TRUE)
      }
    }
    out[i, ] <- x[panel] - center
  }
  rowSums(out, na.rm = TRUE)
}

cluster_bootstrap_auc <- function(y, score, cluster_id = NULL, n_boot = 1000L,
                                  seed = 42L) {
  y <- as.integer(y)
  score <- as.numeric(score)
  n <- length(y)
  if (is.null(cluster_id) || length(cluster_id) != n) cluster_id <- seq_len(n)
  point <- .ad_auc_fast(y, score)
  if (!is.finite(point)) {
    return(list(point_auc = NA_real_, ci_lower = NA_real_,
                ci_upper = NA_real_, n_boot_effective = 0L))
  }
  rows_by_cluster <- split(seq_len(n), factor(cluster_id))
  set.seed(seed)
  draws <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    picked <- sample(seq_along(rows_by_cluster), length(rows_by_cluster),
                     replace = TRUE)
    rows <- unlist(rows_by_cluster[picked], use.names = FALSE)
    draws[b] <- .ad_auc_fast(y[rows], score[rows])
  }
  draws <- draws[is.finite(draws)]
  list(point_auc = point,
       ci_lower = if (length(draws)) unname(stats::quantile(draws, 0.025)) else NA_real_,
       ci_upper = if (length(draws)) unname(stats::quantile(draws, 0.975)) else NA_real_,
       n_boot_effective = length(draws))
}

cluster_bootstrap_delta_auc <- function(y, score_method, score_baseline,
                                        cluster_id = NULL, n_boot = 1000L,
                                        seed = 42L) {
  y <- as.integer(y)
  n <- length(y)
  if (is.null(cluster_id) || length(cluster_id) != n) cluster_id <- seq_len(n)
  point <- .ad_auc_fast(y, score_method) - .ad_auc_fast(y, score_baseline)
  rows_by_cluster <- split(seq_len(n), factor(cluster_id))
  set.seed(seed)
  draws <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    picked <- sample(seq_along(rows_by_cluster), length(rows_by_cluster),
                     replace = TRUE)
    rows <- unlist(rows_by_cluster[picked], use.names = FALSE)
    draws[b] <- .ad_auc_fast(y[rows], score_method[rows]) -
      .ad_auc_fast(y[rows], score_baseline[rows])
  }
  draws <- draws[is.finite(draws)]
  list(point_delta = point,
       ci_lower = if (length(draws)) unname(stats::quantile(draws, 0.025)) else NA_real_,
       ci_upper = if (length(draws)) unname(stats::quantile(draws, 0.975)) else NA_real_,
       n_boot_effective = length(draws))
}

random_effects_pool <- function(theta, var, labels = NULL) {
  ok <- is.finite(theta) & is.finite(var) & var > 0
  theta <- theta[ok]
  var <- var[ok]
  if (!length(theta)) {
    return(list(estimate = NA_real_, se = NA_real_, lo = NA_real_,
                hi = NA_real_, tau2 = NA_real_, I2 = NA_real_, n = 0L))
  }
  w <- 1 / var
  fixed <- sum(w * theta) / sum(w)
  q <- sum(w * (theta - fixed)^2)
  df <- length(theta) - 1L
  c_val <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (length(theta) > 1L && c_val > 0) max(0, (q - df) / c_val) else 0
  w_re <- 1 / (var + tau2)
  est <- sum(w_re * theta) / sum(w_re)
  se <- sqrt(1 / sum(w_re))
  I2 <- if (length(theta) > 1L && q > 0) 100 * max(0, (q - df) / q) else 0
  list(estimate = est, se = se, lo = est - 1.96 * se,
       hi = est + 1.96 * se, tau2 = tau2, I2 = I2, n = length(theta))
}

write_kit_stable_anchors <- function(anchor_fit, path) {
  tab <- anchor_fit$anchor_table
  required <- c("feature", "kit_variance_ratio", "selection_quantile",
                "train_cohorts_csv", "n_train_samples")
  if (is.null(tab) || nrow(tab) == 0L) {
    tab <- data.frame(feature = character(),
                      kit_variance_ratio = numeric(),
                      selection_quantile = numeric(),
                      train_cohorts_csv = character(),
                      n_train_samples = integer(),
                      stringsAsFactors = FALSE)
  }
  tab <- tab[, required, drop = FALSE]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(tab)
}
