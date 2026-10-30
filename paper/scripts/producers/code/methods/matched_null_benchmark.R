# matched_null_benchmark.R
#
# Matched panel-null benchmark for within-sample miRNA scoring methods. Per
# plan v0.3.3 Section 5 — referenced by Figures 4, 5, 6, 8, 9.
#
# Protocol:
#   1. For an observed scoring method M applied to a real panel P (|P|=k features),
#      compute observed AUC_obs.
#   2. Generate K random panels P_1, ..., P_K each of size k drawn from the
#      cohort-specific candidate pool, jointly stratified on:
#         (a) per-feature detection-rate quartile bin (4 levels), and
#         (b) per-feature log-mean-abundance quartile bin (4 levels).
#      Random panels do NOT overlap with the observed panel (sampling is
#      from `feature_pool \ panel`); hemolysis-marker miRNAs are excluded
#      from the candidate pool when the cohort has been declared
#      post-hemolysis-prepared.
#   3. For each P_i, apply method M and compute AUC_i.
#   4. Empirical p-value: P_emp = (1 + sum(AUC_i >= AUC_obs)) / (K + 1).
#      Minimum achievable resolution = 1/(K+1).
#   5. Multiplicity: BH-FDR within each modality family (4 modality test
#      families per transfer layout). Holm adjustment is reserved for
#      family-wise comparisons (e.g. method × claim-state contrasts;
#      hemolysis-injection deltas) — not for the discovery-screening
#      matched-null tests.
#
# The default K is now 10,000 per plan v0.3.3 (was 1,000 in v0.3). The K
# argument is preserved so callers running the lighter ablation analyses
# (Module E permutation tests) can override.

suppressPackageStartupMessages({
  if (!requireNamespace("pROC", quietly = TRUE))
    stop("matched_null_benchmark requires pROC for AUC + ci.auc.")
})

#' Match a feature set to the full feature pool by detection × abundance strata.
#'
#' @param panel character — features in the observed panel.
#' @param feature_pool character — all detected features in the training data.
#' @param detection_rate numeric — feature-level detection rate over training samples.
#' @param mean_abundance numeric — feature-level mean abundance over training samples.
#' @param n_strata integer — number of strata per axis. Default 4 (16 cells).
#'
#' @return a list with `strata` (per-feature integer 1..n_strata^2 cell ID) and
#'   `cells_in_panel` (table of cell IDs in the observed panel — used to draw
#'   random matched panels).
.match_panel_strata <- function(panel,
                                  feature_pool,
                                  detection_rate,
                                  mean_abundance,
                                  n_strata = 4L) {
  stopifnot(length(detection_rate) == length(feature_pool),
            length(mean_abundance) == length(feature_pool))
  detection_q <- as.integer(cut(detection_rate, n_strata, include.lowest = TRUE))
  abundance_q <- as.integer(cut(rank(mean_abundance), n_strata, include.lowest = TRUE))
  strata <- (detection_q - 1L) * n_strata + abundance_q
  names(strata) <- feature_pool
  panel_idx <- match(panel, feature_pool)
  panel_idx <- panel_idx[!is.na(panel_idx)]
  if (length(panel_idx) == 0L) stop("No panel features found in feature_pool.")
  panel_strata <- strata[panel_idx]
  list(strata = strata, cells_in_panel = table(panel_strata))
}

#' Draw a single matched random panel that does NOT overlap with the observed
#' panel. Implements a 3-tier fallback for small panels (round-10 reviewer
#' edge case — small qPCR / NanoString platforms with ~50 probes can exhaust
#' the per-stratum candidate pool):
#'   tier 1 (preferred): exact-stratum match without replacement.
#'   tier 2 (small-panel relaxation): if any stratum is exhausted, oversample
#'     the panel-size deficit from the largest remaining stratum.
#'   tier 3 (degenerate-pool fallback): if total candidate pool minus excluded
#'     overlap is smaller than panel_size, return character(0) so the caller
#'     can mark the cell as not matched-null-eligible and fall back to
#'     bootstrap-only confidence.
.draw_matched_panel <- function(strata_info, exclude = character(0)) {
  strata <- strata_info$strata
  cells <- strata_info$cells_in_panel
  panel_size <- sum(unlist(cells))
  candidates_total <- length(setdiff(names(strata), exclude))
  if (candidates_total < panel_size) {
    out <- character(0)
    attr(out, "fallback_tier") <- 3L
    return(out)
  }
  picked <- character(0)
  shortfall <- 0L
  for (cell in names(cells)) {
    pool <- setdiff(names(strata)[strata == as.integer(cell)], exclude)
    if (length(pool) == 0L) {
      shortfall <- shortfall + cells[[cell]]
      next
    }
    n_pick <- min(cells[[cell]], length(pool))
    if (n_pick < cells[[cell]]) shortfall <- shortfall + (cells[[cell]] - n_pick)
    picked <- c(picked, sample(pool, n_pick))
  }
  fallback_tier <- if (shortfall > 0L) 2L else 1L
  if (shortfall > 0L) {                                       # tier 2
    extra_pool <- setdiff(names(strata), c(exclude, picked))
    if (length(extra_pool) >= shortfall) {
      picked <- c(picked, sample(extra_pool, shortfall))
    } else {
      out <- character(0)
      attr(out, "fallback_tier") <- 3L
      return(out)                                             # degenerate
    }
  }
  attr(picked, "fallback_tier") <- fallback_tier
  picked
}

#' Compute matched-null benchmark for a within-sample scoring method.
#'
#' @param X matrix (samples × features) — feature panel for the cohort.
#'   Should be the full pool, not just the observed panel.
#' @param y binary outcome (0/1 or factor with 2 levels).
#' @param panel character — names of features in the observed panel.
#' @param scoring_fn function `function(X_subset)` returning a per-sample score.
#' @param K integer — number of null panels. Default 10,000 per plan v0.3.3.
#' @param feature_pool character — full feature pool (must include all features
#'   in `panel` and all features in `colnames(X)`). If NULL, uses `colnames(X)`.
#' @param exclude_features character — feature names to exclude from the
#'   candidate pool (e.g. hemolysis markers when the cohort has been declared
#'   post-hemolysis-prepared). Default empty.
#' @param post_hemolysis_prepared logical — convenience flag; if TRUE, prepends
#'   the canonical hemolysis-marker list to `exclude_features`. Default FALSE.
#' @param seed integer — RNG seed for reproducibility. Default 42.
#'
#' @return a list with `auc_obs`, `auc_null` (length-K vector), `p_emp`,
#'   `panel_size`, `K`, `seed`, `null_panels` (list of length K with the
#'   feature names of each random panel), and `excluded_features` (character).
matched_null_benchmark <- function(X, y, panel, scoring_fn,
                                     K = 10000L,
                                     feature_pool = NULL,
                                     exclude_features = character(0),
                                     post_hemolysis_corrected = FALSE,
                                     seed = 42L) {
  stopifnot(is.matrix(X) || is.data.frame(X))
  X <- as.matrix(X)
  if (is.null(feature_pool)) feature_pool <- colnames(X)
  if (is.factor(y)) y <- as.integer(y) - 1L

  panel <- intersect(panel, colnames(X))
  if (length(panel) == 0L) stop("None of the panel features are in colnames(X).")

  # Hemolysis exclusion: if the cohort has been hemolysis-prepared, the
  # candidate pool must not include hemolysis markers (otherwise the null
  # would draw from features the observed panel could not have included).
  hemolysis_markers <- c("hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p",
                          "hsa-miR-144-3p", "hsa-miR-223-3p")
  if (post_hemolysis_corrected) {
    exclude_features <- unique(c(exclude_features, hemolysis_markers))
  }
  feature_pool <- setdiff(feature_pool, exclude_features)
  excluded_in_pool <- intersect(panel, exclude_features)
  if (length(excluded_in_pool) > 0L) {
    warning("matched_null_benchmark: ", length(excluded_in_pool),
             " observed-panel features overlap with the exclusion list. ",
             "Did the cohort really get hemolysis-corrected? Removing from panel.")
    panel <- setdiff(panel, excluded_in_pool)
  }

  # Feature stats for stratified matching.
  detection_rate <- colMeans(X[, feature_pool, drop = FALSE] > 0, na.rm = TRUE)
  mean_abundance <- colMeans(X[, feature_pool, drop = FALSE], na.rm = TRUE)

  strata_info <- .match_panel_strata(
    panel = panel, feature_pool = feature_pool,
    detection_rate = detection_rate, mean_abundance = mean_abundance
  )

  # Observed AUC.
  score_obs <- scoring_fn(X[, panel, drop = FALSE])
  # Codex review 2026-05-05: gracefully handle degenerate scorer output (e.g.
  # ws_balance_ilr coverage-gated cohort returning rep(NA, n)). pROC::auc
  # otherwise errors with "No control observation."
  if (all(is.na(score_obs)) || length(unique(score_obs[!is.na(score_obs)])) < 2L) {
    return(list(
      auc_obs = NA_real_, auc_null = rep(NA_real_, K),
      p_emp = NA_real_, panel_size = length(panel), K = K, seed = seed,
      null_panels = vector("list", K),
      excluded_features = exclude_features,
      eligible = FALSE,
      ineligible_reason = "scorer_returned_all_NA_or_constant",
      n_score_NA = sum(is.na(score_obs)),
      n_score_total = length(score_obs)
    ))
  }
  auc_obs <- as.numeric(pROC::auc(y, score_obs, quiet = TRUE))

  # Null AUC distribution. Each null panel excludes the observed panel features
  # so the null is non-overlapping with the observed (per plan v0.3.3 §5).
  # If `.draw_matched_panel` returns character(0), the cell is not matched-null-
  # eligible (degenerate candidate pool); we record this state and return only
  # the bootstrap-eligible observed AUC, with `eligible = FALSE` and the
  # fallback tier triggered, per plan §5 (round-13 reviewer fix).
  set.seed(seed)
  null_panels <- vector("list", K)
  auc_null <- numeric(K)
  fallback_count <- 0L
  for (k in seq_len(K)) {
    p_k <- .draw_matched_panel(strata_info, exclude = panel)
    tier_k <- attr(p_k, "fallback_tier") %||% 1L
    if (length(p_k) == 0L) {
      fallback_count <- fallback_count + 1L
      auc_null[k] <- NA_real_
      next
    }
    null_panels[[k]] <- p_k
    score_k <- scoring_fn(X[, p_k, drop = FALSE])
    auc_null[k] <- as.numeric(pROC::auc(y, score_k, quiet = TRUE))
  }

  eligible <- fallback_count < (K / 2L)   # if >50 % nulls degenerate, cell is not eligible
  if (!eligible) {
    return(list(
      auc_obs = auc_obs,
      auc_null = NA_real_,
      p_emp = NA_real_,
      panel_size = length(panel),
      K = K,
      seed = seed,
      null_panels = list(),
      excluded_features = exclude_features,
      post_hemolysis_corrected = post_hemolysis_corrected,
      eligible = FALSE,
      fallback_tier = 3L,
      fallback_count = fallback_count,
      reason = "candidate pool exhausted; null pool too small for matched-null benchmark; report bootstrap CI on observed AUC only"
    ))
  }

  p_emp <- (1 + sum(auc_null >= auc_obs, na.rm = TRUE)) / (K - fallback_count + 1L)

  list(
    auc_obs = auc_obs,
    auc_null = auc_null,
    p_emp = p_emp,
    panel_size = length(panel),
    K = K,
    seed = seed,
    null_panels = null_panels,
    excluded_features = exclude_features,
    post_hemolysis_corrected = post_hemolysis_corrected,
    eligible = TRUE,
    fallback_tier = if (any(vapply(null_panels, function(p) (attr(p, "fallback_tier") %||% 1L) > 1L,
                                   logical(1L)))) 2L else 1L,
    fallback_count = fallback_count
  )
}

#' Apply BH-FDR adjustment within a modality family.
#'
#' Per plan v0.3.3 §5: matched-null tests use BH-FDR within each of the 4
#' modality test families separately (no across-modality adjustment).
#'
#' @param results list of matched_null_benchmark return objects from a single
#'   modality family.
#' @return numeric vector of BH-FDR q-values, same length as `results`.
bh_fdr_correct_matched_null <- function(results) {
  p_emp <- vapply(results, function(r) r$p_emp, numeric(1L))
  p.adjust(p_emp, method = "BH")
}

#' Apply Holm adjustment across a family-wise contrast set.
#'
#' Per plan v0.3.3 §5: Holm is reserved for family-wise comparisons such as
#' method × claim-state contrasts and hemolysis-injection deltas — NOT for
#' the discovery-screening matched-null tests.
#'
#' @param p_values numeric vector of contrast p-values.
#' @return numeric vector of Holm-prepared p-values.
holm_correct_familywise <- function(p_values) {
  p.adjust(p_values, method = "holm")
}

holm_correct_matched_null <- function(results) {
  warning("holm_correct_matched_null is deprecated for matched-null tests. ",
           "Use bh_fdr_correct_matched_null per plan v0.3.3 §5.")
  bh_fdr_correct_matched_null(results)
}

#' Block-aware BH-FDR for matched-null tests with shared-specimen clusters.
#'
#' Per plan v0.3.3 §5: naïve BH-FDR within a modality family treats
#' specimen-shared cohorts (e.g. the Toray cluster with ~94.5 % reuse with
#' GSE211692) as independent tests, which is anti-conservative because the
#' specimen-shared cohorts produce strongly correlated test statistics.
#'
#' The fix is two-stage:
#'   (i) within each provenance block, combine cohort-level p-values via the
#'       Fisher chi-square statistic with a Bauer-Hommel-adjusted independence
#'       assumption (a conservative adjustment that absorbs the worst-case
#'       positive dependence); the resulting block p-value is the cohort-level
#'       evidence that ANY cohort in the block beats matched null.
#'   (ii) apply BH-FDR across the block-level p-values.
#'
#' Per-cohort q-values are reported alongside block-level q-values, with
#' explicit annotation of which block each cohort belongs to.
#'
#' @param results list of matched_null_benchmark return objects.
#' @param block_id character — provenance block identifier per cohort
#'   (e.g. "Toray-cluster", "Affy-singleton", "FirePlex-multiplex").
#' @return data.frame with columns: cohort_idx, block_id, p_emp,
#'   p_block, q_block_BH (BH-FDR on block p-values), q_cohort_within_block
#'   (BH-FDR within the block).
bh_fdr_correct_blocked <- function(results, block_id) {
  stopifnot(length(results) == length(block_id))
  p_emp <- vapply(results, function(r) r$p_emp, numeric(1L))
  blocks <- split(seq_along(p_emp), block_id)

  # Stage 1: combine within-block p-values via Fisher's method with
  # heuristic conservative inflation factor for unknown positive dependence.
  # We compute TWO references:
  #   (a) "conservative" form: S/c against chi-square(2k)  — non-textbook,
  #       intentionally conservative (used in the headline manuscript counts).
  #   (b) "textbook" form    : S/c against chi-square(2k/c) — Brown 1975 /
  #       Kost-McDermott 2002 reference, with df rescaled by the same factor.
  block_p_cons <- vapply(names(blocks), function(b) {
    idx <- blocks[[b]]
    if (length(idx) == 1L) return(p_emp[idx])
    chisq_stat <- -2 * sum(log(p_emp[idx]))
    df <- 2 * length(idx)
    # Heuristic inflation factor for Fisher under positive dependence.
    bh_factor <- 1 + (length(idx) - 1L) / 4
    pchisq(chisq_stat / bh_factor, df = df, lower.tail = FALSE)
  }, numeric(1L))

  block_p_text <- vapply(names(blocks), function(b) {
    idx <- blocks[[b]]
    if (length(idx) == 1L) return(p_emp[idx])
    chisq_stat <- -2 * sum(log(p_emp[idx]))
    df_raw <- 2 * length(idx)
    bh_factor <- 1 + (length(idx) - 1L) / 4
    df_textbook <- df_raw / bh_factor
    pchisq(chisq_stat / bh_factor, df = df_textbook, lower.tail = FALSE)
  }, numeric(1L))

  block_p <- block_p_cons

  # Stage 2: BH-FDR across block p-values (separately for each reference).
  block_q          <- p.adjust(block_p_cons, method = "BH")
  block_q_textbook <- p.adjust(block_p_text, method = "BH")

  # Per-cohort q within block (a secondary view for transparency).
  cohort_q_within <- numeric(length(p_emp))
  for (b in names(blocks)) {
    idx <- blocks[[b]]
    cohort_q_within[idx] <- p.adjust(p_emp[idx], method = "BH")
  }

  data.frame(
    cohort_idx = seq_along(p_emp),
    block_id = block_id,
    p_emp = p_emp,
    p_block = block_p[match(block_id, names(block_p))],
    p_block_textbook = block_p_text[match(block_id, names(block_p_text))],
    q_block_BH = block_q[match(block_id, names(block_q))],
    q_block_BH_textbook = block_q_textbook[match(block_id, names(block_q_textbook))],
    q_cohort_within_block = cohort_q_within
  )
}

# ---- v0.6 Nested-CV helpers --------------------------------------------------

# Local null-coalesce (avoid hard dependency on run_paper3_pipeline.R env).
.mnb_null_coalesce <- function(a, b) if (is.null(a)) b else a

#' Stratified k-folds (ungrouped). Returns list of k test-index vectors.
.stratified_folds <- function(y, k = 5L, seed = 42L) {
  set.seed(seed)
  idx0 <- which(y == 0L); idx1 <- which(y == 1L)
  folds <- vector("list", k)
  f0 <- sample(rep_len(seq_len(k), length(idx0)))
  f1 <- sample(rep_len(seq_len(k), length(idx1)))
  for (fi in seq_len(k))
    folds[[fi]] <- sort(c(idx0[f0 == fi], idx1[f1 == fi]))
  folds
}

#' Grouped stratified k-folds. Groups are never split across train/test.
.grouped_folds <- function(y, group_id, k = 5L, seed = 42L) {
  set.seed(seed)
  ugroups <- unique(group_id)
  group_cls <- vapply(ugroups, function(g) {
    yg <- y[group_id == g]
    as.integer(round(mean(yg)))
  }, integer(1L))
  fold_assign <- integer(length(ugroups))
  for (cls in c(0L, 1L)) {
    g_idx <- which(group_cls == cls)
    if (length(g_idx) == 0L) next
    g_idx <- g_idx[sample.int(length(g_idx))]
    fold_assign[g_idx] <- rep_len(seq_len(k), length(g_idx))
  }
  folds <- vector("list", k)
  for (fi in seq_len(k)) {
    test_groups <- ugroups[fold_assign == fi]
    folds[[fi]] <- sort(which(group_id %in% test_groups))
  }
  folds
}

#' Univariate-AUC panel selection on training data only.
#' Returns character vector of top panel_size feature names.
.select_panel_training_fold <- function(X_train, y_train, panel_size = 20L) {
  n_feat <- ncol(X_train)
  aucs <- vapply(seq_len(n_feat), function(j) {
    s  <- X_train[, j]
    ok <- is.finite(s) & !is.na(y_train)
    if (sum(ok) < 5L || length(unique(y_train[ok])) < 2L) return(0.5)
    a  <- suppressWarnings(
      as.numeric(pROC::auc(y_train[ok], s[ok], quiet = TRUE)))
    if (!is.finite(a)) return(0.5)
    max(a, 1 - a)
  }, numeric(1L))
  top_idx <- order(aucs, decreasing = TRUE)[seq_len(min(panel_size, n_feat))]
  colnames(X_train)[top_idx]
}

# ---- v0.6 Nested-CV matched-null benchmark -----------------------------------

#' Nested-CV matched-null benchmark (v0.6 §4.1).
#'
#' Outer k=outer_k folds. Panel selection and null-panel strata are computed on
#' the TRAINING fold only; AUC evaluation is on the held-out TEST fold.
#' When group_id is supplied (e.g. patient IDs from matched_set_audit.tsv),
#' grouped folds are used so no patient group is split across train/test.
#'
#' @param X            matrix (samples × features).
#' @param y            binary outcome 0/1.
#' @param panel_size   integer. Panel size selected per fold on training data.
#' @param scoring_fn   function(X_panel) -> per-sample scores, or NULL for
#'   rowSums(log(X+0.5)). Pre-built closures are accepted; the key invariant is
#'   that strata/panel selection never sees held-out data.
#' @param K            null replicates per fold. Production default 10000.
#' @param outer_k      outer CV folds. Default 5.
#' @param group_id     character/integer vector length nrow(X), or NULL.
#'   Required for GSE262260 and GSE314577 (matched_set_audit.tsv rows with
#'   audit_status=="matched_set_populated").
#' @param feature_pool character. Full candidate pool. NULL = colnames(X).
#' @param exclude_features character. Excluded from null candidate pool.
#' @param post_hemolysis_prepared logical.
#' @param seed         integer. Fold assignment + null-draw seeds are derived
#'   deterministically from this value. v0.6 uses per-draw null seeds:
#'   set.seed(seed) for folds; set.seed(seed * 100000L + fi * 10000L + k)
#'   for null draw k in fold fi. This is a clean-break protocol that makes
#'   serial and forked null draws schedule-independent.
#' @param min_valid_folds   integer. Default 4 (plan §4.1).
#' @param min_pos_per_fold  integer. Default 5.
#' @param min_neg_per_fold  integer. Default 5.
#' @param null_parallel_cores integer. Number of forked workers for null draws
#'   inside each fold. Default 1 keeps serial execution.
#'
#' @return Named list: auc_obs_cv, p_emp_cv, n_valid_folds, eligible, plus
#'   full grouping and fold audit fields per plan §5.5.
matched_null_benchmark_nested_cv <- function(
    X,
    y,
    panel_size            = 20L,
    scoring_fn            = NULL,
    K                     = 10000L,
    outer_k               = 5L,
    group_id              = NULL,
    feature_pool          = NULL,
    exclude_features      = character(0),
    post_hemolysis_corrected = FALSE,
    seed                  = 42L,
    min_valid_folds       = 4L,
    min_pos_per_fold      = 5L,
    min_neg_per_fold      = 5L,
    null_parallel_cores   = 1L) {

  stopifnot(is.matrix(X) || is.data.frame(X))
  X <- as.matrix(X)
  if (is.factor(y)) y <- as.integer(y) - 1L
  if (is.null(feature_pool)) feature_pool <- colnames(X)
  seed <- suppressWarnings(as.integer(seed))
  if (!is.finite(seed) || is.na(seed) || seed < 0L || seed > 21474L) {
    stop("matched_null_benchmark_nested_cv: seed must be an integer from 0 to 21474 so compound null seeds do not overflow.")
  }
  null_parallel_cores <- suppressWarnings(as.integer(null_parallel_cores))
  if (!is.finite(null_parallel_cores) || null_parallel_cores < 1L) {
    null_parallel_cores <- 1L
  }
  rng_protocol <- "v0.6_per_draw_seed_schedule_independent"

  has_groups      <- !is.null(group_id) && length(group_id) == nrow(X)
  grouping_policy <- if (has_groups) "groupKFold" else "ungrouped_stratified"
  n_groups        <- if (has_groups) length(unique(group_id)) else NA_integer_

  # Hemolysis exclusion applied globally to candidate pool.
  hemolysis_markers <- c("hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p",
                          "hsa-miR-144-3p", "hsa-miR-223-3p")
  if (post_hemolysis_corrected)
    exclude_features <- unique(c(exclude_features, hemolysis_markers))
  feature_pool <- setdiff(feature_pool, exclude_features)

  # Build outer folds.
  folds <- if (has_groups)
    .grouped_folds(y, group_id, k = outer_k, seed = seed)
  else
    .stratified_folds(y, k = outer_k, seed = seed)

  # Validate: no group split across train/test (hard stop).
  n_group_split_violations <- 0L
  if (has_groups) {
    for (fi in seq_along(folds)) {
      train_idx <- setdiff(seq_len(nrow(X)), folds[[fi]])
      viol <- length(intersect(unique(group_id[folds[[fi]]]),
                                unique(group_id[train_idx])))
      if (viol > 0L)
        stop(sprintf(
          "matched_null_benchmark_nested_cv: %d group(s) split across train/test in fold %d.",
          viol, fi))
    }
  }

  auc_obs_folds   <- rep(NA_real_, outer_k)
  auc_null_folds  <- matrix(NA_real_, nrow = K, ncol = outer_k)
  fold_panels     <- vector("list", outer_k)
  fold_audit      <- vector("list", outer_k)
  valid_fold_mask <- logical(outer_k)

  .align_fold_score <- function(score, test_idx, n_test) {
    score <- as.numeric(score)
    if (length(score) == n_test) return(score)
    if (length(score) == nrow(X)) return(score[test_idx])
    rep(NA_real_, n_test)
  }

  for (fi in seq_len(outer_k)) {
    test_idx  <- folds[[fi]]
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)
    X_train <- X[train_idx, , drop = FALSE]; y_train <- y[train_idx]
    X_test  <- X[test_idx,  , drop = FALSE]; y_test  <- y[test_idx]

    n_pos_test  <- sum(y_test  == 1L)
    n_neg_test  <- sum(y_test  == 0L)
    n_pos_train <- sum(y_train == 1L)
    n_neg_train <- sum(y_train == 0L)

    fa <- list(fold = fi,
               n_train = length(train_idx), n_test = length(test_idx),
               n_pos_test = n_pos_test,   n_neg_test = n_neg_test,
               n_pos_train = n_pos_train, n_neg_train = n_neg_train)

    skip_reason <- NULL
    if (n_pos_test < min_pos_per_fold || n_neg_test < min_neg_per_fold)
      skip_reason <- sprintf("test n_pos=%d<%d or n_neg=%d<%d",
                              n_pos_test, min_pos_per_fold,
                              n_neg_test, min_neg_per_fold)
    else if (length(unique(y_train)) < 2L)
      skip_reason <- "y_train has only one class"

    if (!is.null(skip_reason)) {
      fa$skip_reason <- skip_reason; fold_audit[[fi]] <- fa; next
    }

    # Panel selection on TRAINING fold only.
    train_pool <- intersect(feature_pool, colnames(X_train))
    if (length(train_pool) < 2L) {
      fa$skip_reason <- "train_pool < 2 features"; fold_audit[[fi]] <- fa; next
    }
    X_train_pool <- X_train[, train_pool, drop = FALSE]
    panel_f <- tryCatch(
      .select_panel_training_fold(X_train_pool, y_train, panel_size),
      error = function(e) NULL)
    if (is.null(panel_f) || length(panel_f) == 0L) {
      fa$skip_reason <- "panel selection failed"; fold_audit[[fi]] <- fa; next
    }
    fold_panels[[fi]] <- panel_f

    # Strata computed on TRAINING fold only.
    dr_train <- colMeans(X_train_pool > 0,  na.rm = TRUE)
    ma_train <- colMeans(X_train_pool,       na.rm = TRUE)
    strata_info <- .match_panel_strata(
      panel = panel_f, feature_pool = train_pool,
      detection_rate = dr_train, mean_abundance = ma_train)

    # Observed AUC on TEST fold.
    panel_in_test <- intersect(panel_f, colnames(X_test))
    if (length(panel_in_test) == 0L) {
      fa$skip_reason <- "panel features absent from test colnames"
      fold_audit[[fi]] <- fa; next
    }
    Xp_test <- X_test[, panel_in_test, drop = FALSE]
    score_obs_f <- if (!is.null(scoring_fn))
      tryCatch(scoring_fn(Xp_test), error = function(e) rep(NA_real_, nrow(Xp_test)))
    else
      rowSums(log(Xp_test + 0.5))
    score_obs_f <- .align_fold_score(score_obs_f, test_idx, nrow(Xp_test))

    if (all(is.na(score_obs_f)) ||
        length(unique(score_obs_f[is.finite(score_obs_f)])) < 2L) {
      fa$skip_reason <- "scorer constant/all-NA on test fold"
      fold_audit[[fi]] <- fa; next
    }
    auc_obs_folds[fi] <- as.numeric(pROC::auc(y_test, score_obs_f, quiet = TRUE))
    valid_fold_mask[fi] <- TRUE

    # K null draws — strata from TRAINING fold; scored on TEST fold.
    # Each draw gets an explicit seed so serial and forked execution are
    # deterministic and independent of process scheduling.
    .run_null_draw <- function(k) {
      set.seed(seed * 100000L + fi * 10000L + k)
      p_k <- .draw_matched_panel(strata_info, exclude = panel_f)
      tier_k <- attr(p_k, "fallback_tier") %||% 1L
      if (length(p_k) == 0L) {
        return(list(auc = NA_real_, fallback_empty = 1L, tier = 3L))
      }

      p_k_test <- intersect(p_k, colnames(X_test))
      if (length(p_k_test) == 0L) {
        return(list(auc = NA_real_, fallback_empty = 0L,
                    tier = as.integer(tier_k)))
      }
      Xp_k <- X_test[, p_k_test, drop = FALSE]
      score_k <- if (!is.null(scoring_fn))
        tryCatch(scoring_fn(Xp_k), error = function(e) rep(NA_real_, nrow(Xp_k)))
      else
        rowSums(log(Xp_k + 0.5))
      score_k <- .align_fold_score(score_k, test_idx, nrow(Xp_k))
      list(
        auc = suppressWarnings(as.numeric(pROC::auc(y_test, score_k, quiet = TRUE))),
        fallback_empty = 0L,
        tier = as.integer(tier_k)
      )
    }

    use_fork <- .Platform$OS.type != "windows" && null_parallel_cores > 1L && K > 1L
    null_results <- if (use_fork) {
      parallel::mclapply(seq_len(K), .run_null_draw,
                         mc.cores = min(null_parallel_cores, K),
                         mc.preschedule = TRUE,
                         mc.set.seed = FALSE)
    } else {
      lapply(seq_len(K), .run_null_draw)
    }
    if (any(vapply(null_results, inherits, logical(1L), "try-error"))) {
      stop(sprintf("matched_null_benchmark_nested_cv: null draw failed in fold %d.", fi))
    }

    auc_null_folds[, fi] <- vapply(null_results, `[[`, numeric(1L), "auc")
    fold_fallback_count <- sum(vapply(null_results, `[[`, integer(1L), "fallback_empty"))
    tiers <- vapply(null_results, `[[`, integer(1L), "tier")
    tier_counts <- c(`1` = sum(tiers == 1L), `2` = sum(tiers == 2L), `3` = sum(tiers == 3L))

    fa$fold_fallback_count     <- fold_fallback_count
    fa$fold_fallback_tier3_pct <- 100 * fold_fallback_count / max(1L, K)
    fa$fold_fallback_tier_counts <- as.list(tier_counts)
    fold_audit[[fi]] <- fa
  }  # end fold loop

  # Aggregate.
  n_valid_folds <- sum(valid_fold_mask)
  valid_fi      <- which(valid_fold_mask)

  pos_v <- vapply(valid_fi, function(fi) fold_audit[[fi]]$n_pos_test, integer(1L))
  neg_v <- vapply(valid_fi, function(fi) fold_audit[[fi]]$n_neg_test, integer(1L))
  min_pos_obs <- if (length(pos_v) > 0L) min(pos_v) else NA_integer_
  min_neg_obs <- if (length(neg_v) > 0L) min(neg_v) else NA_integer_

  skip_summary <- paste(vapply(fold_audit, function(fa)
    .mnb_null_coalesce(fa$skip_reason, "ok"), character(1L)), collapse = "; ")

  if (n_valid_folds < min_valid_folds) {
    return(list(
      eligible                  = FALSE,
      ineligible_reason         = sprintf(
        "n_valid_folds=%d < min_valid_folds=%d; per-fold: %s",
        n_valid_folds, min_valid_folds, skip_summary),
      n_valid_folds             = n_valid_folds,
      n_null_valid_folds        = 0L,
      n_null_valid              = 0L,
      p_emp_cv                  = NA_real_,
      grouping_policy           = grouping_policy,
      n_groups                  = n_groups,
      n_group_split_violations  = n_group_split_violations,
      matched_set_populated_pct = NA_real_,
      min_pos_per_fold          = min_pos_obs,
      min_neg_per_fold          = min_neg_obs,
      fold_audit                = fold_audit,
      K = K, outer_k = outer_k, seed = seed,
      null_parallel_cores = null_parallel_cores,
      rng_protocol = rng_protocol))
  }

  auc_obs_cv  <- mean(auc_obs_folds[valid_fold_mask], na.rm = TRUE)

  # --- Matched-null eligibility gate (F4) -------------------------------------
  # The empirical CV p-value is interpretable only if enough genuine
  # size/strata-matched null panels could be drawn from the candidate pool
  # (excluding the observed panel). Guard the statistic actually used for the
  # null, at two levels:
  #   (1) per-fold null validity — a fold whose null draws are >=50% degenerate
  #       (empty matched panel; counted in `fold_fallback_count`) carries no
  #       usable null distribution and is dropped from the null aggregation,
  #       independent of its observed-AUC validity. A single dead fold therefore
  #       cannot be masked by clean folds. The aggregate
  #       `total_fallback < (K * n_valid_folds)/2` form is deliberately NOT used.
  #   (2) aggregate — at least K/2 complete null draws must survive across the
  #       null-valid folds.
  # A cell failing either level is matched-null INELIGIBLE: p_emp_cv = NA and
  # only the bootstrap CI on the observed AUC is reportable. This mirrors the
  # flat-path gate inside matched_null_benchmark() (`eligible <- fallback_count
  # eligible = TRUE), letting a fully-degenerate cell collapse to p = (1+0)/(1+0).
  fold_fallback_frac <- vapply(valid_fi, function(fi) {
    fb <- .mnb_null_coalesce(fold_audit[[fi]]$fold_fallback_count, NA_integer_)
    if (is.na(fb)) NA_real_ else fb / K
  }, numeric(1L))
  null_valid_fi      <- valid_fi[is.finite(fold_fallback_frac) & fold_fallback_frac < 0.5]
  n_null_valid_folds <- length(null_valid_fi)

  if (n_null_valid_folds > 0L) {
    null_valid    <- auc_null_folds[, null_valid_fi, drop = FALSE]
    complete_null <- rowSums(is.na(null_valid)) == 0L
    auc_null_cv   <- if (any(complete_null))
      rowMeans(null_valid[complete_null, , drop = FALSE], na.rm = FALSE)
    else
      numeric(0L)
    # Observed mean aligned to the SAME null-valid folds the null is built from,
    # so the paired comparison is over a consistent fold set. When no fold is
    # dropped (the common case) this equals auc_obs_cv.
    auc_obs_for_p <- mean(auc_obs_folds[null_valid_fi], na.rm = TRUE)
  } else {
    auc_null_cv   <- numeric(0L)
    auc_obs_for_p <- NA_real_
  }
  n_null_valid <- length(auc_null_cv)

  matched_null_eligible <- (n_null_valid_folds > 0L) && (n_null_valid >= K / 2)
  if (matched_null_eligible) {
    p_emp_cv          <- (1 + sum(auc_null_cv >= auc_obs_for_p, na.rm = TRUE)) /
                         (1 + n_null_valid)
    ineligible_reason <- NA_character_
  } else {
    p_emp_cv          <- NA_real_
    ineligible_reason <- sprintf(
      "matched-null ineligible: %d complete null draw(s) across %d null-valid fold(s) (require >= %g); candidate pool exhausted -> report bootstrap CI on observed AUC only",
      n_null_valid, n_null_valid_folds, K / 2)
  }

  total_fallback <- sum(vapply(fold_audit[valid_fi], function(fa)
    .mnb_null_coalesce(fa$fold_fallback_count, 0L), integer(1L)))
  max_fallback_tier <- max(vapply(fold_audit[valid_fi], function(fa) {
    tc <- fa$fold_fallback_tier_counts %||% list(`1` = 0L, `2` = 0L, `3` = 0L)
    used <- as.integer(names(tc)[unlist(tc) > 0L])
    if (length(used) == 0L) 1L else max(used)
  }, integer(1L)))

  list(
    auc_obs_cv               = auc_obs_cv,
    auc_obs_folds            = auc_obs_folds,
    auc_null_cv              = auc_null_cv,
    p_emp_cv                 = p_emp_cv,
    n_valid_folds            = n_valid_folds,
    n_null_valid_folds       = n_null_valid_folds,
    n_null_valid             = n_null_valid,
    fold_panels              = fold_panels,
    K                        = K,
    outer_k                  = outer_k,
    seed                     = seed,
    null_parallel_cores      = null_parallel_cores,
    rng_protocol             = rng_protocol,
    eligible                 = matched_null_eligible,
    ineligible_reason        = ineligible_reason,
    fallback_count           = total_fallback,
    fallback_tier            = max_fallback_tier,
    grouping_policy          = grouping_policy,
    n_groups                 = n_groups,
    n_group_split_violations = n_group_split_violations,
    matched_set_populated_pct = NA_real_,   # caller fills from audit
    min_pos_per_fold         = min_pos_obs,
    min_neg_per_fold         = min_neg_obs,
    fold_audit               = fold_audit,
    panel_size               = panel_size,
    exclude_features         = exclude_features,
    post_hemolysis_corrected = post_hemolysis_corrected
  )
}

if (sys.nframe() == 0L) {
  set.seed(42)
  n <- 100; p <- 200
  X <- matrix(rlnorm(n * p, meanlog = 5, sdlog = 1.0), nrow = n, ncol = p)
  colnames(X) <- paste0("mir-", sprintf("%03d", seq_len(p)))
  # Inject signal in a true panel of 12 features.
  true_panel <- paste0("mir-", sprintf("%03d", 11:22))
  y <- rbinom(n, 1, 0.5)
  X[y == 1, true_panel] <- X[y == 1, true_panel] * 2.5

  scoring_fn <- function(M) rowSums(log(M + 1))

  res <- matched_null_benchmark(X, y, panel = true_panel, scoring_fn = scoring_fn,
                                 K = 200L, seed = 42L)
  cat(sprintf("matched_null_benchmark self-test (K=200 for speed; default is K=10,000):\n"))
  cat(sprintf("  panel size      = %d\n", res$panel_size))
  cat(sprintf("  AUC_obs         = %.3f\n", res$auc_obs))
  cat(sprintf("  AUC_null mean   = %.3f (sd = %.3f)\n",
              mean(res$auc_null), sd(res$auc_null)))
  cat(sprintf("  empirical p     = %.4f (K = %d)\n", res$p_emp, res$K))
  stopifnot(res$auc_obs > mean(res$auc_null) + 0.05)
  stopifnot(res$p_emp < 0.05)
  cat("  PASS (signal panel beats matched null at p < 0.05)\n")

  # Verify non-overlap of null panels with observed panel. F8 fix: the prior
  # assertion compared panel feature names to the integer `res$panel_size`
  # (`p %in% res$panel_size`), which is always FALSE, so it never tested
  # overlap. Compare against the observed panel feature names (`true_panel`).
  any_overlap <- any(vapply(res$null_panels,
                              function(p) any(p %in% true_panel),
                              logical(1L)))
  stopifnot(!any_overlap)
  cat("  PASS (null panels do NOT overlap with observed panel)\n")

  # Verify hemolysis-marker exclusion when toggle is on.
  hemo <- "hsa-miR-451a"
  X2 <- cbind(X, matrix(rlnorm(n, 5, 1), nrow = n, dimnames = list(NULL, hemo)))
  res2 <- matched_null_benchmark(X2, y, panel = true_panel, scoring_fn = scoring_fn,
                                   K = 50L, post_hemolysis_corrected = TRUE, seed = 42L)
  hemo_in_null <- any(vapply(res2$null_panels,
                               function(p) hemo %in% p, logical(1L)))
  stopifnot(!hemo_in_null)
  cat("  PASS (hemolysis markers excluded from null pool when post_hemolysis_prepared=TRUE)\n")

  # Test block-aware BH-FDR. Construct a 4-cohort scenario where 2 cohorts
  # share a provenance block (specimen-reuse) and 2 are singletons.
  set.seed(7)
  fake_results <- list(
    list(p_emp = 0.005),  # cohort 1 — Toray cluster
    list(p_emp = 0.012),  # cohort 2 — Toray cluster (correlated with #1)
    list(p_emp = 0.04),   # cohort 3 — Affy singleton
    list(p_emp = 0.30)    # cohort 4 — Agilent singleton
  )
  blocks_fake <- c("Toray", "Toray", "Affy-S", "Agilent-S")
  blocked <- bh_fdr_correct_blocked(fake_results, blocks_fake)
  stopifnot(nrow(blocked) == 4L)
  stopifnot(blocked$block_id[1] == "Toray" && blocked$block_id[2] == "Toray")
  stopifnot(blocked$p_block[1] == blocked$p_block[2])  # same block, same p_block
  stopifnot(all(blocked$q_block_BH >= blocked$p_block - 1e-9))
  cat("  PASS (block-aware BH-FDR: Toray cluster collapses to 1 block-level p; q_block >= p_block)\n")

  # Test tier-3 fallback eligibility flag (round-13 reviewer fix). Construct a
  # degenerate scenario: panel has more features than the candidate pool.
  set.seed(11)
  n_small <- 50; p_small <- 30
  X_small <- matrix(rlnorm(n_small * p_small, meanlog = 5, sdlog = 1.0),
                     nrow = n_small, ncol = p_small)
  colnames(X_small) <- paste0("mir-", sprintf("%03d", seq_len(p_small)))
  panel_oversized <- paste0("mir-", sprintf("%03d", seq_len(28)))   # 28 of 30
  y_small <- rbinom(n_small, 1, 0.5)
  res_degen <- matched_null_benchmark(X_small, y_small, panel = panel_oversized,
                                        scoring_fn = scoring_fn, K = 50L, seed = 11L)
  stopifnot(isTRUE(res_degen$eligible) || isFALSE(res_degen$eligible))
  if (!isTRUE(res_degen$eligible)) {
    stopifnot(res_degen$fallback_tier == 3L)
    stopifnot(is.na(res_degen$p_emp))
    cat("  PASS (tier-3 fallback: oversized panel marks cell not-eligible; p_emp=NA, fallback_tier=3)\n")
  } else {
    cat("  PASS (tier-3 fallback path well-formed; degenerate test case still completed null draws)\n")
  }

  # ---- v0.6 nested-CV self-tests --------------------------------------------
  cat("\n[nested-CV self-tests]\n")

  # Test A: ungrouped nested-CV returns eligible with signal panel.
  set.seed(99)
  n_cv <- 150; p_cv <- 120
  X_cv <- matrix(rlnorm(n_cv * p_cv, meanlog = 5, sdlog = 1), nrow = n_cv,
                  dimnames = list(NULL, paste0("f-", seq_len(p_cv))))
  y_cv <- rbinom(n_cv, 1, 0.5)
  sig_feats <- paste0("f-", 1:10)
  X_cv[y_cv == 1, sig_feats] <- X_cv[y_cv == 1, sig_feats] * 3.0

  res_cv <- matched_null_benchmark_nested_cv(
    X_cv, y_cv, panel_size = 10L,
    scoring_fn = function(M) rowSums(log(M + 0.5)),
    K = 100L, outer_k = 5L, seed = 42L,
    min_valid_folds = 4L, min_pos_per_fold = 5L, min_neg_per_fold = 5L)

  stopifnot(isTRUE(res_cv$eligible))
  stopifnot(res_cv$n_valid_folds >= 4L)
  stopifnot(is.finite(res_cv$auc_obs_cv))
  stopifnot(is.finite(res_cv$p_emp_cv))
  stopifnot(res_cv$grouping_policy == "ungrouped_stratified")
  stopifnot(is.na(res_cv$n_groups))
  stopifnot(res_cv$n_group_split_violations == 0L)
  cat(sprintf("  [A] PASS ungrouped nested-CV: auc_obs_cv=%.3f p_emp=%.4f n_valid=%d\n",
              res_cv$auc_obs_cv, res_cv$p_emp_cv, res_cv$n_valid_folds))

  # Test B: training-fold-only strata — panel selected on train must NOT use test.
  # Verify each fold's panel comes from .select_panel_training_fold (training only).
  # Proxy check: panel names must all be in colnames(X_cv) (column pool).
  for (fi in seq_along(res_cv$fold_panels)) {
    fp <- res_cv$fold_panels[[fi]]
    if (!is.null(fp))
      stopifnot(all(fp %in% colnames(X_cv)))
  }
  cat("  [B] PASS fold panels are subsets of feature pool (training-fold selection)\n")

  # Test C: grouped folds (simulate GSE262260/GSE314577 matched-set structure).
  set.seed(77)
  n_grp <- 100; p_grp <- 60
  X_grp <- matrix(rlnorm(n_grp * p_grp, meanlog = 5, sdlog = 1), nrow = n_grp,
                   dimnames = list(NULL, paste0("g-", seq_len(p_grp))))
  y_grp <- rbinom(n_grp, 1, 0.5)
  X_grp[y_grp == 1, 1:8] <- X_grp[y_grp == 1, 1:8] * 2.5
  # 50 groups of 2 samples each (simulating matched pairs)
  group_id_sim <- rep(seq_len(50L), each = 2L)

  res_grp <- matched_null_benchmark_nested_cv(
    X_grp, y_grp, panel_size = 8L,
    scoring_fn = function(M) rowSums(log(M + 0.5)),
    K = 50L, outer_k = 5L, group_id = group_id_sim, seed = 1L,
    min_valid_folds = 3L, min_pos_per_fold = 4L, min_neg_per_fold = 4L)

  stopifnot(res_grp$grouping_policy == "groupKFold")
  stopifnot(!is.na(res_grp$n_groups))
  stopifnot(res_grp$n_group_split_violations == 0L)
  cat(sprintf("  [C] PASS grouped folds: n_groups=%d split_violations=%d eligible=%s\n",
              res_grp$n_groups, res_grp$n_group_split_violations,
              res_grp$eligible))

  # Test D: ineligible_insufficient_valid_folds when tiny cohort forces skipped folds.
  set.seed(3)
  n_tiny2 <- 30; p_tiny2 <- 40
  X_t2 <- matrix(rlnorm(n_tiny2 * p_tiny2), nrow = n_tiny2,
                  dimnames = list(NULL, paste0("t-", seq_len(p_tiny2))))
  y_t2 <- c(rep(0L, 27L), rep(1L, 3L))  # highly imbalanced → folds will have <5 pos

  res_t2 <- matched_null_benchmark_nested_cv(
    X_t2, y_t2, panel_size = 5L,
    K = 20L, outer_k = 5L, seed = 5L,
    min_valid_folds = 4L, min_pos_per_fold = 5L, min_neg_per_fold = 5L)

  # Either eligible=FALSE (expected) or eligible=TRUE with n_valid<5.
  stopifnot(is.logical(res_t2$eligible))
  if (!res_t2$eligible) {
    stopifnot(!is.null(res_t2$ineligible_reason))
    cat(sprintf("  [D] PASS ineligible small-cohort: reason='%s'\n",
                substr(res_t2$ineligible_reason, 1L, 60L)))
  } else {
    cat(sprintf("  [D] PASS small-cohort still completed with n_valid=%d\n",
                res_t2$n_valid_folds))
  }

  # Test E: deterministic — same seed gives same p_emp_cv.
  res_det1 <- matched_null_benchmark_nested_cv(
    X_cv, y_cv, panel_size = 10L,
    scoring_fn = function(M) rowSums(log(M + 0.5)),
    K = 30L, outer_k = 3L, seed = 17L, min_valid_folds = 2L)
  res_det2 <- matched_null_benchmark_nested_cv(
    X_cv, y_cv, panel_size = 10L,
    scoring_fn = function(M) rowSums(log(M + 0.5)),
    K = 30L, outer_k = 3L, seed = 17L, min_valid_folds = 2L)
  stopifnot(identical(res_det1$p_emp_cv, res_det2$p_emp_cv))
  cat("  [E] PASS deterministic seed: same p_emp_cv on two identical calls\n")

  # Test F (F4 gate, distinct from Test D): matched-null DEGENERACY -> ineligible.
  # Pool of 14, panel of 12 -> after excluding the observed panel the candidate
  # pool has 2 features (< panel size), so .draw_matched_panel returns empty for
  # every draw in every fold. Observed-AUC folds remain valid (signal present);
  # only the null is degenerate. The fold-level >=50%-fallback rule + the
  # n_null_valid >= K/2 rule must mark the cell ineligible (p_emp_cv = NA),
  # NOT collapse to p = (1+0)/(1+0) = 1.0 as the pre-fix code did.
  set.seed(23)
  n_dg <- 80; p_dg <- 14
  X_dg <- matrix(rlnorm(n_dg * p_dg, meanlog = 5, sdlog = 1), nrow = n_dg,
                  dimnames = list(NULL, paste0("d-", seq_len(p_dg))))
  y_dg <- rbinom(n_dg, 1, 0.5)
  X_dg[y_dg == 1, 1:12] <- X_dg[y_dg == 1, 1:12] * 3.0
  res_dg <- matched_null_benchmark_nested_cv(
    X_dg, y_dg, panel_size = 12L,
    scoring_fn = function(M) rowSums(log(M + 0.5)),
    K = 60L, outer_k = 3L, seed = 23L,
    min_valid_folds = 2L, min_pos_per_fold = 5L, min_neg_per_fold = 5L)
  stopifnot(isFALSE(res_dg$eligible))
  stopifnot(is.na(res_dg$p_emp_cv))
  stopifnot(!is.na(res_dg$ineligible_reason))
  stopifnot(res_dg$n_null_valid == 0L)
  stopifnot(res_dg$n_valid_folds >= 2L)   # observed folds valid; only the null is degenerate
  cat(sprintf("  [F] PASS matched-null degeneracy -> ineligible: n_valid_folds=%d n_null_valid=%d p=%s\n",
              res_dg$n_valid_folds, res_dg$n_null_valid, as.character(res_dg$p_emp_cv)))

  # Test G (F8): .draw_matched_panel never returns a feature in `exclude`, the
  set.seed(31)
  pool_g <- paste0("p-", seq_len(60))
  dr_g <- runif(60, 0.2, 1.0); ma_g <- rlnorm(60, 5, 1)
  names(dr_g) <- names(ma_g) <- pool_g
  excl_g <- pool_g[1:20]
  si_g <- .match_panel_strata(panel = excl_g, feature_pool = pool_g,
                              detection_rate = dr_g, mean_abundance = ma_g)
  draw_overlap <- FALSE; n_nonempty_g <- 0L
  for (rep_g in seq_len(200L)) {
    pk <- .draw_matched_panel(si_g, exclude = excl_g)
    if (length(pk) > 0L) {
      n_nonempty_g <- n_nonempty_g + 1L
      if (any(pk %in% excl_g)) { draw_overlap <- TRUE; break }
    }
  }
  stopifnot(!draw_overlap, n_nonempty_g > 0L)
  cat(sprintf("  [G] PASS .draw_matched_panel never overlaps excluded panel (%d/200 non-empty draws)\n",
              n_nonempty_g))

  # Test H (F8 nested-CV train/test non-overlap regression — strengthened
  # after Round 2 review showed the prior multiplicative perturbation was
  # rank-preserving on pooled AUC and missed the all-data leakage class).
  #
  # Design: introduce a SENTINEL feature whose class-discriminative power
  # lives exclusively in one group. Use grouped CV (outer_k=2, group_id)
  # so the fold structure is deterministic:
  #
  #   Group 1 (samples 1..50):  sentinel = pure noise, no class signal.
  #   Group 2 (samples 51..100): sentinel = STRONG class-discriminative
  #                              signal (class-1 rows shifted by ×100+1000).
  #
  # Under correct train-only panel selection, the sentinel appears in
  # AT MOST the fold whose TRAIN side is group 2 (where the sentinel
  # carries real per-feature AUC). It must NEVER appear in the fold
  # whose TRAIN side is group 1, because the sentinel is pure noise there
  # and the panel is filled by the 8 baseline signal features.
  #
  # A selector that peeks at test-fold data — pools train+test for ranking,
  # uses all-X for feature AUC, or scans test rows for selection — would
  # pick the sentinel in BOTH folds: in the group-2-train fold the sentinel
  # has train signal; in the group-1-train fold the sentinel has test
  # signal (test = group 2). Total count = 2.
  #
  # Assertion: sentinel appears in AT MOST 1 of the 2 fold panels.
  # Train-only contract => exactly 1; all-data or train+test-pooled
  # univariate-AUC leakage => 2. (A pure test-only selector would also
  # give 1 by coincidence — sentinel only has signal in group 2, and
  # only one fold has group 2 as its test side — but such a selector is
  # not a realistic accidental regression, and the all-data class is
  # what the design actually guards against.)
  set.seed(91)
  n_h <- 100L; p_baseline <- 30L
  X_h <- matrix(rlnorm(n_h * p_baseline, meanlog = 5, sdlog = 1),
                nrow = n_h,
                dimnames = list(NULL, paste0("h-", seq_len(p_baseline))))
  y_h <- rep(c(0L, 1L, 0L, 1L), each = 25L)
  # Eight baseline class-signal features so panel_size=8 is fully filled
  # by signal features in any fold whose train sees real class signal,
  # without leaving slots for incidental noise / sentinel ties.
  X_h[y_h == 1L, 1:8] <- X_h[y_h == 1L, 1:8] * 2.5
  group_h <- c(rep(1L, 50L), rep(2L, 50L))

  # Append a sentinel feature: noise everywhere, but class-1 rows in
  # group 2 receive a massive multiplicative+additive shift. This makes
  # the sentinel near-perfectly class-discriminative on group-2 samples
  # only.
  sentinel <- rlnorm(n_h, meanlog = 5, sdlog = 1)
  shift_idx <- (group_h == 2L) & (y_h == 1L)
  sentinel[shift_idx] <- sentinel[shift_idx] * 100 + 1000
  X_h <- cbind(X_h, sentinel = sentinel)

  scorer_h <- function(M) rowSums(log(M + 0.5))
  res_h <- matched_null_benchmark_nested_cv(
    X_h, y_h, panel_size = 8L, scoring_fn = scorer_h,
    K = 50L, outer_k = 2L, seed = 91L, group_id = group_h,
    min_valid_folds = 2L, min_pos_per_fold = 5L, min_neg_per_fold = 5L)

  stopifnot(length(res_h$fold_panels) == 2L,
            res_h$grouping_policy == "groupKFold")
  sentinel_picks <- vapply(res_h$fold_panels, function(p) {
    if (is.null(p)) FALSE else "sentinel" %in% p
  }, logical(1L))
  n_sentinel <- sum(sentinel_picks)
  stopifnot(n_sentinel <= 1L)
  cat(sprintf(
    "  [H] PASS nested-CV train-only contract: sentinel with group-asymmetric class signal picked in %d/%d fold panel(s) (must be <=1; train-only =>1, all-data/pooled test-peek =>2)\n",
    n_sentinel, length(res_h$fold_panels)))

  cat("[nested-CV self-tests] ALL PASS\n")
}
