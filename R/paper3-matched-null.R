#' @title Paper 3 matched-null benchmark for panel-vs-random-panel AUC inference
#'
#' @description
#' Generalised matched-null benchmark introduced in Paper 3 (Module A
#' validation framework; Stawiski et al., in preparation). Unlike the simpler
#' \code{\link{os_panel_null_benchmark}} in \code{panel-gates.R}, this
#' implementation stratifies random-panel draws by per-feature detection-rate
#' and log-mean-abundance quartile bins, explicitly excludes hemolysis-marker
#' miRNAs when \code{post_hemolysis_corrected = TRUE}, and provides a
#' three-tier fallback protocol for small platforms with exhausted candidate
#' pools. BH-FDR and block-aware BH-FDR correction helpers are also provided
#' for multi-modality family testing.
#'
#' Methods provided:
#' \itemize{
#'   \item \code{\link{paper3_matched_null_benchmark}}: main benchmark function.
#'   \item \code{\link{paper3_hanley_mcneil_auc_ci}}: closed-form
#'     Hanley-McNeil 1982 95\% AUC confidence interval helper.
#'   \item \code{\link{paper3_bh_fdr_correct_matched_null}}: BH-FDR within a
#'     modality family.
#'   \item \code{\link{paper3_holm_correct_familywise}}: Holm correction for
#'     family-wise method × claim-state contrasts.
#'   \item \code{\link{paper3_bh_fdr_correct_blocked}}: block-aware two-stage
#'     BH-FDR for specimen-shared cohort clusters.
#' }
#'
#' @references
#' Benjamini Y, Hochberg Y. (1995) Controlling the False Discovery Rate:
#' A Practical and Powerful Approach to Multiple Testing.
#' \emph{Journal of the Royal Statistical Society Series B} 57(1): 289–300.
#'
#' Stawiski K. (in preparation) Provenance-aware within-sample scoring for
#' circulating-microRNA biomarkers across cancers and platforms (Paper 3 of
#' the OmicSelector programme; Nature Methods target).
#'
#' @name paper3-matched-null
NULL


# ----------------------------------------------------------------------------
# .paper3_match_panel_strata — internal helper
# ----------------------------------------------------------------------------

#' @noRd
.paper3_match_panel_strata <- function(panel,
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


# ----------------------------------------------------------------------------
# .paper3_draw_matched_panel — internal helper
# ----------------------------------------------------------------------------

#' @noRd
.paper3_draw_matched_panel <- function(strata_info, exclude = character(0)) {
  strata <- strata_info$strata
  cells <- strata_info$cells_in_panel
  panel_size <- sum(unlist(cells))
  candidates_total <- length(setdiff(names(strata), exclude))
  if (candidates_total < panel_size) return(character(0))   # tier 3
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
  if (shortfall > 0L) {                                       # tier 2
    extra_pool <- setdiff(names(strata), c(exclude, picked))
    if (length(extra_pool) >= shortfall) {
      picked <- c(picked, sample(extra_pool, shortfall))
    } else {
      return(character(0))                                    # degenerate
    }
  }
  picked
}


# ----------------------------------------------------------------------------
# paper3_matched_null_benchmark
# ----------------------------------------------------------------------------

#' @title Paper 3 matched-null benchmark for within-sample miRNA panel scoring
#'
#' @description
#' For an observed scoring method \code{scoring_fn} applied to a panel of
#' size \code{k}, generates \code{K} random panels of the same size drawn from
#' the cohort feature pool, jointly stratified on per-feature detection-rate and
#' log-mean-abundance quartile bins (16 strata total). The empirical p-value
#' measures whether the observed panel AUC exceeds the matched-null distribution.
#'
#' Name note: the simpler (non-stratified) package predecessor is
#' \code{\link{os_panel_null_benchmark}} (\code{panel-gates.R}). This function
#' uses the name \code{paper3_matched_null_benchmark} to avoid collision.
#'
#' @param X Matrix (samples \eqn{\times} features). Must have column names.
#'   Should be the full feature pool, not just the observed panel.
#' @param y Binary outcome (0/1 or factor with 2 levels).
#' @param panel Character vector of feature names in the observed panel.
#' @param scoring_fn Function \code{function(X_subset)} returning a per-sample
#'   numeric score.
#' @param K Integer — number of null panels. Default 10,000 (per plan v0.3.3).
#' @param feature_pool Character vector — full feature pool. Defaults to
#'   \code{colnames(X)}.
#' @param exclude_features Character — feature names to exclude from the
#'   candidate pool (e.g., hemolysis markers). Default \code{character(0)}.
#' @param post_hemolysis_corrected Logical. If \code{TRUE}, prepends the
#'   canonical hemolysis-marker list to \code{exclude_features}. Default
#'   \code{FALSE}.
#' @param seed Integer — RNG seed. Default 42.
#' @param auc_ci_method Character; one of \code{"hanley_mcneil"} (default;
#'   closed-form Hanley-McNeil 1982 standard error), \code{"delong"}
#'   (\code{pROC::ci.auc(method="delong")}; requires \code{pROC} >= 1.16),
#'   or \code{"none"} (skip CI). The 95\% CI is reported alongside the
#'   observed AUC as \code{auc_obs_ci_lo} / \code{auc_obs_ci_hi}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{auc_obs}}{Observed panel AUC.}
#'     \item{\code{auc_obs_ci_lo}, \code{auc_obs_ci_hi}}{95\% confidence
#'       interval for the observed AUC; NA when
#'       \code{auc_ci_method = "none"} or computation fails.}
#'     \item{\code{auc_null}}{Length-\code{K} numeric vector of null AUCs
#'       (NA for degenerate draws).}
#'     \item{\code{p_emp}}{Empirical p-value, or \code{NA_real_} if not
#'       eligible.}
#'     \item{\code{panel_size}}{Number of features in the observed panel.}
#'     \item{\code{K}, \code{seed}}{Parameters as supplied.}
#'     \item{\code{null_panels}}{List of length \code{K} of random feature
#'       name sets.}
#'     \item{\code{excluded_features}}{Final exclusion list applied.}
#'     \item{\code{eligible}}{Logical: whether matched-null p-value is valid.}
#'     \item{\code{fallback_tier}}{1, 2, or 3 (3 = not eligible).}
#'     \item{\code{fallback_count}}{Number of draws that fell back to tier 3.}
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X <- matrix(rlnorm(100 * 50, meanlog = 5, sdlog = 1), nrow = 100,
#'             dimnames = list(NULL, paste0("mir-", seq_len(50))))
#' y <- rbinom(100, 1, 0.5)
#' panel <- paste0("mir-", 1:10)
#' res <- paper3_matched_null_benchmark(X, y, panel,
#'          scoring_fn = function(M) rowSums(log(M + 1)),
#'          K = 200L)
#' res$p_emp
#' }
#'
#' @references
#' Stawiski K. (in preparation) Paper 3 of the OmicSelector programme.
#'
#' @export
paper3_matched_null_benchmark <- function(X, y, panel, scoring_fn,
                                           K = 10000L,
                                           feature_pool = NULL,
                                           exclude_features = character(0),
                                           post_hemolysis_corrected = FALSE,
                                           seed = 42L,
                                           auc_ci_method = c("hanley_mcneil",
                                                             "delong",
                                                             "none")) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("paper3_matched_null_benchmark requires pROC for AUC computation.")
  }
  auc_ci_method <- match.arg(auc_ci_method)
  stopifnot(is.matrix(X) || is.data.frame(X))
  X <- as.matrix(X)
  if (is.null(feature_pool)) feature_pool <- colnames(X)
  if (is.factor(y)) y <- as.integer(y) - 1L

  panel <- intersect(panel, colnames(X))
  if (length(panel) == 0L) stop("None of the panel features are in colnames(X).")

  hemolysis_markers <- c("hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p",
                          "hsa-miR-144-3p", "hsa-miR-223-3p")
  if (post_hemolysis_corrected) {
    exclude_features <- unique(c(exclude_features, hemolysis_markers))
  }
  feature_pool <- setdiff(feature_pool, exclude_features)
  excluded_in_pool <- intersect(panel, exclude_features)
  if (length(excluded_in_pool) > 0L) {
    warning("paper3_matched_null_benchmark: ", length(excluded_in_pool),
             " observed-panel features overlap with the exclusion list. Removing from panel.")
    panel <- setdiff(panel, excluded_in_pool)
  }

  detection_rate <- colMeans(X[, feature_pool, drop = FALSE] > 0, na.rm = TRUE)
  mean_abundance <- colMeans(X[, feature_pool, drop = FALSE], na.rm = TRUE)

  strata_info <- .paper3_match_panel_strata(
    panel = panel, feature_pool = feature_pool,
    detection_rate = detection_rate, mean_abundance = mean_abundance
  )

  score_obs <- scoring_fn(X[, panel, drop = FALSE])
  auc_obs <- as.numeric(pROC::auc(y, score_obs, quiet = TRUE))

  auc_ci <- .paper3_auc_ci(y, score_obs, auc_obs, auc_ci_method)
  auc_obs_ci_lo <- auc_ci[1L]
  auc_obs_ci_hi <- auc_ci[2L]

  set.seed(seed)
  null_panels <- vector("list", K)
  auc_null <- numeric(K)
  fallback_count <- 0L
  for (k in seq_len(K)) {
    p_k <- .paper3_draw_matched_panel(strata_info, exclude = panel)
    if (length(p_k) == 0L) {
      fallback_count <- fallback_count + 1L
      auc_null[k] <- NA_real_
      next
    }
    null_panels[[k]] <- p_k
    score_k <- scoring_fn(X[, p_k, drop = FALSE])
    auc_null[k] <- as.numeric(pROC::auc(y, score_k, quiet = TRUE))
  }

  eligible <- fallback_count < (K / 2L)
  if (!eligible) {
    return(list(
      auc_obs = auc_obs,
      auc_obs_ci_lo = auc_obs_ci_lo,
      auc_obs_ci_hi = auc_obs_ci_hi,
      auc_ci_method = auc_ci_method,
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
      reason = "candidate pool exhausted; report bootstrap CI on observed AUC only"
    ))
  }

  p_emp <- (1 + sum(auc_null >= auc_obs, na.rm = TRUE)) / (K - fallback_count + 1L)

  list(
    auc_obs = auc_obs,
    auc_obs_ci_lo = auc_obs_ci_lo,
    auc_obs_ci_hi = auc_obs_ci_hi,
    auc_ci_method = auc_ci_method,
    auc_null = auc_null,
    p_emp = p_emp,
    panel_size = length(panel),
    K = K,
    seed = seed,
    null_panels = null_panels,
    excluded_features = exclude_features,
    post_hemolysis_corrected = post_hemolysis_corrected,
    eligible = TRUE,
    fallback_tier = if (fallback_count > 0L) 2L else 1L,
    fallback_count = fallback_count
  )
}


# ----------------------------------------------------------------------------
# .paper3_auc_ci -- shared internal AUC CI helper
# ----------------------------------------------------------------------------

#' @noRd
.paper3_auc_ci <- function(y, scores, auc_obs,
                            method = c("hanley_mcneil", "delong", "none")) {
  method <- match.arg(method)
  if (method == "none" || !is.finite(auc_obs)) {
    return(c(NA_real_, NA_real_))
  }
  y_int <- if (is.factor(y)) as.integer(y) - 1L else as.integer(y)
  if (length(unique(y_int[!is.na(y_int)])) < 2L) {
    return(c(NA_real_, NA_real_))
  }

  if (method == "delong") {
    if (!requireNamespace("pROC", quietly = TRUE)) {
      return(c(NA_real_, NA_real_))
    }
    ci <- tryCatch(
      pROC::ci.auc(y, scores, method = "delong", quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(ci) || length(ci) < 3L) return(c(NA_real_, NA_real_))
    return(c(as.numeric(ci[1L]), as.numeric(ci[3L])))
  }

  n_pos <- sum(y_int == 1L, na.rm = TRUE)
  n_neg <- sum(y_int == 0L, na.rm = TRUE)
  ci <- paper3_hanley_mcneil_auc_ci(auc_obs, n_pos, n_neg)
  c(unname(ci[["ci_lo"]]), unname(ci[["ci_hi"]]))
}


#' @title Hanley-McNeil 1982 AUC confidence interval
#'
#' @description
#' Computes the closed-form Hanley-McNeil 1982 standard error and Wald
#' confidence interval for a binary-classification AUC from the AUC estimate
#' and the number of positive and negative samples. This is the manuscript's
#' default per-cell AUC interval helper; it is conservative for cross-validated
#' estimates because it does not model fold dependence.
#'
#' @param auc Numeric scalar in [0, 1].
#' @param n_pos Integer count of positive/case samples.
#' @param n_neg Integer count of negative/control samples.
#' @param conf_level Numeric confidence level in (0, 1). Default 0.95.
#'
#' @return Named numeric vector with \code{ci_lo}, \code{ci_hi}, and
#'   \code{se}. Returns \code{NA} values when inputs are not estimable.
#'
#' @references
#' Hanley JA, McNeil BJ. (1982) The meaning and use of the area under a
#' receiver operating characteristic (ROC) curve. \emph{Radiology}
#' 143(1): 29-36.
#'
#' @examples
#' paper3_hanley_mcneil_auc_ci(auc = 0.72, n_pos = 50, n_neg = 50)
#'
#' @export
paper3_hanley_mcneil_auc_ci <- function(auc, n_pos, n_neg,
                                         conf_level = 0.95) {
  bad <- !is.numeric(auc) || length(auc) != 1L || !is.finite(auc) ||
    auc < 0 || auc > 1 ||
    !is.numeric(n_pos) || length(n_pos) != 1L || !is.finite(n_pos) ||
    !is.numeric(n_neg) || length(n_neg) != 1L || !is.finite(n_neg) ||
    !is.numeric(conf_level) || length(conf_level) != 1L ||
    !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1
  if (bad || n_pos < 1L || n_neg < 1L) {
    return(c(ci_lo = NA_real_, ci_hi = NA_real_, se = NA_real_))
  }

  n_pos <- as.integer(n_pos)
  n_neg <- as.integer(n_neg)
  Q1 <- auc / (2 - auc)
  Q2 <- 2 * auc^2 / (1 + auc)
  var_auc <- (auc * (1 - auc) +
                (n_pos - 1L) * (Q1 - auc^2) +
                (n_neg - 1L) * (Q2 - auc^2)) /
    (n_pos * n_neg)
  if (!is.finite(var_auc) || var_auc < 0) {
    return(c(ci_lo = NA_real_, ci_hi = NA_real_, se = NA_real_))
  }
  se_auc <- sqrt(var_auc)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  c(
    ci_lo = max(0, auc - z * se_auc),
    ci_hi = min(1, auc + z * se_auc),
    se = se_auc
  )
}


# ----------------------------------------------------------------------------
# BH-FDR and Holm helpers
# ----------------------------------------------------------------------------

#' @title BH-FDR correction within a matched-null modality family
#'
#' @description
#' Applies Benjamini-Hochberg FDR correction within a single modality family
#' of matched-null tests. Per plan v0.3.3 section 5: matched-null tests use
#' BH-FDR within each of the 4 modality families separately (not across modalities).
#'
#' @param results List of \code{\link{paper3_matched_null_benchmark}} return
#'   objects from a single modality family.
#'
#' @return Numeric vector of BH-FDR q-values, same length as \code{results}.
#'
#' @export
paper3_bh_fdr_correct_matched_null <- function(results) {
  p_emp <- vapply(results, function(r) r$p_emp, numeric(1L))
  stats::p.adjust(p_emp, method = "BH")
}


#' @title Holm correction for family-wise method contrasts
#'
#' @description
#' Applies Holm step-down correction across a family-wise contrast set
#' (method \eqn{\times} claim-state contrasts; hemolysis-injection deltas).
#' Per plan v0.3.3 section 5: Holm is reserved for these family-wise comparisons,
#' NOT for discovery-screening matched-null tests (which use BH-FDR).
#'
#' @param p_values Numeric vector of contrast p-values.
#'
#' @return Numeric vector of Holm-corrected p-values.
#'
#' @export
paper3_holm_correct_familywise <- function(p_values) {
  stats::p.adjust(p_values, method = "holm")
}


#' @title Block-aware two-stage BH-FDR for specimen-shared cohort clusters
#'
#' @description
#' Naïve BH-FDR within a modality family treats specimen-shared cohorts (e.g.,
#' Toray clusters with \eqn{\sim 94.5\%} specimen reuse) as independent tests,
#' which is anti-conservative due to strongly correlated test statistics. This
#' function implements a two-stage correction:
#' \enumerate{
#'   \item Within each provenance block, combine cohort-level p-values via
#'     Fisher's chi-square statistic. Two reference distributions are
#'     reported side by side:
#'     \itemize{
#'       \item \emph{conservative reference}: \eqn{S/c} against
#'         \eqn{\chi^2_{2k}}; intentionally conservative; used for the
#'         manuscript's headline counts.
#'       \item \emph{textbook reference} (Brown 1975 / Kost-McDermott 2002):
#'         \eqn{S/c} against \eqn{\chi^2_{2k/c}} with \eqn{df} rescaled by
#'         the same inflation factor.
#'     }
#'   \item Apply BH-FDR across block-level p-values, separately for each
#'     reference.
#' }
#'
#' The inflation factor is \eqn{c(k) = 1 + (k - 1)\rho} with \eqn{k} the
#' number of cohorts in the block and \eqn{\rho} the assumed within-block
#' correlation (default \eqn{\rho = 0.25} reproduces the manuscript's headline).
#'
#' @param results List of \code{\link{paper3_matched_null_benchmark}} return
#'   objects.
#' @param block_id Character vector --- provenance block identifier per
#'   cohort (e.g., \code{"Toray-cluster"}, \code{"Affy-singleton"}).
#' @param rho Numeric in [0, 1]. Assumed within-block correlation; controls
#'   the inflation factor \eqn{c(k) = 1 + (k - 1)\rho}. Default 0.25
#'   reproduces the manuscript headline. Set to 0 to recover vanilla Fisher.
#'
#' @return A \code{data.frame} with columns: \code{cohort_idx},
#'   \code{block_id}, \code{p_emp}, \code{p_block} (conservative),
#'   \code{p_block_textbook}, \code{q_block_BH} (conservative),
#'   \code{q_block_BH_textbook}, \code{q_cohort_within_block}.
#'
#' @export
paper3_bh_fdr_correct_blocked <- function(results, block_id, rho = 0.25) {
  stopifnot(length(results) == length(block_id),
            is.numeric(rho), length(rho) == 1L,
            !is.na(rho), rho >= 0, rho <= 1)
  p_emp <- vapply(results, function(r) r$p_emp, numeric(1L))
  blocks <- split(seq_along(p_emp), block_id)

  # Conservative reference: S/c against chi-square(2k).
  block_p_cons <- vapply(names(blocks), function(b) {
    idx <- blocks[[b]]
    if (length(idx) == 1L) return(p_emp[idx])
    chisq_stat <- -2 * sum(log(p_emp[idx]))
    df <- 2 * length(idx)
    bh_factor <- 1 + (length(idx) - 1L) * rho
    stats::pchisq(chisq_stat / bh_factor, df = df, lower.tail = FALSE)
  }, numeric(1L))

  # Textbook Brown 1975 / Kost-McDermott 2002 reference:
  # S/c against chi-square(2k/c).
  block_p_text <- vapply(names(blocks), function(b) {
    idx <- blocks[[b]]
    if (length(idx) == 1L) return(p_emp[idx])
    chisq_stat <- -2 * sum(log(p_emp[idx]))
    df_raw <- 2 * length(idx)
    bh_factor <- 1 + (length(idx) - 1L) * rho
    df_textbook <- df_raw / bh_factor
    stats::pchisq(chisq_stat / bh_factor, df = df_textbook, lower.tail = FALSE)
  }, numeric(1L))

  block_q          <- stats::p.adjust(block_p_cons, method = "BH")
  block_q_textbook <- stats::p.adjust(block_p_text, method = "BH")

  cohort_q_within <- numeric(length(p_emp))
  for (b in names(blocks)) {
    idx <- blocks[[b]]
    cohort_q_within[idx] <- stats::p.adjust(p_emp[idx], method = "BH")
  }

  data.frame(
    cohort_idx           = seq_along(p_emp),
    block_id             = block_id,
    p_emp                = p_emp,
    p_block              = block_p_cons[match(block_id, names(block_p_cons))],
    p_block_textbook     = block_p_text[match(block_id, names(block_p_text))],
    q_block_BH           = block_q[match(block_id, names(block_q))],
    q_block_BH_textbook  = block_q_textbook[match(block_id, names(block_q_textbook))],
    q_cohort_within_block = cohort_q_within
  )
}


# ----------------------------------------------------------------------------
# 5-fold nested-CV matched-null benchmark (Patch 4)
# ----------------------------------------------------------------------------

#' @noRd
.paper3_null_coalesce <- function(a, b) if (is.null(a)) b else a

#' Stratified k-folds (ungrouped). Returns list of k test-index vectors.
#' @noRd
.paper3_stratified_folds <- function(y, k = 5L, seed = 42L) {
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
#' @noRd
.paper3_grouped_folds <- function(y, group_id, k = 5L, seed = 42L) {
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
#' @noRd
.paper3_select_panel_training_fold <- function(X_train, y_train, panel_size = 20L) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop(".paper3_select_panel_training_fold requires pROC.")
  }
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

#' Internal variant of .paper3_draw_matched_panel that also propagates the
#' fallback tier (1=exact, 2=oversampled, 3=degenerate). Mirrors the
#' auxiliary script at OmicSelector_paper/code/methods/matched_null_benchmark.R.
#' @noRd
.paper3_draw_matched_panel_with_tier <- function(strata_info, exclude = character(0)) {
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
  if (shortfall > 0L) {
    extra_pool <- setdiff(names(strata), c(exclude, picked))
    if (length(extra_pool) >= shortfall) {
      picked <- c(picked, sample(extra_pool, shortfall))
    } else {
      out <- character(0)
      attr(out, "fallback_tier") <- 3L
      return(out)
    }
  }
  attr(picked, "fallback_tier") <- fallback_tier
  picked
}


#' @title 5-fold nested-CV matched-null benchmark
#'
#' @description
#' Nested cross-validation variant of \code{\link{paper3_matched_null_benchmark}}.
#' Outer folds are stratified by outcome; when \code{group_id} is supplied
#' (e.g., matched-set patient IDs from a provenance audit), grouped folds
#' guarantee no patient group is split across train/test. For each outer fold
#' a panel of size \code{panel_size} is selected on the TRAINING fold by
#' univariate AUC, and the matched-null strata are computed on the training
#' fold; the held-out test fold is then used to compute both the observed
#' AUC and \code{K} matched-null AUCs.
#'
#' This is the engine that produced the manuscript's per-cell numbers in
#' Supplementary Table S1 / Figure 3 / Table 2.
#'
#' @param X Numeric matrix (samples \eqn{\times} features) with column names.
#' @param y Binary outcome (0/1 or factor with two levels).
#' @param panel_size Integer. Panel size selected per outer fold on training
#'   data. Default 20.
#' @param scoring_fn Function \code{function(X_panel)} returning per-sample
#'   scores; if \code{NULL} the default \code{rowSums(log(X + 0.5))} is used.
#'   The closure must NEVER look at held-out data.
#' @param K Integer. Null replicates per fold. Default 10000.
#' @param outer_k Integer. Outer CV folds. Default 5.
#' @param group_id Vector of length \code{nrow(X)}, or \code{NULL}. When
#'   provided, grouped K-fold is used and groups are never split.
#' @param feature_pool Character. Full candidate pool. \code{NULL} uses
#'   \code{colnames(X)}.
#' @param exclude_features Character. Features excluded from the null
#'   candidate pool.
#' @param post_hemolysis_corrected Logical. If \code{TRUE} prepends the
#'   canonical hemolysis-marker list to \code{exclude_features}.
#' @param seed Integer in [0, 21474]. Fold-assignment and per-draw null
#'   seeds are derived deterministically from this value.
#' @param min_valid_folds Integer. Default 4 (manuscript / plan v0.6 sec 4.1).
#' @param min_pos_per_fold Integer. Default 5.
#' @param min_neg_per_fold Integer. Default 5.
#' @param null_parallel_cores Integer. Number of forked workers per fold.
#'   Default 1 (serial).
#' @param auc_ci_method Character; one of \code{"hanley_mcneil"} (default),
#'   \code{"delong"} (uses pooled fold-level scores via
#'   \code{pROC::ci.auc(method = "delong")}), or \code{"none"}. Returned as
#'   \code{auc_obs_ci_lo} / \code{auc_obs_ci_hi} on the cohort-level
#'   pooled-fold prediction.
#'
#' @return Named list including \code{auc_obs_cv}, \code{auc_obs_ci_lo},
#'   \code{auc_obs_ci_hi}, \code{p_emp_cv}, \code{n_valid_folds},
#'   \code{eligible}, \code{fold_panels}, and a per-fold audit. See the
#'   manuscript Methods §"Statistical evaluation".
#'
#' @examples
#' \dontrun{
#' set.seed(99)
#' X <- matrix(rlnorm(150 * 120, 5, 1), nrow = 150,
#'             dimnames = list(NULL, paste0("f-", seq_len(120))))
#' y <- rbinom(150, 1, 0.5)
#' X[y == 1, paste0("f-", 1:10)] <- X[y == 1, paste0("f-", 1:10)] * 3
#' res <- paper3_matched_null_benchmark_cv(X, y, panel_size = 10L,
#'                                          K = 100L, outer_k = 5L)
#' res$p_emp_cv
#' }
#'
#' @export
paper3_matched_null_benchmark_cv <- function(
    X,
    y,
    panel_size               = 20L,
    scoring_fn               = NULL,
    K                        = 10000L,
    outer_k                  = 5L,
    group_id                 = NULL,
    feature_pool             = NULL,
    exclude_features         = character(0),
    post_hemolysis_corrected = FALSE,
    seed                     = 42L,
    min_valid_folds          = 4L,
    min_pos_per_fold         = 5L,
    min_neg_per_fold         = 5L,
    null_parallel_cores      = 1L,
    auc_ci_method            = c("hanley_mcneil", "delong", "none")) {

  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("paper3_matched_null_benchmark_cv requires pROC.")
  }
  auc_ci_method <- match.arg(auc_ci_method)
  stopifnot(is.matrix(X) || is.data.frame(X))
  X <- as.matrix(X)
  if (is.factor(y)) y <- as.integer(y) - 1L
  if (is.null(feature_pool)) feature_pool <- colnames(X)
  seed <- suppressWarnings(as.integer(seed))
  if (!is.finite(seed) || is.na(seed) || seed < 0L || seed > 21474L) {
    stop("paper3_matched_null_benchmark_cv: seed must be an integer from 0 to 21474.")
  }
  null_parallel_cores <- suppressWarnings(as.integer(null_parallel_cores))
  if (!is.finite(null_parallel_cores) || null_parallel_cores < 1L) {
    null_parallel_cores <- 1L
  }
  rng_protocol <- "v0.6_per_draw_seed_schedule_independent"

  has_groups      <- !is.null(group_id) && length(group_id) == nrow(X)
  grouping_policy <- if (has_groups) "groupKFold" else "ungrouped_stratified"
  n_groups        <- if (has_groups) length(unique(group_id)) else NA_integer_

  hemolysis_markers <- c("hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p",
                          "hsa-miR-144-3p", "hsa-miR-223-3p")
  if (post_hemolysis_corrected) {
    exclude_features <- unique(c(exclude_features, hemolysis_markers))
  }
  feature_pool <- setdiff(feature_pool, exclude_features)

  folds <- if (has_groups)
    .paper3_grouped_folds(y, group_id, k = outer_k, seed = seed)
  else
    .paper3_stratified_folds(y, k = outer_k, seed = seed)

  n_group_split_violations <- 0L
  if (has_groups) {
    for (fi in seq_along(folds)) {
      train_idx <- setdiff(seq_len(nrow(X)), folds[[fi]])
      viol <- length(intersect(unique(group_id[folds[[fi]]]),
                                unique(group_id[train_idx])))
      if (viol > 0L)
        stop(sprintf(
          "paper3_matched_null_benchmark_cv: %d group(s) split across train/test in fold %d.",
          viol, fi))
    }
  }

  auc_obs_folds   <- rep(NA_real_, outer_k)
  auc_null_folds  <- matrix(NA_real_, nrow = K, ncol = outer_k)
  fold_panels     <- vector("list", outer_k)
  fold_audit      <- vector("list", outer_k)
  fold_test_idx   <- vector("list", outer_k)
  fold_test_score <- vector("list", outer_k)
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

    train_pool <- intersect(feature_pool, colnames(X_train))
    if (length(train_pool) < 2L) {
      fa$skip_reason <- "train_pool < 2 features"; fold_audit[[fi]] <- fa; next
    }
    X_train_pool <- X_train[, train_pool, drop = FALSE]
    panel_f <- tryCatch(
      .paper3_select_panel_training_fold(X_train_pool, y_train, panel_size),
      error = function(e) NULL)
    if (is.null(panel_f) || length(panel_f) == 0L) {
      fa$skip_reason <- "panel selection failed"; fold_audit[[fi]] <- fa; next
    }
    fold_panels[[fi]] <- panel_f

    dr_train <- colMeans(X_train_pool > 0,  na.rm = TRUE)
    ma_train <- colMeans(X_train_pool,       na.rm = TRUE)
    strata_info <- .paper3_match_panel_strata(
      panel = panel_f, feature_pool = train_pool,
      detection_rate = dr_train, mean_abundance = ma_train)

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
    fold_test_idx[[fi]]   <- test_idx
    fold_test_score[[fi]] <- score_obs_f

    .run_null_draw <- function(k) {
      set.seed(seed * 100000L + fi * 10000L + k)
      p_k <- .paper3_draw_matched_panel_with_tier(strata_info, exclude = panel_f)
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
      stop(sprintf("paper3_matched_null_benchmark_cv: null draw failed in fold %d.", fi))
    }

    auc_null_folds[, fi] <- vapply(null_results, `[[`, numeric(1L), "auc")
    fold_fallback_count <- sum(vapply(null_results, `[[`, integer(1L), "fallback_empty"))
    tiers <- vapply(null_results, `[[`, integer(1L), "tier")
    tier_counts <- c(`1` = sum(tiers == 1L), `2` = sum(tiers == 2L), `3` = sum(tiers == 3L))

    fa$fold_fallback_count       <- fold_fallback_count
    fa$fold_fallback_tier3_pct   <- 100 * fold_fallback_count / max(1L, K)
    fa$fold_fallback_tier_counts <- as.list(tier_counts)
    fold_audit[[fi]] <- fa
  }

  n_valid_folds <- sum(valid_fold_mask)
  valid_fi      <- which(valid_fold_mask)

  pos_v <- vapply(valid_fi, function(fi) fold_audit[[fi]]$n_pos_test, integer(1L))
  neg_v <- vapply(valid_fi, function(fi) fold_audit[[fi]]$n_neg_test, integer(1L))
  min_pos_obs <- if (length(pos_v) > 0L) min(pos_v) else NA_integer_
  min_neg_obs <- if (length(neg_v) > 0L) min(neg_v) else NA_integer_

  skip_summary <- paste(vapply(fold_audit, function(fa)
    .paper3_null_coalesce(fa$skip_reason, "ok"), character(1L)), collapse = "; ")

  if (n_valid_folds < min_valid_folds) {
    return(list(
      eligible                  = FALSE,
      ineligible_reason         = sprintf(
        "n_valid_folds=%d < min_valid_folds=%d; per-fold: %s",
        n_valid_folds, min_valid_folds, skip_summary),
      n_valid_folds             = n_valid_folds,
      grouping_policy           = grouping_policy,
      n_groups                  = n_groups,
      n_group_split_violations  = n_group_split_violations,
      matched_set_populated_pct = NA_real_,
      min_pos_per_fold          = min_pos_obs,
      min_neg_per_fold          = min_neg_obs,
      auc_obs_ci_lo             = NA_real_,
      auc_obs_ci_hi             = NA_real_,
      auc_ci_method             = auc_ci_method,
      fold_audit                = fold_audit,
      K = K, outer_k = outer_k, seed = seed,
      null_parallel_cores       = null_parallel_cores,
      rng_protocol              = rng_protocol))
  }

  auc_obs_cv  <- mean(auc_obs_folds[valid_fold_mask], na.rm = TRUE)
  null_valid  <- auc_null_folds[, valid_fold_mask, drop = FALSE]
  complete_null <- rowSums(is.na(null_valid)) == 0L
  auc_null_cv <- if (any(complete_null))
    rowMeans(null_valid[complete_null, , drop = FALSE], na.rm = FALSE)
  else
    numeric(0L)
  n_null_valid <- sum(!is.na(auc_null_cv))
  p_emp_cv    <- (1 + sum(auc_null_cv >= auc_obs_cv, na.rm = TRUE)) /
                 (1 + n_null_valid)

  total_fallback <- sum(vapply(fold_audit[valid_fi], function(fa)
    .paper3_null_coalesce(fa$fold_fallback_count, 0L), integer(1L)))
  max_fallback_tier <- max(vapply(fold_audit[valid_fi], function(fa) {
    tc <- fa$fold_fallback_tier_counts %||% list(`1` = 0L, `2` = 0L, `3` = 0L)
    used <- as.integer(names(tc)[unlist(tc) > 0L])
    if (length(used) == 0L) 1L else max(used)
  }, integer(1L)))

  # Pooled-fold per-cell AUC CI. Concatenate held-out scores in fold order.
  pooled_y      <- unlist(lapply(valid_fi,
                                   function(fi) y[fold_test_idx[[fi]]]))
  pooled_scores <- unlist(lapply(valid_fi,
                                   function(fi) fold_test_score[[fi]]))
  if (auc_ci_method == "delong") {
    auc_ci <- .paper3_auc_ci(pooled_y, pooled_scores, auc_obs_cv, "delong")
  } else if (auc_ci_method == "hanley_mcneil") {
    auc_ci <- .paper3_auc_ci(pooled_y, pooled_scores, auc_obs_cv, "hanley_mcneil")
  } else {
    auc_ci <- c(NA_real_, NA_real_)
  }
  auc_obs_ci_lo <- auc_ci[1L]
  auc_obs_ci_hi <- auc_ci[2L]

  list(
    auc_obs_cv               = auc_obs_cv,
    auc_obs_ci_lo            = auc_obs_ci_lo,
    auc_obs_ci_hi            = auc_obs_ci_hi,
    auc_ci_method            = auc_ci_method,
    auc_obs_folds            = auc_obs_folds,
    auc_null_cv              = auc_null_cv,
    p_emp_cv                 = p_emp_cv,
    n_valid_folds            = n_valid_folds,
    fold_panels              = fold_panels,
    K                        = K,
    outer_k                  = outer_k,
    seed                     = seed,
    null_parallel_cores      = null_parallel_cores,
    rng_protocol             = rng_protocol,
    eligible                 = TRUE,
    ineligible_reason        = NA_character_,
    fallback_count           = total_fallback,
    fallback_tier            = max_fallback_tier,
    grouping_policy          = grouping_policy,
    n_groups                 = n_groups,
    n_group_split_violations = n_group_split_violations,
    matched_set_populated_pct = NA_real_,
    min_pos_per_fold         = min_pos_obs,
    min_neg_per_fold         = min_neg_obs,
    fold_audit               = fold_audit,
    panel_size               = panel_size,
    exclude_features         = exclude_features,
    post_hemolysis_corrected = post_hemolysis_corrected
  )
}
