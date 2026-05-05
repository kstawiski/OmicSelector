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
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{auc_obs}}{Observed panel AUC.}
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
                                           seed = 42L) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("paper3_matched_null_benchmark requires pROC for AUC computation.")
  }
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
# BH-FDR and Holm helpers
# ----------------------------------------------------------------------------

#' @title BH-FDR correction within a matched-null modality family
#'
#' @description
#' Applies Benjamini-Hochberg FDR correction within a single modality family
#' of matched-null tests. Per plan v0.3.3 \S5: matched-null tests use BH-FDR
#' within each of the 4 modality families separately (not across modalities).
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
#' Per plan v0.3.3 \S5: Holm is reserved for these family-wise comparisons,
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
#'     Fisher's chi-square statistic with a conservative inflation factor for
#'     positive dependence (Bauer-Hommel heuristic).
#'   \item Apply BH-FDR across block-level p-values.
#' }
#'
#' @param results List of \code{\link{paper3_matched_null_benchmark}} return
#'   objects.
#' @param block_id Character vector — provenance block identifier per cohort
#'   (e.g., \code{"Toray-cluster"}, \code{"Affy-singleton"}).
#'
#' @return A \code{data.frame} with columns: \code{cohort_idx}, \code{block_id},
#'   \code{p_emp}, \code{p_block}, \code{q_block_BH}, \code{q_cohort_within_block}.
#'
#' @export
paper3_bh_fdr_correct_blocked <- function(results, block_id) {
  stopifnot(length(results) == length(block_id))
  p_emp <- vapply(results, function(r) r$p_emp, numeric(1L))
  blocks <- split(seq_along(p_emp), block_id)

  block_p <- vapply(names(blocks), function(b) {
    idx <- blocks[[b]]
    if (length(idx) == 1L) return(p_emp[idx])
    chisq_stat <- -2 * sum(log(p_emp[idx]))
    df <- 2 * length(idx)
    bh_factor <- 1 + (length(idx) - 1L) / 4
    stats::pchisq(chisq_stat / bh_factor, df = df, lower.tail = FALSE)
  }, numeric(1L))

  block_q <- stats::p.adjust(block_p, method = "BH")

  cohort_q_within <- numeric(length(p_emp))
  for (b in names(blocks)) {
    idx <- blocks[[b]]
    cohort_q_within[idx] <- stats::p.adjust(p_emp[idx], method = "BH")
  }

  data.frame(
    cohort_idx = seq_along(p_emp),
    block_id = block_id,
    p_emp = p_emp,
    p_block = block_p[match(block_id, names(block_p))],
    q_block_BH = block_q[match(block_id, names(block_q))],
    q_cohort_within_block = cohort_q_within
  )
}
