#' @title Single-sample within-sample compositional methods (Module A, P1)
#'
#' @description
#' Five Module-A within-sample compositional transforms introduced in the
#' OmicSelector single-sample scoring bank. All methods are single-sample and
#' reference-cohort-free:
#' each sample is normalised, scored, or projected using only its own panel
#' values, with no requirement for an external training cohort to exist at
#' deployment time. This is the property that makes within-sample methods
#' suitable for clinical translation: a diagnostic test built on these
#' primitives can be deployed by computing a per-patient score from the
#' patient's own panel measurement, without re-fitting any normalization on
#' a population reference.
#'
#' Methods provided:
#' \itemize{
#'   \item \code{\link{ws_rclr_trimmed}}: robust trimmed centred log-ratio
#'     with optional zero-detection-only convention; tolerant to dominant
#'     contaminating species (haemolysis miRNAs) and to small-head-of-pool
#'     compositions typical of plasma/serum miRNA panels.
#'   \item \code{\link{ws_balance_ilr}}: orthonormal isometric log-ratio
#'     balances on a frozen sequential binary partition encoding circulating
#'     miRNA biology (red-blood-cell vs rest, platelet vs rest, let-7 family
#'     within-cluster, miR-17~92, miR-200, miR-371~373, dominant-vs-tail
#'     abundance balance).
#'   \item \code{\link{ws_alr_pivot}}: additive log-ratio with a frozen
#'     pivot pool of putatively low-variance circulating miRNAs serving as
#'     the geometric-mean denominator (anchors the transformation against
#'     a stable internal reference rather than choosing a single feature
#'     as denominator).
#'   \item \code{\link{ws_mad_logratio}}: median-centred log-ratio with
#'     optional per-sample MAD scaling; not Aitchison-valid but a robust
#'     descriptive transform for downstream learners that do not require
#'     zero-sum compositional coordinates.
#'   \item \code{\link{ws_dominance_score}}: panel-internal scalar QC
#'     dominance score reporting the fraction of total panel signal carried
#'     by the top-K features; flags pre-analytical failures (haemolysis,
#'     incomplete extraction) when calibrated against a Tier R reference
#'     distribution.
#' }
#'
#' Helper utilities are provided alongside: \code{\link{ws_default_sbp}}
#' returns the curated v1 sequential-binary-partition tree for circulating
#' miRNA panels, and \code{\link{ws_default_pivot_pool}} returns the curated
#' v1 ALR pivot pool (six low-variance miRNAs; miR-92a-3p excluded because
#' recent serum/plasma literature treats it as erythrocyte-derived).
#'
#' @details
#' Background: circulating miRNA biomarkers were proposed as minimally
#' invasive cancer markers in 2008 (Mitchell et al., PNAS). In the
#' intervening years no clinically actionable cancer-screening test has
#' emerged from this source despite thousands of "candidate panel" papers.
#' The reasons are predominantly technical (batch effects, hemolysis
#' contamination, lack of standards, reference-cohort dependency, cross-cohort
#' leakage, absence of provenance/claim-gating) rather than biological.
#' These five within-sample compositional methods target two of those failure
#' modes directly: reference-cohort dependency (each method scores a single
#' sample using only its own panel values) and contamination by dominant
#' erythrocyte-derived species (miR-451a, miR-486-5p, miR-16, miR-144-3p,
#' miR-92a-3p - all explicitly handled by either exclusion from the rCLR
#' centering subset or by having dedicated balances in the SBP tree).
#'
#' Module A is one of six modules in the OmicSelector v2.4 toolkit. The
#' other modules - Module B (hemolysis correction), Module C (frozen batch
#' correction), Module D (outlier detection and conformal claim-gating),
#' Module E (multi-cancer / domain generalisation), and Module F
#' (self-supervised pretraining and ablation) - are documented in the
#' single-sample scoring bank and the per-module help pages.
#'
#' @references
#' Aitchison J. (1986) The Statistical Analysis of Compositional Data. Chapman & Hall.
#'
#' Egozcue J. J., Pawlowsky-Glahn V. (2005) Groups of parts and their balances
#' in compositional data analysis. \emph{Mathematical Geology} 37(7): 795-828.
#'
#' Vandeputte D., Kathagen G., D'hoe K., et al. (2017) Quantitative microbiome
#' profiling links gut community variation to microbial load. \emph{Nature}
#' 551(7681): 507-511. (Origin of the "detected_only" rCLR convention.)
#'
#' Mitchell P. S., Parkin R. K., Kroh E. M., et al. (2008) Circulating
#' microRNAs as stable blood-based markers for cancer detection.
#' \emph{Proceedings of the National Academy of Sciences USA} 105(30):
#' 10513-10518.
#'
#' @name singlesample-within-sample
NULL


# ----------------------------------------------------------------------------
# ws_rclr_trimmed
# ----------------------------------------------------------------------------

#' @title Robust trimmed centred log-ratio for compositional miRNA panels
#'
#' @description
#' Computes a robust trimmed centred log-ratio (rCLR) on a single sample or a
#' samples x features matrix. The centering subset excludes a configurable
#' set of contaminating miRNAs (default: the haemolysis and platelet-activation
#' panel) and trims the top-\code{trim_upper} and bottom-\code{trim_lower}
#' fractions of the remaining features by abundance before computing the
#' geometric mean used as the CLR denominator. Falls back to a global CLR
#' (with warning) when the centering subset shrinks below
#' \code{min_centering_size}.
#'
#' @param x Numeric vector or matrix (samples x features). For a matrix the
#'   transformation is applied row-wise.
#' @param pseudocount Additive pseudocount to avoid log(0). Default 1e-6 of
#'   the sample sum.
#' @param trim_upper Fraction of highest-abundance features to drop from the
#'   centering subset. Default 0.10.
#' @param trim_lower Fraction of lowest-abundance features to drop from the
#'   centering subset. Default 0.05.
#' @param exclude_features Named vector or character vector of features to
#'   always drop from the centering subset (e.g., known haemolysis markers).
#'   Default = c("hsa-miR-451a","hsa-miR-16-5p","hsa-miR-486-5p","hsa-miR-144-3p","hsa-miR-223-3p").
#' @param zero_policy "pseudocount" (default) replaces zeros with the
#'   pseudocount; "detected_only" excludes zeros from the centering subset
#'   (rCLR convention; Vandeputte et al. 2017).
#' @param min_centering_size Minimum number of features that must remain after
#'   trim + exclusion to compute the rCLR. If fewer remain, falls back to
#'   global CLR with a warning. Default 8.
#'
#' @return Same shape as \code{x}; the rCLR-transformed values. When invoked
#'   on a single sample the return value carries attributes
#'   \code{centering_features}, \code{n_centering}, and \code{pseudocount}
#'   for downstream auditability.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x <- rlnorm(50, meanlog = 5, sdlog = 1.2)
#' names(x) <- c(paste0("hsa-miR-", sprintf("%03d", seq_len(45))),
#'               "hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p",
#'               "hsa-miR-144-3p", "hsa-miR-223-3p")
#' z <- ws_rclr_trimmed(x)
#' attr(z, "n_centering")  # number of features used for centering
#' }
#'
#' @export
ws_rclr_trimmed <- function(x,
                             pseudocount = NULL,
                             trim_upper = 0.10,
                             trim_lower = 0.05,
                             exclude_features = c("hsa-miR-451a", "hsa-miR-16-5p",
                                                   "hsa-miR-486-5p", "hsa-miR-144-3p",
                                                   "hsa-miR-223-3p"),
                             zero_policy = c("pseudocount", "detected_only"),
                             min_centering_size = 8L) {
  zero_policy <- match.arg(zero_policy)

  if (is.matrix(x) || is.data.frame(x)) {
    out <- t(apply(x, 1L, function(row) {
      ws_rclr_trimmed(row,
                      pseudocount = pseudocount,
                      trim_upper = trim_upper,
                      trim_lower = trim_lower,
                      exclude_features = exclude_features,
                      zero_policy = zero_policy,
                      min_centering_size = min_centering_size)
    }))
    dimnames(out) <- dimnames(x)
    return(out)
  }

  if (length(x) < min_centering_size) {
    stop("ws_rclr_trimmed: panel size (", length(x), ") < min_centering_size (",
         min_centering_size, "). Refusing to compute on degenerate input.")
  }

  feat_names <- names(x)
  if (is.null(feat_names)) feat_names <- as.character(seq_along(x))

  if (is.null(pseudocount)) pseudocount <- 1e-6 * sum(x, na.rm = TRUE)
  if (pseudocount <= 0) pseudocount <- .Machine$double.eps

  x_pos <- x + pseudocount
  total <- sum(x_pos, na.rm = TRUE)
  p     <- x_pos / total
  log_p <- log(p)

  is_excluded <- feat_names %in% exclude_features
  if (zero_policy == "detected_only") {
    is_zero <- x <= 0 | is.na(x)
  } else {
    is_zero <- rep(FALSE, length(x))
  }
  candidates <- !is_excluded & !is_zero

  ranks <- rank(-log_p, ties.method = "first")
  n_cand <- sum(candidates)
  if (n_cand < min_centering_size) {
    warning("ws_rclr_trimmed: candidate centering set (", n_cand,
            ") < min_centering_size; falling back to global CLR.")
    return(log_p - mean(log_p, na.rm = TRUE))
  }

  upper_cut <- floor(trim_upper * n_cand)
  lower_cut <- floor(trim_lower * n_cand)
  cand_ranks_sorted <- sort(ranks[candidates])
  keep_rank_lo <- cand_ranks_sorted[upper_cut + 1L]
  keep_rank_hi <- cand_ranks_sorted[n_cand - lower_cut]
  in_centering <- candidates & ranks >= keep_rank_lo & ranks <= keep_rank_hi

  if (sum(in_centering) < min_centering_size) {
    warning("ws_rclr_trimmed: centering subset shrank below threshold; falling back to global CLR.")
    return(log_p - mean(log_p, na.rm = TRUE))
  }

  centering <- mean(log_p[in_centering], na.rm = TRUE)
  z <- log_p - centering
  names(z) <- feat_names
  attr(z, "centering_features") <- feat_names[in_centering]
  attr(z, "n_centering") <- sum(in_centering)
  attr(z, "pseudocount") <- pseudocount
  z
}


# ----------------------------------------------------------------------------
# fit_ws_rclr_panel / score_ws_rclr_panel
# ----------------------------------------------------------------------------

#' Fit a training-frozen trimmed-rCLR signed-panel score
#'
#' @description
#' Fits the canonical trimmed-rCLR reference used by the OmicSelector
#' single-sample benchmark. Each training specimen is transformed independently
#' with [ws_rclr_trimmed()]. Features are ranked by the absolute standardized
#' case-control mean difference in the training data, the first `panel_size`
#' features are retained, and each feature receives the sign of its training
#' case-minus-control median difference. No held-out value or outcome is used.
#'
#' @param X_train Numeric training matrix with specimens in rows and named
#'   features in columns.
#' @param y_train Binary training labels coded 0/1.
#' @param meta_train Unused; accepted for the canonical roster fit contract.
#' @param hp Optional named list. Supported fields are `panel_size` (default
#'   20), `pseudocount`, `trim_upper` (0.10), `trim_lower` (0.05),
#'   `exclude_features`, `zero_policy` (`"pseudocount"`), and
#'   `min_centering_size` (8).
#'
#' @return A `ws_rclr_panel_model` containing the frozen feature panel, signs,
#'   training effects, and transform settings.
#' @seealso [score_ws_rclr_panel()], [ws_rclr_trimmed()]
#' @export
fit_ws_rclr_panel <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  if (is.data.frame(X_train)) X_train <- as.matrix(X_train)
  if (!is.matrix(X_train) || !is.numeric(X_train) || nrow(X_train) < 4L ||
      ncol(X_train) < 2L) {
    stop("fit_ws_rclr_panel: X_train must be a numeric samples-by-features matrix with at least four rows and two columns.",
         call. = FALSE)
  }
  if (is.null(colnames(X_train)) || any(!nzchar(colnames(X_train))) ||
      anyDuplicated(colnames(X_train))) {
    stop("fit_ws_rclr_panel: X_train requires unique non-empty feature names.",
         call. = FALSE)
  }
  y <- as.integer(y_train)
  if (length(y) != nrow(X_train) || anyNA(y) || !setequal(unique(y), 0:1)) {
    stop("fit_ws_rclr_panel: y_train must contain aligned binary 0/1 labels with both classes present.",
         call. = FALSE)
  }
  if (!is.null(meta_train) && nrow(as.data.frame(meta_train)) != nrow(X_train)) {
    stop("fit_ws_rclr_panel: meta_train must have one row per training specimen.",
         call. = FALSE)
  }
  if (!is.list(hp)) stop("fit_ws_rclr_panel: hp must be a named list.", call. = FALSE)
  if (length(hp) && (is.null(names(hp)) || any(!nzchar(names(hp))))) {
    stop("fit_ws_rclr_panel: all hp entries must be named.", call. = FALSE)
  }
  defaults <- list(
    panel_size = 20L,
    pseudocount = NULL,
    trim_upper = 0.10,
    trim_lower = 0.05,
    exclude_features = c("hsa-miR-451a", "hsa-miR-16-5p",
                         "hsa-miR-486-5p", "hsa-miR-144-3p",
                         "hsa-miR-223-3p"),
    zero_policy = "pseudocount",
    min_centering_size = 8L
  )
  unknown <- setdiff(names(hp), names(defaults))
  if (length(unknown)) {
    stop("fit_ws_rclr_panel: unknown hp field(s): ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  settings <- utils::modifyList(defaults, hp)
  settings$panel_size <- as.integer(settings$panel_size)
  settings$min_centering_size <- as.integer(settings$min_centering_size)
  if (length(settings$panel_size) != 1L || is.na(settings$panel_size) ||
      settings$panel_size < 1L) {
    stop("fit_ws_rclr_panel: hp$panel_size must be a positive integer.",
         call. = FALSE)
  }
  transform_args <- settings[setdiff(names(settings), "panel_size")]
  Z_train <- do.call(ws_rclr_trimmed, c(list(x = X_train), transform_args))
  effects <- vapply(colnames(Z_train), function(feature) {
    x <- Z_train[, feature]
    ok <- is.finite(x)
    denom <- stats::sd(x[ok])
    if (sum(ok) < 5L || !is.finite(denom) || denom <= 0) return(0)
    abs(mean(x[ok & y == 1L]) - mean(x[ok & y == 0L])) / denom
  }, numeric(1L))
  order_idx <- order(effects, decreasing = TRUE)
  panel <- colnames(Z_train)[order_idx][seq_len(min(settings$panel_size,
                                                     ncol(Z_train)))]
  signs <- vapply(panel, function(feature) {
    difference <- stats::median(Z_train[y == 1L, feature], na.rm = TRUE) -
      stats::median(Z_train[y == 0L, feature], na.rm = TRUE)
    if (is.finite(difference) && difference < 0) -1 else 1
  }, numeric(1L))
  model <- list(
    panel = panel,
    weights = signs,
    effects = effects[panel],
    panel_size = length(panel),
    transform_args = transform_args,
    feature_universe = colnames(X_train)
  )
  class(model) <- "ws_rclr_panel_model"
  model
}

#' Score specimens with a frozen trimmed-rCLR signed panel
#'
#' @description
#' Applies [ws_rclr_trimmed()] independently to each incoming specimen and sums
#' the transformed values over the training-frozen feature panel after applying
#' the training-frozen signs. Extra features are ignored. Missing panel features
#' fail closed rather than being imputed from a scored batch.
#'
#' @param model A `ws_rclr_panel_model` returned by [fit_ws_rclr_panel()].
#' @param X Numeric specimens-by-features matrix.
#' @param meta Unused; accepted for the canonical roster score contract.
#'
#' @return One finite numeric score per specimen.
#' @export
score_ws_rclr_panel <- function(model, X, meta = NULL) {
  if (!inherits(model, "ws_rclr_panel_model")) {
    stop("score_ws_rclr_panel: model must come from fit_ws_rclr_panel().",
         call. = FALSE)
  }
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X) || nrow(X) < 1L) {
    stop("score_ws_rclr_panel: X must be a numeric samples-by-features matrix.",
         call. = FALSE)
  }
  if (!is.null(meta) && nrow(as.data.frame(meta)) != nrow(X)) {
    stop("score_ws_rclr_panel: meta must have one row per specimen.",
         call. = FALSE)
  }
  missing <- setdiff(model$panel, colnames(X))
  if (length(missing)) {
    stop("score_ws_rclr_panel: missing frozen panel feature(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  Z <- do.call(ws_rclr_trimmed, c(list(x = X), model$transform_args))
  score <- as.numeric(Z[, model$panel, drop = FALSE] %*% model$weights)
  if (length(score) != nrow(X) || any(!is.finite(score))) {
    stop("score_ws_rclr_panel: non-finite score produced.", call. = FALSE)
  }
  score
}


# ----------------------------------------------------------------------------
# ws_balance_ilr - frozen sequential binary partition for ILR balances
# ----------------------------------------------------------------------------

#' @title Default circulating-miRNA sequential binary partition (v1)
#'
#' @description
#' Returns the v1 frozen sequential binary partition (SBP) for circulating
#' miRNA panels used by \code{\link{ws_balance_ilr}}. Each entry of the
#' returned list is a named partition with components \code{numerator} and
#' \code{denominator} listing miRNA features. Several partitions use
#' sentinel denominators \code{"ALL_NON_RBC"} and \code{"ALL_NON_PLT"} that
#' are resolved at runtime as "every panel feature not in the numerator,"
#' and one partition uses the sentinel \code{"TOP_K_BY_ABUNDANCE"} with a
#' \code{k} field, resolved at runtime by sorting the panel.
#'
#' The partition is fixed at v1 and any future revision will be tagged
#' \code{circulating_v2}, etc., with a separate justification log under
#' the single-sample scoring bank's sequential-binary-partition justification notes.
#'
#' @param version Character. Currently only \code{"circulating_v1"} is
#'   supported.
#'
#' @return Named list of balance partitions.
#' @export
ws_default_sbp <- function(version = "circulating_v1") {
  if (version != "circulating_v1") stop("Unknown SBP version: ", version)

  list(
    rbc_vs_rest = list(
      numerator   = c("hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p", "hsa-miR-144-3p"),
      denominator = "ALL_NON_RBC"
    ),
    platelet_vs_rest = list(
      numerator   = c("hsa-miR-223-3p", "hsa-miR-126-3p"),
      denominator = "ALL_NON_PLT"
    ),
    let7_a_vs_let7_g = list(
      numerator   = c("hsa-let-7a-5p", "hsa-let-7a-3p"),
      denominator = c("hsa-let-7g-5p", "hsa-let-7g-3p")
    ),
    let7_canonical_vs_b_d_f = list(
      numerator   = c("hsa-let-7a-5p","hsa-let-7c-5p","hsa-let-7e-5p","hsa-let-7i-5p"),
      denominator = c("hsa-let-7b-5p","hsa-let-7d-5p","hsa-let-7f-5p")
    ),
    mir17_vs_mir92 = list(
      numerator   = c("hsa-miR-17-5p","hsa-miR-18a-5p","hsa-miR-19a-3p","hsa-miR-19b-3p","hsa-miR-20a-5p"),
      denominator = c("hsa-miR-92a-3p")
    ),
    mir200_vs_mir141 = list(
      numerator   = c("hsa-miR-200a-3p","hsa-miR-200b-3p","hsa-miR-200c-3p","hsa-miR-429"),
      denominator = c("hsa-miR-141-3p","hsa-miR-141-5p")
    ),
    mir371_373_vs_rest_germcell = list(
      numerator   = c("hsa-miR-371a-3p","hsa-miR-372-3p","hsa-miR-373-3p"),
      denominator = c("hsa-miR-302a-3p","hsa-miR-302b-3p","hsa-miR-302c-3p","hsa-miR-302d-3p")
    ),
    dominant_vs_tail = list(
      numerator   = "TOP_K_BY_ABUNDANCE",
      denominator = "TAIL_BY_ABUNDANCE",
      k           = 8L
    )
  )
}


#' @title Isometric log-ratio balances on a frozen partition tree
#'
#' @description
#' Computes orthonormal isometric log-ratio (ILR) balances on a frozen
#' sequential binary partition. For each partition (G+, G-) the balance is
#' \deqn{b = \sqrt{\frac{|G+|\,|G-|}{|G+| + |G-|}} \log\!\left(\frac{\mathrm{gmean}(x_{G+})}{\mathrm{gmean}(x_{G-})}\right).}
#' The default partition (\code{\link{ws_default_sbp}}) encodes circulating
#' miRNA biology: erythrocyte-vs-rest, platelet-vs-rest, let-7 family
#' within-cluster, miR-17~92, miR-200, miR-371~373, and a dominant-vs-tail
#' abundance balance computed at runtime.
#'
#' Tolerant to feature absence: features missing from the panel are simply
#' removed from the numerator/denominator before computation; a balance
#' becomes NA only if either side becomes empty.
#'
#' The matrix/data.frame path is vectorized for deployment throughput. For
#' \code{aggregate = "gmean"}, it uses row-wise means of the logged feature
#' matrix and is regression-locked to the prior scalar row-loop semantics,
#' including large \code{"ALL_NON_RBC"} / \code{"ALL_NON_PLT"} denominators
#' and structurally generated SBPs. The same path also preserves feature
#' names for one-column matrices and avoids the prior multi-column
#' data.frame recursion failure.
#'
#' Backwards-compatibility note: the pre-existing
#' \code{"TOP_K_BY_ABUNDANCE"} behavior when the panel size is less than or
#' equal to \code{k} is intentionally preserved here for frozen within-sample
#' result identity. Treating that case as a semantic empty-tail denominator
#' would change existing scores and is reserved for a separate fix.
#'
#' @param x Numeric vector with miRNA names, OR a matrix (samples x features).
#'   Must have feature names - names are used to match the partition entries.
#' @param balances Named list of balances, each with components
#'   \code{numerator} and \code{denominator}. Defaults to
#'   \code{ws_default_sbp()}.
#' @param pseudocount Additive pseudocount before logging. Default 1e-6 of
#'   sample sum.
#' @param aggregate "gmean" (default) or "trimmed_gmean" (10\% trim each end).
#' @param min_balance_coverage Numeric in [0, 1]. Minimum fraction of the
#'   eight balances that must be computable from the input panel for the
#'   sample / cohort to be deemed eligible. When fewer than this fraction
#'   of balances have both numerator and denominator features present, the
#'   returned vector / matrix is filled with NA and carries
#'   the attribute \code{coverage_failed} set to TRUE. Default 0.8
#'   (within-sample compositional methods).
#'
#' @return Named numeric vector of balance values (single sample), or
#'   matrix samples x balances (matrix input). Carries attributes
#'   \code{coverage} (numeric in [0, 1]) and \code{coverage_failed}
#'   (logical) indicating whether the 80\% gate was satisfied.
#'
#' @examples
#' \dontrun{
#' x <- rlnorm(60, meanlog = 5, sdlog = 1.2)
#' names(x) <- ws_default_pivot_pool()[1:6]
#' ws_balance_ilr(x)  # NA where biology dictionary names don't match panel
#' }
#'
#' @export
ws_balance_ilr <- function(x,
                            balances = ws_default_sbp(),
                            pseudocount = NULL,
                            aggregate = c("gmean", "trimmed_gmean"),
                            min_balance_coverage = 0.8) {
  aggregate <- match.arg(aggregate)
  if (!is.numeric(min_balance_coverage) ||
      length(min_balance_coverage) != 1L ||
      is.na(min_balance_coverage) ||
      min_balance_coverage < 0 || min_balance_coverage > 1) {
    stop("ws_balance_ilr: min_balance_coverage must be a single value in [0, 1].")
  }

  if (is.matrix(x) || is.data.frame(x)) {
    x <- as.matrix(x)
    feat_names <- colnames(x)
    if (is.null(feat_names)) stop("x must have feature names (miRNA IDs).")

    balance_names <- names(balances)
    out <- matrix(NA_real_, nrow = nrow(x), ncol = length(balance_names),
                  dimnames = list(rownames(x), balance_names))
    if (is.null(pseudocount)) {
      pc <- 1e-6 * rowSums(x, na.rm = TRUE)
      pc[pc <= 0] <- .Machine$double.eps
    } else {
      pc <- pseudocount
      if (pc <= 0) pc <- .Machine$double.eps
    }
    log_x <- log(x + pc)

    for (j in seq_along(balance_names)) {
      b <- balances[[balance_names[[j]]]]
      num <- b$numerator
      den <- b$denominator

      if (length(num) == 1L && num == "TOP_K_BY_ABUNDANCE") {
        k <- if (!is.null(b$k)) b$k else 8L
        for (i in seq_len(nrow(x))) {
          ord <- order(x[i, ], decreasing = TRUE)
          num_idx <- ord[seq_len(min(k, length(x[i, ])))]
          den_idx <- ord[(min(k, length(x[i, ])) + 1L):length(ord)]
          if (length(num_idx) == 0L || length(den_idx) == 0L) next
          nN <- length(num_idx); nD <- length(den_idx)
          coef <- sqrt(nN * nD / (nN + nD))
          if (aggregate == "gmean") {
            num_mean <- mean(log_x[i, num_idx], na.rm = TRUE)
            den_mean <- mean(log_x[i, den_idx], na.rm = TRUE)
          } else {
            num_mean <- mean(log_x[i, num_idx], trim = 0.10, na.rm = TRUE)
            den_mean <- mean(log_x[i, den_idx], trim = 0.10, na.rm = TRUE)
          }
          out[i, j] <- coef * (num_mean - den_mean)
        }
        next
      }

      if (length(den) == 1L && den == "ALL_NON_RBC") {
        num_idx <- match(num, feat_names); num_idx <- num_idx[!is.na(num_idx)]
        den_idx <- setdiff(seq_along(feat_names), num_idx)
      } else if (length(den) == 1L && den == "ALL_NON_PLT") {
        num_idx <- match(num, feat_names); num_idx <- num_idx[!is.na(num_idx)]
        den_idx <- setdiff(seq_along(feat_names), num_idx)
      } else {
        num_idx <- match(num, feat_names); num_idx <- num_idx[!is.na(num_idx)]
        den_idx <- match(den, feat_names); den_idx <- den_idx[!is.na(den_idx)]
      }

      if (length(num_idx) == 0L || length(den_idx) == 0L) next

      nN <- length(num_idx); nD <- length(den_idx)
      coef <- sqrt(nN * nD / (nN + nD))
      if (aggregate == "gmean") {
        num_mean <- rowMeans(log_x[, num_idx, drop = FALSE], na.rm = TRUE)
        den_mean <- rowMeans(log_x[, den_idx, drop = FALSE], na.rm = TRUE)
      } else {
        num_mean <- apply(log_x[, num_idx, drop = FALSE], 1L, function(v) {
          mean(v, trim = 0.10, na.rm = TRUE)
        })
        den_mean <- apply(log_x[, den_idx, drop = FALSE], 1L, function(v) {
          mean(v, trim = 0.10, na.rm = TRUE)
        })
      }
      out[, j] <- coef * (num_mean - den_mean)
    }

    coverage <- if (ncol(out) == 0L) {
      rep(1, nrow(out))
    } else {
      unname(rowMeans(is.finite(out)))
    }
    coverage_failed <- coverage + sqrt(.Machine$double.eps) < min_balance_coverage
    if (any(coverage_failed)) out[coverage_failed, ] <- NA_real_
    rownames(out) <- rownames(x)
    attr(out, "coverage") <- coverage
    attr(out, "coverage_failed") <- coverage_failed
    attr(out, "min_balance_coverage") <- min_balance_coverage
    return(out)
  }

  feat_names <- names(x)
  if (is.null(feat_names)) stop("x must have feature names (miRNA IDs).")

  if (is.null(pseudocount)) pseudocount <- 1e-6 * sum(x, na.rm = TRUE)
  if (pseudocount <= 0) pseudocount <- .Machine$double.eps

  log_x <- log(x + pseudocount)

  agg_fn <- if (aggregate == "gmean") {
    function(v) mean(v, na.rm = TRUE)
  } else {
    function(v) mean(v, trim = 0.10, na.rm = TRUE)
  }

  out <- vapply(names(balances), function(bn) {
    b <- balances[[bn]]
    num <- b$numerator
    den <- b$denominator

    if (length(num) == 1L && num == "TOP_K_BY_ABUNDANCE") {
      k <- if (!is.null(b$k)) b$k else 8L
      ord <- order(x, decreasing = TRUE)
      num_idx <- ord[seq_len(min(k, length(x)))]
      den_idx <- ord[(min(k, length(x)) + 1L):length(ord)]
    } else if (length(den) == 1L && den == "ALL_NON_RBC") {
      num_idx <- match(num, feat_names); num_idx <- num_idx[!is.na(num_idx)]
      den_idx <- setdiff(seq_along(feat_names), num_idx)
    } else if (length(den) == 1L && den == "ALL_NON_PLT") {
      num_idx <- match(num, feat_names); num_idx <- num_idx[!is.na(num_idx)]
      den_idx <- setdiff(seq_along(feat_names), num_idx)
    } else {
      num_idx <- match(num, feat_names); num_idx <- num_idx[!is.na(num_idx)]
      den_idx <- match(den, feat_names); den_idx <- den_idx[!is.na(den_idx)]
    }

    if (length(num_idx) == 0L || length(den_idx) == 0L) return(NA_real_)

    nN <- length(num_idx); nD <- length(den_idx)
    coef <- sqrt(nN * nD / (nN + nD))
    coef * (agg_fn(log_x[num_idx]) - agg_fn(log_x[den_idx]))
  }, numeric(1L))

  coverage <- if (length(out) == 0L) 1 else mean(is.finite(out))
  if (coverage + sqrt(.Machine$double.eps) < min_balance_coverage) {
    out[] <- NA_real_
    attr(out, "coverage_failed") <- TRUE
  } else {
    attr(out, "coverage_failed") <- FALSE
  }
  attr(out, "coverage") <- coverage
  attr(out, "min_balance_coverage") <- min_balance_coverage
  out
}


# ----------------------------------------------------------------------------
# ws_alr_pivot
# ----------------------------------------------------------------------------

#' @title Default ALR pivot pool for circulating-miRNA panels (v1)
#'
#' @description
#' Returns the curated v1 pivot pool used as the geometric-mean denominator
#' by \code{\link{ws_alr_pivot}}. Six low-variance circulating miRNAs with
#' no documented contamination from haemolysis (miR-451a/16/486/144/223) or
#' platelet activation (miR-223 family). The pool is fixed at v1; any
#' revision will be tagged \code{circulating_v2}, etc., with a separate
#' justification log under
#' the single-sample scoring bank's pivot-pool justification notes.
#'
#' Notably absent from the v1 pool: miR-92a-3p (recent serum/plasma
#' literature treats it as erythrocyte-derived; codex Round 1 review of
#' the single-sample scoring bank plan v0.2 flagged its prior inclusion as a concern).
#'
#' @return Character vector of six miRNA IDs.
#' @export
ws_default_pivot_pool <- function() {
  c(
    "hsa-miR-103a-3p",
    "hsa-miR-191-5p",
    "hsa-miR-26a-5p",
    "hsa-miR-30c-5p",
    "hsa-let-7g-5p",
    "hsa-miR-93-5p"
  )
}


#' @title Additive log-ratio with a frozen pivot pool
#'
#' @description
#' Within-sample additive log-ratio (ALR) transformation using a frozen
#' pivot pool of housekeeping-like miRNAs as the geometric-mean denominator.
#' \deqn{x_{\mathrm{alr},i} = \log(x_i + c) - \log(\mathrm{gmean}(x_{\mathrm{pivot}} + c))}
#' Anchors the transformation against a stable internal reference rather
#' than choosing a single feature as denominator (the classical Aitchison
#' ALR).
#'
#' Fail-closed by default: refuses to compute when fewer than
#' \code{min_pivot_present} pivots are present in the input panel, unless
#' the operator explicitly opts in to a global rCLR fallback. Silent
#' substitution of all-feature centering for a curated denominator can
#' mask haemolysis or platelet contamination - flagged by codex Round 1
#' review.
#'
#' @details
#' \strong{Backwards-compatibility / single-sample-pipeline note (v2.4.0).}
#' The package default is \code{allow_global_fallback = FALSE} (fail-closed):
#' the function errors out rather than silently substituting a global rCLR
#' centering when fewer than \code{min_pivot_present} pivots are available.
#' This is the conservative choice for end users running interactive
#' analyses, because it prevents haemolysis-contaminated cells from being
#' silently rescued by an all-feature centering that includes the very
#' contaminants they were meant to flag.
#'
#' The single-sample pipeline (within-sample compositional methods) describes
#' the fallback as having been "activated" with the
#' affected cells explicitly flagged in the per-cell results. The pipeline
#' opts in via \code{allow_global_fallback = TRUE} on every call. This
#' alignment is intentional: the package keeps the safer default for ad-hoc
#' use, while the single-sample pipeline has the explicit logging discipline
#' to track every cell that took the fallback path.
#'
#' @param x Samples x features matrix (rows = samples). Must have column
#'   names to identify pivot features.
#' @param pivot_features Character vector of feature names to use as pivot
#'   pool. If \code{NULL}, defaults to \code{ws_default_pivot_pool()}.
#' @param pseudocount Additive pseudocount. Default 0.5.
#' @param allow_global_fallback Logical. If \code{TRUE}, falls back to a
#'   global rCLR (all features) when fewer than \code{min_pivot_present}
#'   pivots are present. Default \code{FALSE} (fail-closed).
#' @param min_pivot_present Minimum pivot features required. Default 3.
#' @param verbose Logical. If \code{TRUE} and fallback is taken, emits a
#'   message.
#'
#' @return Same shape as \code{x}; the ALR-transformed values, with
#'   attributes \code{pivot_features}, \code{pseudocount}, and
#'   \code{fallback_used}.
#'
#' @export
ws_alr_pivot <- function(x,
                         pivot_features = NULL,
                         pseudocount = 0.5,
                         allow_global_fallback = FALSE,
                         min_pivot_present = 3L,
                         verbose = FALSE) {
  if (!is.matrix(x)) x <- as.matrix(x)
  if (is.null(colnames(x))) {
    stop("ws_alr_pivot requires named features (column names) to identify pivot pool")
  }
  if (any(x < 0, na.rm = TRUE)) {
    stop("ws_alr_pivot: negative values not allowed (compositional input)")
  }
  if (pseudocount <= 0) stop("ws_alr_pivot: pseudocount must be > 0")
  if (is.null(pivot_features)) pivot_features <- ws_default_pivot_pool()

  pivot_present <- intersect(pivot_features, colnames(x))
  fallback_used <- FALSE
  if (length(pivot_present) < min_pivot_present) {
    if (!allow_global_fallback) {
      stop("ws_alr_pivot: only ", length(pivot_present), " of ", length(pivot_features),
            " pivot features present (need >=", min_pivot_present,
            "). Pass allow_global_fallback=TRUE to opt into global rCLR fallback.")
    }
    if (verbose) message("ws_alr_pivot: <", min_pivot_present,
                          " pivot features present - falling back to global rCLR")
    pivot_present <- colnames(x)
    fallback_used <- TRUE
  }

  x_pos <- x + pseudocount
  pivot_mat <- x_pos[, pivot_present, drop = FALSE]
  log_pivot_geomean <- rowMeans(log(pivot_mat), na.rm = TRUE)

  out <- log(x_pos) - log_pivot_geomean
  attr(out, "pivot_features") <- pivot_present
  attr(out, "pseudocount") <- pseudocount
  attr(out, "fallback_used") <- fallback_used
  out
}


# ----------------------------------------------------------------------------
# ws_mad_logratio
# ----------------------------------------------------------------------------

#' @title Median-centred log-ratio with optional MAD scaling
#'
#' @description
#' Robust descriptive transform on log-abundances: subtracts the per-sample
#' median log-abundance (more robust to outlier features than mean-centering
#' as in CLR/rCLR) and optionally divides by the per-sample median absolute
#' deviation (MAD; limits leverage from individual high-magnitude features).
#'
#' This is NOT an Aitchison-valid simplex projection. It is a robust
#' descriptive transform on log-abundances, suitable as input to learners
#' that do not require the zero-sum compositional coordinate constraint.
#'
#' @param x Samples x features matrix (rows = samples).
#' @param pseudocount Additive pseudocount before logging. Default 0.5.
#' @param scale_by_mad Logical. If \code{TRUE}, divides each sample by its
#'   per-sample MAD after median-centering. Default \code{TRUE}.
#'
#' @return Same shape as \code{x}; the median-centred (MAD-scaled)
#'   log-abundance values, with attributes \code{centering} and
#'   \code{scaling}.
#'
#' @export
ws_mad_logratio <- function(x,
                             pseudocount = 0.5,
                             scale_by_mad = TRUE) {
  if (!is.matrix(x)) x <- as.matrix(x)
  if (any(x < 0, na.rm = TRUE)) {
    stop("ws_mad_logratio: negative values not allowed (compositional input)")
  }
  if (pseudocount <= 0) stop("ws_mad_logratio: pseudocount must be > 0")
  log_x <- log(x + pseudocount)

  per_sample_median <- apply(log_x, 1L, median, na.rm = TRUE)
  if (any(is.na(per_sample_median))) {
    stop("ws_mad_logratio: ", sum(is.na(per_sample_median)),
          " sample(s) have all-NA log values - cannot compute median")
  }
  centered <- sweep(log_x, 1L, per_sample_median, FUN = "-")

  if (scale_by_mad) {
    per_sample_mad <- apply(log_x, 1L, mad, na.rm = TRUE)
    per_sample_mad[is.na(per_sample_mad) | per_sample_mad < 1e-6] <- 1
    out <- sweep(centered, 1L, per_sample_mad, FUN = "/")
  } else {
    out <- centered
  }
  attr(out, "centering") <- "median_log"
  attr(out, "scaling") <- if (scale_by_mad) "mad" else "none"
  out
}


# ----------------------------------------------------------------------------
# ws_dominance_score, fit_dominance_threshold, ws_dominance_flag
# ----------------------------------------------------------------------------

#' @title Panel-internal dominance score (within-sample QC)
#'
#' @description
#' Reports the fraction of total panel signal carried by the top-K features
#' in a single sample, on the raw abundance scale. Dimensionless,
#' scale-invariant, in [K/N, 1]. Sera with extreme dominance (e.g. top-2
#' features carrying > 60% of signal) often have pre-analytical issues
#' (massive haemolysis, incomplete extraction).
#'
#' This is a SCALAR QC metric, not a normalization peer. It is designed for
#' use as an input to claim-gating in Module D, not as a Module-A
#' normalization method comparable to rCLR / ILR / ALR / MAD-logratio.
#'
#' @param x Samples x features matrix.
#' @param top_k Integer; number of top features to include in the
#'   dominance numerator. Default 2.
#' @param pseudocount Additive pseudocount. Default 0.5.
#'
#' @return Numeric vector of length \code{nrow(x)} with values in [0, 1].
#'
#' @export
ws_dominance_score <- function(x,
                                top_k = 2L,
                                pseudocount = 0.5) {
  if (!is.matrix(x)) x <- as.matrix(x)
  if (any(x < 0, na.rm = TRUE)) {
    stop("ws_dominance_score: negative values not allowed (compositional input)")
  }
  if (pseudocount < 0) stop("ws_dominance_score: pseudocount must be >= 0")
  if (top_k < 1L) stop("ws_dominance_score: top_k must be >= 1")
  x_pos <- x + pseudocount

  apply(x_pos, 1L, function(row) {
    s <- sort(row, decreasing = TRUE, na.last = NA)
    if (length(s) == 0L) return(NA_real_)
    sum(s[seq_len(min(top_k, length(s)))]) / sum(s)
  })
}


#' @title Calibrate a cohort-specific dominance threshold
#'
#' @description
#' Calibrates a cohort-specific threshold for \code{ws_dominance_flag} as
#' the (1 - alpha) quantile of \code{ws_dominance_score} computed on a
#' training reference cohort. Recommended over the universal 0.60 default,
#' which is provided only for backwards compatibility.
#'
#' @param x_train Samples x features matrix of the training reference
#'   cohort (typically Tier R, control samples only).
#' @param top_k Integer; passed to \code{ws_dominance_score}. Default 2.
#' @param alpha Numeric in (0, 1). Threshold is the (1 - alpha) quantile of
#'   the training dominance distribution, i.e. samples in the top
#'   \code{alpha} fraction of the training distribution are flagged.
#'   Default 0.05.
#' @param pseudocount Numeric; passed to \code{ws_dominance_score}.
#'   Default 0.5.
#'
#' @return List with components \code{threshold}, \code{alpha},
#'   \code{top_k}, \code{reference_n}, \code{reference_quantiles}.
#' @export
fit_dominance_threshold <- function(x_train, top_k = 2L,
                                     alpha = 0.05, pseudocount = 0.5) {
  scores <- ws_dominance_score(x_train, top_k = top_k, pseudocount = pseudocount)
  thr <- as.numeric(stats::quantile(scores, probs = 1 - alpha, na.rm = TRUE))
  list(threshold = thr, alpha = alpha, top_k = top_k,
       reference_n = sum(!is.na(scores)),
       reference_quantiles = stats::quantile(scores, c(0.5, 0.9, 0.95, 0.99), na.rm = TRUE))
}


#' @title Flag samples whose dominance score exceeds a calibrated threshold
#'
#' @description
#' Convenience wrapper: returns a logical vector of length \code{nrow(x)}
#' indicating which samples have a dominance score exceeding
#' \code{threshold}. Use \code{\link{fit_dominance_threshold}} to calibrate
#' the threshold from a training reference cohort.
#'
#' @inheritParams ws_dominance_score
#' @param threshold Numeric. Default 0.60.
#'
#' @return Logical vector of length \code{nrow(x)}.
#' @export
ws_dominance_flag <- function(x,
                                top_k = 2L,
                                threshold = 0.60,
                                pseudocount = 0.5) {
  ws_dominance_score(x, top_k = top_k, pseudocount = pseudocount) > threshold
}
