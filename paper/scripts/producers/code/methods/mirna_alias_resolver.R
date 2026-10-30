# mirna_alias_resolver.R
#
# Pipeline-side wrapper that loads the package's MIMAT alias resolver
# (paper3-mirna-name-resolver.R, copied into mirna_alias_resolver_core.R)
# and adds a scale-aware collapse helper that the run_paper3_pipeline.R
#
# Codex review (2026-05-05) rejected naive sum-collapse on log-scale data:
# - Toray 3D-Gene deposits and most public microarray matrices are
#   already log_2 normalised. Summing log-scale probe values is wrong.
# - Multiple Toray probes targeting the same canonical miRNA encode
#   replicate measurement, not disjoint molecular components, so a
#   median (on log-scale) or sum (on abundance scale) is correct
#   depending on the scale state.
#
# This file is sourced by code/run_paper3_pipeline.R after the
# preprocess_inverse_log routing decision but before Module A methods.

source(file.path("code/methods/mirna_alias_resolver_core.R"))


# ---------------------------------------------------------------------------
# resolve_and_collapse_features
# ---------------------------------------------------------------------------

#' Rename feature columns of a samples × features matrix to canonical
#' miRBase v22.1 mature names where possible, then collapse duplicate
#' canonical columns with a scale-aware aggregator.
#'
#' @param expr_per_sample Matrix with samples in rows and features in
#'   columns (`colnames` are platform feature IDs — MIMAT, hsa-miR,
#'   probe IDs, etc.).
#' @param scale_state One of "log" (default; median collapse on log scale),
#'   "abundance" (sum collapse; opt-in for absolute-count data), or
#'   "auto" (infer from sign of values; reference fallback — not recommended
#'   for Tier A / Track A Tier R where scale is known from declared metadata).
#' @return Matrix with the same row count, fewer or equal column count,
#'   `colnames` set to canonical names where resolved (with original
#'   names retained for unresolved features). Carries attributes:
#'   `resolved_mask`, `scale_state`, `n_resolved`, `n_unresolved`,
#'   `n_collapsed_groups`, `collapse_method`.
resolve_and_collapse_features <- function(expr_per_sample,
                                          scale_state = c("log", "abundance",
                                                           "auto")) {
  scale_state <- match.arg(scale_state)
  if (!is.matrix(expr_per_sample)) expr_per_sample <- as.matrix(expr_per_sample)
  if (is.null(colnames(expr_per_sample))) {
    stop("resolve_and_collapse_features: matrix must have column names")
  }

  if (scale_state == "auto") {
    has_neg <- any(expr_per_sample < 0, na.rm = TRUE)
    scale_state <- if (isTRUE(has_neg)) "log" else "abundance"
  }
  collapse_method <- if (scale_state == "log") "median_log" else "sum_abundance"

  feats <- colnames(expr_per_sample)
  resolved <- resolve_mirna_aliases(feats,
                                    target_namespace = "mirna_name",
                                    keep_unresolved  = TRUE,
                                    verbose          = FALSE)
  resolved_mask <- !is.na(resolved)
  final_names   <- ifelse(resolved_mask, resolved, feats)

  groups <- split(seq_along(final_names), final_names)
  collapsed_groups <- 0L

  if (max(lengths(groups)) == 1L) {
    out <- expr_per_sample
    colnames(out) <- final_names
    attr(out, "resolved_mask")     <- resolved_mask
    attr(out, "scale_state")       <- scale_state
    attr(out, "collapse_method")   <- collapse_method
    attr(out, "n_resolved")        <- sum(resolved_mask)
    attr(out, "n_unresolved")      <- sum(!resolved_mask)
    attr(out, "n_collapsed_groups")<- 0L
    attr(out, "n_input_features") <- length(feats)
    attr(out, "n_output_features")<- ncol(out)
    return(out)
  }

  out_cols <- vector("list", length(groups))
  out_resolved_mask <- logical(length(groups))
  for (g in seq_along(groups)) {
    idx <- groups[[g]]
    if (length(idx) == 1L) {
      out_cols[[g]] <- expr_per_sample[, idx, drop = TRUE]
    } else {
      sub <- expr_per_sample[, idx, drop = FALSE]
      if (scale_state == "log") {
        out_cols[[g]] <- apply(sub, 1L, function(v) median(v, na.rm = TRUE))
      } else {
        out_cols[[g]] <- rowSums(sub, na.rm = TRUE)
      }
      collapsed_groups <- collapsed_groups + 1L
    }
    out_resolved_mask[g] <- all(resolved_mask[idx])
  }
  out <- do.call(cbind, out_cols)
  colnames(out) <- names(groups)
  rownames(out) <- rownames(expr_per_sample)
  storage.mode(out) <- "double"

  attr(out, "resolved_mask")      <- out_resolved_mask
  attr(out, "scale_state")        <- scale_state
  attr(out, "collapse_method")    <- collapse_method
  attr(out, "n_resolved")         <- sum(resolved_mask)
  attr(out, "n_unresolved")       <- sum(!resolved_mask)
  attr(out, "n_collapsed_groups") <- collapsed_groups
  attr(out, "n_input_features")  <- length(feats)
  attr(out, "n_output_features") <- ncol(out)
  out
}


# ---------------------------------------------------------------------------
# ws_balance_ilr_coverage
# ---------------------------------------------------------------------------

#' Compute SBP-balance coverage statistics for a (canonical) feature set.
#'
#' The biology-frozen sequential binary partition uses canonical mature
#' miRBase v22.1 names. For each balance, count whether the numerator
#' and denominator sides have at least one feature in the canonical
#' subset of the cohort matrix.
#'
#' @param canonical_features character vector of canonical names (subset
#'   of cohort feature names that resolved to the mirna_name namespace).
#' @param sbp the SBP list as returned by `ws_default_sbp()`.
#' @return list with `n_balances`, `n_balances_with_coverage`,
#'   `coverage_fraction`, and per-balance details.
ws_balance_ilr_coverage <- function(canonical_features, sbp) {
  results <- lapply(seq_along(sbp), function(i) {
    b <- sbp[[i]]
    num <- b$numerator
    den <- b$denominator
    num_match <- sum(num %in% canonical_features)
    den_match <- if (length(den) == 1L && den %in% c("ALL_NON_RBC", "ALL_NON_PLT",
                                                       "TOP_K_BY_ABUNDANCE",
                                                       "TAIL_BY_ABUNDANCE"))
      length(canonical_features) - num_match else sum(den %in% canonical_features)
    list(name = names(sbp)[i] %||% paste0("balance_", i),
         num_n = num_match,
         den_n = den_match,
         valid = (num_match >= 1L) && (den_match >= 1L))
  })
  n_valid <- sum(vapply(results, function(r) isTRUE(r$valid), logical(1L)))
  list(n_balances = length(sbp),
       n_balances_with_coverage = n_valid,
       coverage_fraction = n_valid / max(1L, length(sbp)),
       per_balance = results)
}

# Local null-coalesce so this file is sourcing-safe before run_paper3_pipeline.R
# defines its own.
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  set.seed(42L)
  feats <- c("MIMAT0001631", "MIMAT0000062",   # canonical hsa-miR-451a, hsa-let-7a-5p
             "hsa-miR-451a", "hsa-let-7a-5p",   # already canonical → duplicate after resolve
             "junk_probe_001", "junk_probe_002")
  log_x <- matrix(rnorm(20L * length(feats), mean = 0, sd = 1.5),
                   nrow = 20L, ncol = length(feats),
                   dimnames = list(paste0("S", seq_len(20L)), feats))
  out_log <- resolve_and_collapse_features(log_x, scale_state = "log")
  cat(sprintf("[self-test] log-scale: in=%d → out=%d, collapsed_groups=%d, scale=%s, method=%s\n",
              ncol(log_x), ncol(out_log), attr(out_log, "n_collapsed_groups"),
              attr(out_log, "scale_state"), attr(out_log, "collapse_method")))
  stopifnot(attr(out_log, "n_collapsed_groups") == 2L)
  stopifnot(attr(out_log, "collapse_method") == "median_log")

  abs_x <- 2 ^ log_x
  out_abs <- resolve_and_collapse_features(abs_x, scale_state = "abundance")
  cat(sprintf("[self-test] abundance: in=%d → out=%d, collapsed_groups=%d, scale=%s, method=%s\n",
              ncol(abs_x), ncol(out_abs), attr(out_abs, "n_collapsed_groups"),
              attr(out_abs, "scale_state"), attr(out_abs, "collapse_method")))
  stopifnot(attr(out_abs, "n_collapsed_groups") == 2L)
  stopifnot(attr(out_abs, "collapse_method") == "sum_abundance")

  # Default (log) collapses by median without needing explicit scale_state
  out_default <- resolve_and_collapse_features(log_x)
  cat(sprintf("[self-test] default (log): scale=%s, method=%s\n",
              attr(out_default, "scale_state"), attr(out_default, "collapse_method")))
  stopifnot(attr(out_default, "scale_state") == "log")
  stopifnot(attr(out_default, "collapse_method") == "median_log")

  out_auto <- resolve_and_collapse_features(log_x, scale_state = "auto")
  cat(sprintf("[self-test] auto on log: detected scale=%s\n",
              attr(out_auto, "scale_state")))
  stopifnot(attr(out_auto, "scale_state") == "log")

  cat(sprintf("[self-test] resolver wrapper PASS\n"))
}
