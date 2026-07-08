#' Paper 3 method-bank registry
#'
#' Returns the package-facing implementation map for the OmicSelector manuscript
#' method bank. The registry is intentionally executable metadata: each row
#' names the exported fit, predict, score, or helper function that implements
#' the manuscript method or auxiliary protocol.
#'
#' @return A data frame with one row per method or auxiliary protocol.
#' @export
singlesample_method_bank <- function() {
  rows <- list(
    .p3_row("A", "rCLR trimmed", NA, "ws_rclr_trimmed", "evaluable"),
    .p3_row("A", "ILR balance dictionary", NA, "ws_balance_ilr", "available_user_panels"),
    .p3_row("A", "ALR pivot pool", NA, "ws_alr_pivot", "evaluable"),
    .p3_row("A", "MAD log-ratio", NA, "ws_mad_logratio", "evaluable"),
    .p3_row("A", "dominance score", "fit_dominance_threshold", "ws_dominance_score", "evaluable"),

    .p3_row("B", "frozen RUV", "fit_frozen_ruv", "apply_frozen_ruv", "evaluable"),
    .p3_row("B", "robust PCA residual", "fit_robust_pca_residual", "apply_robust_pca_residual", "evaluable"),
    .p3_row("B", "hemolysis residual removal", "fit_hemolysis_rr", "apply_hemolysis_rr", "evaluable"),

    .p3_row("C", "block fixed-effect rCLR", "fit_block_fe_rclr", "score_block_fe_rclr", "non_evaluable_loco_corpus"),
    .p3_row("C", "cohort-z then pool", "fit_cohort_z_then_pool_score", "score_cohort_z_then_pool_score", "non_evaluable_loco_corpus"),
    .p3_row("C", "cluster-robust ensemble", NA, "score_cluster_robust_ensemble", "non_evaluable_loco_corpus"),
    .p3_row("C", "mixed-effects scorer", "fit_mixed_effects_scorer", "score_mixed_effects_scorer", "non_evaluable_loco_corpus"),

    .p3_row("D", "kit-stratified rCLR", "fit_kit_stratified_rclr", "score_kit_stratified_rclr", "evaluable"),
    .p3_row("D", "kit fixed-effect ALR", "fit_kit_fe_adjusted_alr", "score_kit_fe_adjusted_alr", "evaluable"),
    .p3_row("D", "kit-orthogonal ILR", "fit_kit_orthogonal_ilr", "score_kit_orthogonal_ilr", "evaluable"),
    .p3_row("D", "kit-residual MAD", "fit_kit_residual_mad", "score_kit_residual_mad", "evaluable"),
    .p3_row("D", "ComBat-seq transductive comparator", NA, "score_combat_seq_then_rclr", "exploratory_comparator"),

    .p3_row("E", "kit-stable anchor rCLR", "fit_kit_stable_anchors", "score_kit_stable_anchor_rclr", "evaluable"),
    .p3_row("E", "hemolysis-kit dual anchor", "fit_hemolysis_kit_dual_anchor", "score_hemolysis_kit_dual_anchor", "evaluable"),
    .p3_row("E", "RIN-weighted reference", "fit_rin_weighted_reference", "score_rin_weighted_score", "evaluable"),

    .p3_row("F", "technology-stratified rCLR", "fit_tech_stratified_rclr", "score_tech_stratified_rclr", "evaluable"),
    .p3_row("F", "cross-technology harmonized rCLR", "fit_cross_tech_harmonized_rclr", "score_cross_tech_harmonized_rclr", "evaluable"),
    .p3_row("F", "technology-residualized ALR", "fit_tech_residualized_alr", "score_tech_residualized_alr", "evaluable"),

    .p3_row("G", "biofluid-stratified rCLR", "fit_biofluid_stratified_rclr", "score_biofluid_stratified_rclr", "evaluable"),
    .p3_row("G", "biofluid-stable anchor rCLR", "fit_biofluid_anchor_rclr", "score_biofluid_anchor_rclr", "evaluable"),
    .p3_row("G", "biofluid-residualized ALR", "fit_biofluid_residualized_alr", "score_biofluid_residualized_alr", "evaluable"),
    .p3_row("G", "blood-cell-marker exclusion rCLR", NA, "score_platelet_exclusion_rclr", "evaluable"),

    .p3_row("H", "DANN with kit covariate", "os_fit_dann_kit_extended", "os_score_dann_kit_extended", "evaluable"),
    .p3_row("H", "kit-conditional VAE", "os_fit_kit_conditional_vae", "os_score_kit_conditional_vae", "evaluable"),
    .p3_row("H", "per-kit ICP", "os_fit_icp_per_kit", "os_score_icp_per_kit", "evaluable"),

    .p3_row("I", "Group-DRO kit by biofluid", "os_fit_group_dro_scorer", "os_score_group_dro_scorer", "attempted"),
    .p3_row("I", "Sinkhorn OT cohort barycenter", "os_fit_sinkhorn_ot_scorer", "os_score_sinkhorn_ot_scorer", "attempted"),
    .p3_row("I", "TabPFN-v2 learner", "make_tabpfn_learner", NA, "deferred_interactive_license_token"),

    .p3_row("J", "qPCR non-detect imputation", "qpcr_nondetect_impute", "qpcr_nondetect_lod_fallback", "auxiliary"),
    .p3_row("J", "subtype-stratified matched random panels", NA, "singlesample_matched_null_benchmark_cv", "auxiliary"),
    .p3_row("J", "LOCO splits with same-block exclusion", NA, "singlesample_make_loco_splits", "auxiliary"),
    .p3_row("J", "LOCTO splits with same-block exclusion", NA, "singlesample_make_locto_splits", "auxiliary"),
    .p3_row("J", "provenance pre-flight gate", NA, "os_provenance_preflight", "auxiliary"),
    .p3_row("J", "block-aware Brown/Kost BH", NA, "singlesample_bh_fdr_correct_blocked", "auxiliary"),
    .p3_row("J", "MIMAT-to-hsa-miR namespace lookup", "resolve_mimat_to_hsa", "resolve_hsa_to_mimat", "auxiliary")
  )
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  exports <- tryCatch(getNamespaceExports("OmicSelector"), error = function(e) character())
  funcs <- unique(stats::na.omit(c(out$fit_function, out$score_function)))
  existing <- vapply(funcs, exists, logical(1), envir = parent.env(environment()), mode = "function")
  exported <- funcs %in% exports
  out$functions_exist <- vapply(seq_len(nrow(out)), function(i) {
    f <- stats::na.omit(c(out$fit_function[[i]], out$score_function[[i]]))
    all(f %in% funcs[existing])
  }, logical(1))
  out$functions_exported <- vapply(seq_len(nrow(out)), function(i) {
    f <- stats::na.omit(c(out$fit_function[[i]], out$score_function[[i]]))
    all(f %in% funcs[exported])
  }, logical(1))
  out <- .singlesample_annotate_bank_roster(out)
  out
}

# Annotate the export-check bank with the Amendment #4 roster routing metadata.
# Rows are matched to the roster by score function; auxiliary / non-roster rows
# (e.g. family J protocols, non-evaluable family C, the batchwise out-N Sinkhorn
# scorer) receive NA. Backward compatible: existing columns are unchanged and no
# rows are added or removed.
.singlesample_annotate_bank_roster <- function(out) {
  ros <- tryCatch(singlesample_method_roster(), error = function(e) NULL)
  add_cols <- c("method_id", "estimand", "role", "single_sample_gate",
                "negative_control", "row_source", "lopbo_mechanism")
  if (is.null(ros)) {
    for (col in add_cols) {
      out[[col]] <- if (col == "negative_control") NA else NA_character_
    }
    return(out)
  }
  idx <- match(out$score_function, ros$score_fn)
  out$method_id         <- ros$method_id[idx]
  out$estimand          <- ros$estimand[idx]
  out$role              <- ros$role[idx]
  out$single_sample_gate <- ifelse(is.na(idx), NA_character_, "required")
  out$negative_control  <- ros$is_negative_control[idx]
  out$row_source        <- ros$row_source[idx]
  out$lopbo_mechanism   <- ros$lopbo_mechanism[idx]
  out
}

.p3_row <- function(family, method, fit_function, score_function, status) {
  data.frame(
    family = family,
    method = method,
    fit_function = fit_function %||% NA_character_,
    score_function = score_function %||% NA_character_,
    status = status,
    stringsAsFactors = FALSE
  )
}

#' Assert that Paper 3 method-bank functions are present and exported
#'
#' @param bank Optional registry returned by [singlesample_method_bank()].
#' @param allow_deferred Logical. If `TRUE`, deferred methods such as TabPFN-v2
#'   are not treated as failures.
#' @return Invisibly returns the registry if all required functions are exported.
#' @export
singlesample_assert_method_bank_exports <- function(bank = singlesample_method_bank(),
                                              allow_deferred = TRUE) {
  required <- bank
  if (isTRUE(allow_deferred)) {
    required <- required[!grepl("^deferred", required$status), , drop = FALSE]
  }
  bad <- required[!required$functions_exist | !required$functions_exported, , drop = FALSE]
  if (nrow(bad) > 0L) {
    stop("Paper 3 method-bank export check failed for: ",
         paste(paste0(bad$family, "/", bad$method), collapse = "; "),
         call. = FALSE)
  }
  invisible(bank)
}

#' Build leave-one-cohort-out splits with same-block exclusion
#'
#' @param sample_meta Data frame with one row per sample.
#' @param cohort_col Column identifying cohorts.
#' @param provenance_block_col Optional column identifying specimen-provenance
#'   blocks. When present, samples in the held-out block are excluded from
#'   training.
#' @param min_train_cohorts Minimum number of training cohorts required for an
#'   evaluable split.
#' @return A list of split objects with train/test indices and status.
#' @export
singlesample_make_loco_splits <- function(sample_meta,
                                    cohort_col = "cohort",
                                    provenance_block_col = "provenance_block",
                                    min_train_cohorts = 1L) {
  .p3_validate_split_meta(sample_meta, cohort_col)
  meta <- as.data.frame(sample_meta, stringsAsFactors = FALSE)
  cohorts <- unique(as.character(meta[[cohort_col]]))
  has_block <- provenance_block_col %in% names(meta)
  lapply(cohorts, function(held) {
    test <- which(as.character(meta[[cohort_col]]) == held)
    train <- setdiff(seq_len(nrow(meta)), test)
    held_block <- NA_character_
    if (has_block) {
      held_blocks <- unique(as.character(meta[[provenance_block_col]][test]))
      held_blocks <- held_blocks[!is.na(held_blocks) & nzchar(held_blocks)]
      held_block <- paste(held_blocks, collapse = ",")
      train <- train[!(as.character(meta[[provenance_block_col]][train]) %in% held_blocks)]
    }
    n_train_cohorts <- length(unique(as.character(meta[[cohort_col]][train])))
    status <- if (n_train_cohorts >= min_train_cohorts) "evaluable" else "not_evaluable_insufficient_train_cohorts"
    list(held_out = held, held_out_block = held_block, train = train,
         test = test, n_train_cohorts = n_train_cohorts, status = status)
  })
}

#' Build leave-one-cancer-type-out splits with same-block exclusion
#'
#' @param sample_meta Data frame with one row per sample.
#' @param cancer_col Column identifying cancer type or endpoint family.
#' @inheritParams singlesample_make_loco_splits
#' @return A list of split objects with train/test indices and status.
#' @export
singlesample_make_locto_splits <- function(sample_meta,
                                     cancer_col = "cancer_type",
                                     cohort_col = "cohort",
                                     provenance_block_col = "provenance_block",
                                     min_train_cohorts = 1L) {
  .p3_validate_split_meta(sample_meta, cancer_col)
  .p3_validate_split_meta(sample_meta, cohort_col)
  meta <- as.data.frame(sample_meta, stringsAsFactors = FALSE)
  cancer_types <- unique(as.character(meta[[cancer_col]]))
  has_block <- provenance_block_col %in% names(meta)
  lapply(cancer_types, function(held) {
    test <- which(as.character(meta[[cancer_col]]) == held)
    train <- setdiff(seq_len(nrow(meta)), test)
    held_blocks <- character()
    if (has_block) {
      held_blocks <- unique(as.character(meta[[provenance_block_col]][test]))
      held_blocks <- held_blocks[!is.na(held_blocks) & nzchar(held_blocks)]
      train <- train[!(as.character(meta[[provenance_block_col]][train]) %in% held_blocks)]
    }
    n_train_cohorts <- length(unique(as.character(meta[[cohort_col]][train])))
    status <- if (n_train_cohorts >= min_train_cohorts) "evaluable" else "not_evaluable_insufficient_train_cohorts"
    list(held_out_cancer_type = held, held_out_blocks = held_blocks,
         train = train, test = test, n_train_cohorts = n_train_cohorts,
         status = status)
  })
}

.p3_validate_split_meta <- function(sample_meta, col) {
  if (!is.data.frame(sample_meta)) {
    stop("sample_meta must be a data frame.", call. = FALSE)
  }
  if (!(col %in% names(sample_meta))) {
    stop("sample_meta is missing column: ", col, call. = FALSE)
  }
  invisible(TRUE)
}

#' Advanced Paper 3 learned and transport scorers
#'
#' Thin exported aliases for the learned kit-aware, Group-DRO, and Sinkhorn-OT
#' scorer implementations used by the OmicSelector manuscript. The aliases keep
#' experimental method names under the `os_` namespace while preserving the
#' train/apply separation in the underlying fit and score functions.
#'
#' @param ... Arguments passed to the underlying method-specific fit, transform,
#'   apply, or score function.
#' @name singlesample_advanced_scorers
NULL

#' @rdname singlesample_advanced_scorers
#' @export
os_fit_dann_kit_extended <- function(...) fit_dann_kit_extended(...)

#' @rdname singlesample_advanced_scorers
#' @export
os_score_dann_kit_extended <- function(...) score_dann_kit_extended(...)

#' @rdname singlesample_advanced_scorers
#' @export
os_transform_dann_kit_extended <- function(...) transform_dann_kit_extended(...)

#' @rdname singlesample_advanced_scorers
#' @export
os_fit_kit_conditional_vae <- function(...) fit_kit_conditional_vae(...)

#' @rdname singlesample_advanced_scorers
#' @export
os_score_kit_conditional_vae <- function(...) score_kit_conditional_vae(...)

#' @rdname singlesample_advanced_scorers
#' @export
os_transform_kit_conditional_vae <- function(...) transform_kit_conditional_vae(...)

#' @rdname singlesample_advanced_scorers
#' @export
os_fit_icp_per_kit <- function(...) fit_icp_per_kit(...)

#' @rdname singlesample_advanced_scorers
#' @export
os_score_icp_per_kit <- function(...) score_icp_per_kit(...)

#' @rdname singlesample_advanced_scorers
#' @export
os_predict_icp_per_kit_details <- function(...) predict_icp_per_kit_details(...)

#' @rdname singlesample_advanced_scorers
#' @export
os_fit_group_dro_scorer <- function(...) fit_group_dro_scorer(...)

#' @rdname singlesample_advanced_scorers
#' @export
os_score_group_dro_scorer <- function(...) score_group_dro_scorer(...)

#' @rdname singlesample_advanced_scorers
#' @export
os_fit_sinkhorn_ot_scorer <- function(...) fit_sinkhorn_ot_scorer(...)

#' @rdname singlesample_advanced_scorers
#' @export
os_apply_sinkhorn_ot_scorer <- function(...) apply_sinkhorn_ot_scorer(...)

#' @rdname singlesample_advanced_scorers
#' @export
os_score_sinkhorn_ot_scorer <- function(...) score_sinkhorn_ot_scorer(...)
