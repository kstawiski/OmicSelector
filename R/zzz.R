# Define the null-coalescing operator %||%
# Returns lhs if not NULL, otherwise rhs
# This is commonly provided by rlang, but we define it locally to avoid
# adding rlang as a hard dependency just for this operator.
`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("mlr3", quietly = TRUE)) {
    return(invisible())
  }

  .paper1_register_learner(
    "classif.ratio_cnn",
    function() LearnerClassifRatioCNN$new()
  )
  .paper1_register_learner(
    "classif.clr_mlp",
    function() LearnerClassifCLRMLP$new()
  )
  .paper1_register_learner(
    "classif.codacore",
    function() LearnerClassifCoDaCoRe$new()
  )
  register_gof_filters()
  register_coda_feature_selection_filters()
  invisible()
}

# Suppress R CMD check notes for NSE and data.table symbols.
utils::globalVariables(
  c(
    ".",
    ".BY",
    ".N",
    ".SD",
    "..cols",
    "..feat_cols",
    "..feature_names",
    ".data",
    ".grouped_folds",
    ".stratified_folds",
    ":=",
    "AUC_re_pooled",
    "N",
    "accession",
    "anticoagulant",
    "anticoagulant_raw",
    "auc",
    "auc_ci_high",
    "auc_ci_low",
    "auc_se",
    "average_auc",
    "average_auc_ci_high",
    "average_auc_ci_low",
    "average_auc_se",
    "axis",
    "axis_label",
    "benchmark_component",
    "biofluid",
    "body_fluid",
    "cancer_type",
    "cancer_type_canonical",
    "cohort",
    "comparison_label",
    "comparison_type",
    "correlated",
    "dist_utopia",
    "dropout_loss",
    "epsilon",
    "extraction_kit",
    "extraction_kit_raw",
    "features",
    "fifelse",
    "frequency",
    "held_out",
    "k",
    "k_norm",
    "k_tb",
    "kit_biofluid_group",
    "kit_family",
    "kit_family_for_scoring",
    "kit_label",
    "kit_product",
    "kit_source",
    "kit_stratum_mean_auc",
    "label",
    "lambda",
    "leakage_flag",
    "learner_id",
    "library_kit",
    "library_kit_assignment_status",
    "lift",
    "lift_ci_hi",
    "lift_ci_lo",
    "lift_vs_rCLR_baseline",
    "max_cohort_atoms",
    "mean_k",
    "mean_metric",
    "metric_norm",
    "method",
    "method_label",
    "modality",
    "modality_family",
    "n_barycenter_atoms",
    "n_common_features",
    "n_features",
    "n_folds",
    "n_group_neg",
    "n_group_pos",
    "n_group_samples",
    "n_samples",
    "n_test",
    "n_test_neg",
    "n_test_pos",
    "n_train_samples",
    "neglog10p",
    "nogueira_index",
    "observed_auc",
    "pareto",
    "original_metric",
    "panel_size",
    "pars_norm",
    "perf_norm",
    "plausibility_flag",
    "present_in",
    "private",
    "provenance_block",
    "q_block_BH",
    "rCLR_baseline",
    "rclr_baseline_auc",
    "rclr_loco_auc",
    "re_auc",
    "re_ci_high",
    "re_ci_low",
    "re_se",
    "row_type",
    "sample_id",
    "score",
    "sd_metric",
    "se_metric",
    "selected",
    "selected_features",
    "self",
    "sinkhorn_backend",
    "sinkhorn_ot_auc",
    "smrnaseq_protocol",
    "stab_norm",
    "stability",
    "stability_tb",
    "status",
    "status_reason",
    "stratum",
    "subscore.id",
    "super",
    "technology",
    "track",
    "train_technologies",
    "truth",
    "var",
    "var_auc",
    "var_lift",
    "var_rclr_baseline",
    "variable",
    "worst_group_auc",
    "y"
  )
)

# cache test
