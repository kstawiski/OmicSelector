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
    "..cols",
    "..feature_names",
    ".data",
    ":=",
    "correlated",
    "dist_utopia",
    "dropout_loss",
    "features",
    "fifelse",
    "frequency",
    "k",
    "k_norm",
    "k_tb",
    "label",
    "learner_id",
    "mean_k",
    "mean_metric",
    "metric_norm",
    "n_features",
    "n_folds",
    "neglog10p",
    "nogueira_index",
    "pareto",
    "original_metric",
    "pars_norm",
    "perf_norm",
    "private",
    "score",
    "sd_metric",
    "se_metric",
    "selected",
    "selected_features",
    "self",
    "stab_norm",
    "stability",
    "stability_tb",
    "subscore.id",
    "super",
    "truth",
    "var",
    "variable"
  )
)

# cache test
