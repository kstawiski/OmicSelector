# Define the null-coalescing operator %||%
# Returns lhs if not NULL, otherwise rhs
# This is commonly provided by rlang, but we define it locally to avoid
# adding rlang as a hard dependency just for this operator.
`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

# Suppress R CMD check notes for NSE and data.table symbols.
utils::globalVariables(
  c(
    ".",
    "..cols",
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
    "var",
    "variable"
  )
)

# cache test
