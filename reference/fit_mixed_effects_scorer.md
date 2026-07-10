# Fit disease weights from a cohort/provenance mixed-effects model

Fits one per-feature model on the training pool: \`log(feature_count) ~
disease + (1\|cohort) + (1\|provenance_block)\`. Terms with only one
level are omitted. If \`lme4\` cannot fit a feature, a fixed-effect
\`lm()\` fallback is attempted and recorded in the fit audit.

## Usage

``` r
fit_mixed_effects_scorer(
  train_expr,
  train_meta,
  panel,
  outcome_col = "disease",
  cohort_col = "cohort",
  block_col = "provenance_block",
  pseudocount = NULL
)
```

## Arguments

- train_expr:

  Numeric matrix, training samples x features.

- train_meta:

  Data frame aligned to \`train_expr\`.

- panel:

  Character vector of panel features.

- outcome_col:

  Binary disease label column in \`train_meta\`.

- cohort_col:

  Training cohort column.

- block_col:

  Training provenance-block column.

- pseudocount:

  Optional additive pseudocount before log transform.

## Value

A fit object consumed by \`predict_mixed_effects_scorer()\`.
