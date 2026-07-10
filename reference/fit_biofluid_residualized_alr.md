# Fit biofluid-residualized ALR

Fits one per-feature training model, \`log(feature_count) ~ disease +
biofluid + (1 \| cohort)\`, subtracts the estimated biofluid fixed
effect at prediction time, and scores the panel as an ALR against a
frozen pivot pool. If \`lme4\` is unavailable or a feature model fails,
an \`lm()\` fallback is used and recorded.

## Usage

``` r
fit_biofluid_residualized_alr(
  train_expr,
  train_meta,
  panel,
  biofluid_col = "biofluid",
  outcome_col = "disease",
  cohort_col = "cohort",
  pivot_features = NULL,
  pseudocount = NULL,
  min_pivot_present = 3L
)
```

## Arguments

- train_expr:

  Numeric matrix, training samples x features.

- train_meta:

  Data frame aligned to \`train_expr\`.

- panel:

  Character vector of panel features.

- biofluid_col:

  Column identifying biofluid class.

- outcome_col:

  Binary disease label column in \`train_meta\`.

- cohort_col:

  Cohort column in \`train_meta\`.

- pivot_features:

  Optional ALR pivot features.

- pseudocount:

  Optional additive pseudocount before log transform.

- min_pivot_present:

  Minimum pivot features required for prediction.
