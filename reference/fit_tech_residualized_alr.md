# Fit technology-residualized ALR

Per feature, the training model is \`log(feature_count) ~ disease +
technology + (1\|cohort)\` when \`lme4\` is available and more than one
cohort exists. Prediction subtracts the learned technology fixed effect
and then computes an ALR-style panel score.

## Usage

``` r
fit_tech_residualized_alr(
  train_expr,
  train_meta,
  panel,
  technology_col = "technology",
  outcome_col = "disease",
  cohort_col = "cohort",
  pivot_features = .ta_default_pivot_pool(),
  anchor_features = character(),
  min_pivot_present = 3L,
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

- technology_col:

  Column identifying technology class.

- outcome_col:

  Binary disease label column in \`train_meta\`.

- cohort_col:

  Cohort column in \`train_meta\`.

- pivot_features:

  Candidate ALR pivot features.

- anchor_features:

  Optional train-selected fallback pivots.

- min_pivot_present:

  Minimum pivots required for prediction.

- pseudocount:

  Optional additive pseudocount before log transform.

## Value

Fit object consumed by \`predict_tech_residualized_alr()\`.
