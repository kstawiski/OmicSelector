# Score technology-residualized ALR

Score technology-residualized ALR

## Usage

``` r
score_tech_residualized_alr(
  expr,
  meta,
  panel,
  technology_col = "technology",
  fit = NULL,
  train_expr = NULL,
  train_meta = NULL,
  outcome_col = "disease",
  cohort_col = "cohort",
  pivot_features = .ta_default_pivot_pool(),
  anchor_features = character(),
  min_pivot_present = 3L,
  pseudocount = NULL
)
```

## Arguments

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- panel:

  Character vector of panel features.

- technology_col:

  Column identifying technology class.

- fit:

  Optional fit from \`fit_tech_stratified_rclr()\`.

- train_expr:

  Optional training matrix when \`fit\` is NULL.

- train_meta:

  Optional training metadata when \`fit\` is NULL.

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

Numeric score vector.
