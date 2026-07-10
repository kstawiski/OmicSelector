# Score biofluid-residualized ALR

Score biofluid-residualized ALR

## Usage

``` r
score_biofluid_residualized_alr(
  expr,
  meta,
  panel,
  biofluid_col = "biofluid",
  fit = NULL,
  feature_weights = NULL,
  train_expr = NULL,
  train_meta = NULL,
  outcome_col = "disease",
  cohort_col = "cohort",
  pseudocount = NULL,
  ...
)
```

## Arguments

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- panel:

  Character vector of panel features.

- biofluid_col:

  Column in \`meta\` with canonical or raw biofluid.

- fit:

  Optional fit. If NULL, the fit is built on \`expr\` and \`meta\`.

- feature_weights:

  Optional named numeric panel weights.

- train_expr:

  Optional training matrix when \`fit\` is NULL.

- train_meta:

  Optional training metadata when \`fit\` is NULL.

- outcome_col:

  Binary disease label column in \`train_meta\`.

- cohort_col:

  Cohort column in \`train_meta\`.

- pseudocount:

  Optional additive pseudocount before log transform.

- ...:

  Additional arguments passed to \`fit_biofluid_stratified_rclr()\`.
