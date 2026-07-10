# Score with the mixed-effects scorer

If \`fit\` is supplied, this function predicts from the existing fit and
does not read disease labels from \`meta\`. If \`fit\` is NULL,
\`train_expr\`, \`train_meta\`, and \`outcome_col\` are required and the
model is fit on that training pool before prediction.

## Usage

``` r
score_mixed_effects_scorer(
  expr,
  meta,
  panel,
  cohort_col = "cohort",
  block_col = "provenance_block",
  fit = NULL,
  train_expr = NULL,
  train_meta = NULL,
  outcome_col = "disease",
  pseudocount = NULL
)
```

## Arguments

- expr:

  Numeric matrix, test samples x features.

- meta:

  Data frame aligned to \`expr\`.

- panel:

  Character vector of panel features.

- cohort_col:

  Cohort column.

- block_col:

  Provenance-block column.

- fit:

  Optional fit from \`fit_mixed_effects_scorer()\`.

- train_expr:

  Optional training expression matrix when \`fit\` is NULL.

- train_meta:

  Optional training metadata when \`fit\` is NULL.

- outcome_col:

  Binary disease label column in \`train_meta\`.

- pseudocount:

  Optional additive pseudocount before log transform.

## Value

Numeric score vector.
