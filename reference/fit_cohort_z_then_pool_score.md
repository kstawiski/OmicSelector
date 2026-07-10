# Fit cohort-wise robust z-score moments

Fit cohort-wise robust z-score moments

## Usage

``` r
fit_cohort_z_then_pool_score(
  expr,
  meta,
  panel,
  cohort_col = "cohort",
  pseudocount = NULL
)
```

## Arguments

- expr:

  Numeric matrix, samples x features, training or scoring pool.

- meta:

  Data frame aligned to \`expr\`.

- panel:

  Character vector of panel features.

- cohort_col:

  Column in \`meta\` identifying cohorts.

- pseudocount:

  Optional additive pseudocount before log transform.

## Value

A fit object consumed by \`predict_cohort_z_then_pool_score()\`.
