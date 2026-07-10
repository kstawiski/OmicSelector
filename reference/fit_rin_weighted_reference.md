# Fit a RIN-weighted frozen reference profile on training samples only.

Fit a RIN-weighted frozen reference profile on training samples only.

## Usage

``` r
fit_rin_weighted_reference(
  expr,
  meta,
  rin_col = NULL,
  cohort_col = "cohort",
  pseudocount = NULL
)
```

## Arguments

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- rin_col:

  Optional RNA-integrity-number column.

- cohort_col:

  Column identifying training cohorts.

- pseudocount:

  Optional additive pseudocount before log transform.
