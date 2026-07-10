# Score baseline robust CLR for provenance-aware benchmark comparisons

Score baseline robust CLR for provenance-aware benchmark comparisons

## Usage

``` r
score_baseline_rclr(
  expr,
  meta,
  panel,
  cohort_col = "cohort",
  block_col = "provenance_block",
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

- cohort_col:

  Column in \`meta\` identifying cohorts.

- block_col:

  Column in \`meta\` identifying provenance blocks.

- pseudocount:

  Optional additive pseudocount before log transform.

## Value

Numeric score vector.
