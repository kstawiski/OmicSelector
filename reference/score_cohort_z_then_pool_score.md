# Score a panel after cohort-wise robust z-standardization

Within each cohort, all panel features are log-transformed, centered by
the cohort median, scaled by cohort MAD, and summed. When \`fit\` lacks
a new test cohort, the default is to compute unsupervised moments from
that test cohort; no disease labels are used for this fallback.

## Usage

``` r
score_cohort_z_then_pool_score(
  expr,
  meta,
  panel,
  cohort_col = "cohort",
  block_col = "provenance_block",
  fit = NULL,
  allow_new_group_moments = TRUE,
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

  Column in \`meta\` identifying provenance blocks; accepted for a
  common scorer signature but not used by this method.

- fit:

  Optional object from \`fit_cohort_z_then_pool_score()\`.

- allow_new_group_moments:

  Logical; if TRUE, new cohorts are centered on their own unlabeled
  samples.

- pseudocount:

  Optional additive pseudocount before log transform.

## Value

Numeric score vector, one value per sample.
