# Score with provenance-block fixed-effect rCLR centering

Score with provenance-block fixed-effect rCLR centering

## Usage

``` r
score_block_fe_rclr(
  expr,
  meta,
  panel,
  cohort_col = "cohort",
  block_col = "provenance_block",
  fit = NULL,
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

  Column in \`meta\` identifying cohorts; accepted for a common scorer
  signature but not used by this method.

- block_col:

  Column in \`meta\` identifying provenance blocks.

- fit:

  Optional object from \`fit_block_fe_rclr()\`.

- pseudocount:

  Optional additive pseudocount before log transform.

## Value

Numeric score vector, one value per sample. Attribute
\`fallback_samples\` counts samples assigned to the training global
center because their block was absent from the fit.
