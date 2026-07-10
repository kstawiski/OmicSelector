# Fit the hemolysis plus kit-stable dual-anchor denominator.

Fit the hemolysis plus kit-stable dual-anchor denominator.

## Usage

``` r
fit_hemolysis_kit_dual_anchor(
  expr,
  meta,
  kit_col = NULL,
  cohort_col = "cohort",
  anchor_fit = NULL,
  pseudocount = NULL,
  ...
)
```

## Arguments

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- kit_col:

  Optional kit column.

- cohort_col:

  Column identifying training cohorts.

- anchor_fit:

  Optional precomputed \`kit_stable_anchor_fit\`.

- pseudocount:

  Optional additive pseudocount before log transform.

- ...:

  Additional arguments passed to \`fit_kit_stable_anchors()\`.
