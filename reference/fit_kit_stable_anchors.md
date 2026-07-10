# Fit kit-stable denominator anchors on training data only.

Lower kit_variance_ratio means lower between-kit signal relative to
within-kit noise. The returned anchor table records the training cohorts
and denominator used by downstream scoring. Evaluation data must not be
included.

## Usage

``` r
fit_kit_stable_anchors(
  expr,
  meta,
  kit_col = NULL,
  cohort_col = "cohort",
  anchor_fraction = 0.05,
  min_anchors = 8L,
  min_samples_per_kit = 5L,
  pseudocount = NULL
)
```

## Arguments

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- kit_col:

  Optional kit column. When NULL, a kit-like column is inferred.

- cohort_col:

  Column identifying training cohorts.

- anchor_fraction:

  Fraction of eligible features selected as anchors.

- min_anchors:

  Minimum number of anchors to select.

- min_samples_per_kit:

  Minimum samples per kit stratum.

- pseudocount:

  Optional additive pseudocount before log transform.
