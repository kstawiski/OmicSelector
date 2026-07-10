# Fit kit-orthogonal ILR-like residualization.

The training CLR feature space is used to estimate the dominant kit axis
between the two most populated known strata. Scores are computed after
removing the projection onto that axis. If fewer than two strata are
available, scoring falls back to the global CLR panel sum and the fit
status records the limitation.

## Usage

``` r
fit_kit_orthogonal_ilr(
  expr_matrix,
  sample_meta,
  kit_label_col,
  biofluid_col = "biofluid",
  min_samples_per_stratum = 5L,
  pseudocount = NULL
)
```

## Arguments

- expr_matrix:

  Numeric matrix, samples x features, counts or non-negative abundance
  values.

- sample_meta:

  data.frame with one row per sample.

- kit_label_col:

  Column in \`sample_meta\` containing library-kit family.
  Missing/unknown kit values fall back to \`biofluid_col\`, then global.

- biofluid_col:

  Optional fallback column. Default "biofluid".

- min_samples_per_stratum:

  Minimum training samples needed for a stratum-specific centering set.

- pseudocount:

  Optional additive pseudocount before log transform.

## Value

Fit object consumed by \`score_kit_orthogonal_ilr()\`.
