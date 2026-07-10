# Fit kit-residual MAD scoring.

Per feature, log-CPM values are robustly centered and scaled within
kit/biofluid strata using median and MAD. Unknown kit falls back to
biofluid, then global. This optional method is included as a negative or
robustness comparator.

## Usage

``` r
fit_kit_residual_mad(
  expr_matrix,
  sample_meta,
  kit_label_col,
  biofluid_col = "biofluid",
  min_samples_per_stratum = 5L,
  pseudocount = 0.5,
  mad_floor = 0.001
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

  Additive pseudocount before log-CPM transformation.

- mad_floor:

  Minimum robust scale used to avoid division by zero.

## Value

Fit object consumed by \`score_kit_residual_mad()\`.
