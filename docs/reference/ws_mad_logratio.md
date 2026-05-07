# Median-centred log-ratio with optional MAD scaling

Robust descriptive transform on log-abundances: subtracts the per-sample
median log-abundance (more robust to outlier features than
mean-centering as in CLR/rCLR) and optionally divides by the per-sample
median absolute deviation (MAD; limits leverage from individual
high-magnitude features).

This is NOT an Aitchison-valid simplex projection. It is a robust
descriptive transform on log-abundances, suitable as input to learners
that do not require the zero-sum compositional coordinate constraint.

## Usage

``` r
ws_mad_logratio(x, pseudocount = 0.5, scale_by_mad = TRUE)
```

## Arguments

- x:

  Samples x features matrix (rows = samples).

- pseudocount:

  Additive pseudocount before logging. Default 0.5.

- scale_by_mad:

  Logical. If `TRUE`, divides each sample by its per-sample MAD after
  median-centering. Default `TRUE`.

## Value

Same shape as `x`; the median-centred (MAD-scaled) log-abundance values,
with attributes `centering` and `scaling`.
