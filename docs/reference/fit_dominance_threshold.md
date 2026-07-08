# Calibrate a cohort-specific dominance threshold

Calibrates a cohort-specific threshold for `ws_dominance_flag` as the
(1 - alpha) quantile of `ws_dominance_score` computed on a training
reference cohort. Recommended over the universal 0.60 default, which is
provided only for backwards compatibility.

## Usage

``` r
fit_dominance_threshold(x_train, top_k = 2L, alpha = 0.05, pseudocount = 0.5)
```

## Arguments

- x_train:

  Samples x features matrix of the training reference cohort (typically
  Tier R, control samples only).

- top_k:

  Integer; passed to `ws_dominance_score`. Default 2.

- alpha:

  Numeric in (0, 1). Threshold is the (1 - alpha) quantile of the
  training dominance distribution, i.e. samples in the top `alpha`
  fraction of the training distribution are flagged. Default 0.05.

- pseudocount:

  Numeric; passed to `ws_dominance_score`. Default 0.5.

## Value

List with components `threshold`, `alpha`, `top_k`, `reference_n`,
`reference_quantiles`.
