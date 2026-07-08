# Matched Random-Panel Null Benchmark

Benchmarks a locked panel statistic against random panels of the same
size. The scorer receives `X` subset to a feature set and must return
one score per sample. Null panels are sampled without replacement and
are not allowed to overlap the locked panel.

## Usage

``` r
os_panel_null_benchmark(
  X,
  y,
  panel_features,
  scorer,
  n_draws = 1000L,
  feature_pool = NULL,
  seed = 1L,
  metric = "auc",
  ...
)
```

## Arguments

- X:

  Numeric matrix or data frame with named columns.

- y:

  Binary outcome vector.

- panel_features:

  Character vector of locked-panel feature names.

- scorer:

  Function with signature `scorer(X_sub, features)`.

- n_draws:

  Number of random panels.

- feature_pool:

  Candidate null features. Defaults to all non-panel columns.

- seed:

  Random seed.

- metric:

  Currently only `"auc"`.

- ...:

  Additional arguments passed to `scorer`.

## Value

Object of class `"os_panel_null"`.
