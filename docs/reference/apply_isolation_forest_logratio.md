# Apply isolation forest to new samples and flag anomalies

Scores new samples using a fitted isolation forest
([`fit_isolation_forest_logratio`](https://kstawiski.github.io/OmicSelector/reference/fit_isolation_forest_logratio.md)),
computes empirical p-values relative to the training-score distribution,
and flags samples exceeding the training \\(1 - \texttt{fpr\\target})\\
quantile.

## Usage

``` r
apply_isolation_forest_logratio(fit, x_test, fpr_target = 0.05)
```

## Arguments

- fit:

  An `isolation_forest_logratio_fit` object.

- x_test:

  Numeric matrix (samples \\\times\\ features). All training features
  must be present.

- fpr_target:

  Numeric in (0, 1). Empirical FPR target for the flagging threshold.
  Default 0.05.

## Value

List with `score` (anomaly score in \\\[0,1\]\\; higher = more
anomalous), `pvalue_emp`, `threshold`, `flagged`.
