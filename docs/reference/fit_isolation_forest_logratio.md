# Fit a pure-R isolation forest on rCLR log-ratio inputs

Builds a pure-R isolation forest (Liu et al., 2008) on rCLR-transformed
compositional inputs. Anomalies have shorter expected path lengths to
isolation-tree leaves. The decision threshold is calibrated on training
scores to control the empirical false-positive rate at `fpr_target`
(applied at scoring time in
[`apply_isolation_forest_logratio`](https://kstawiski.github.io/OmicSelector/reference/apply_isolation_forest_logratio.md)).

This implementation uses no external dependencies and is intended for
robustness and auditability. For production use on large panels, a C++
backend (e.g., `isotree`) would be faster.

## Usage

``` r
fit_isolation_forest_logratio(
  x_train,
  n_trees = 100L,
  sample_size = 256L,
  max_depth = NULL,
  pseudocount = 0.5,
  seed = 42L
)
```

## Arguments

- x_train:

  Numeric matrix (samples \\\times\\ features).

- n_trees:

  Integer - number of isolation trees. Default 100.

- sample_size:

  Integer - subsample size per tree. Default 256.

- max_depth:

  Integer or `NULL`. Maximum tree depth. Default
  `ceiling(log2(min(sample_size, n)))`.

- pseudocount:

  Numeric. rCLR pseudocount. Default 0.5.

- seed:

  Integer - RNG seed. Default 42.

## Value

Object of class `isolation_forest_logratio_fit`.

## References

Liu FT, Ting KM, Zhou Z-H. (2008) Isolation Forest. *ICDM*.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
x_train <- matrix(abs(rnorm(80 * 20, 100, 30)), nrow = 80,
                  dimnames = list(NULL, paste0("miR-", seq_len(20))))
fit <- fit_isolation_forest_logratio(x_train, n_trees = 20L, sample_size = 32L)
} # }
```
