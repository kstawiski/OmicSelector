# Fit a conformal anomaly detector from a Tier R healthy reference cohort

Fits a conformal anomaly detector using k-nearest-neighbour distance as
the non-conformity measure. A held-out calibration fraction of the
healthy reference cohort is used to calibrate the conformal p-value,
providing a distribution-free, finite-sample FPR guarantee at level
`alpha` (Vovk et al., 2005; Romano et al., 2019).

## Usage

``` r
fit_conformal_anomaly(
  X_ref,
  representation = identity,
  k = 10L,
  calibration_split = 0.3,
  seed = 42L
)
```

## Arguments

- X_ref:

  Numeric matrix (samples \\\times\\ features) of healthy reference
  samples.

- representation:

  Function to project samples to a latent space. Default `identity` (use
  raw features).

- k:

  Integer - number of nearest neighbours. Default 10.

- calibration_split:

  Fraction of `X_ref` held out for calibration. Default 0.30.

- seed:

  Integer - RNG seed for the train/calibration split. Default 42.

## Value

Object of class `conformal_anomaly_fit`.

## References

Vovk V, Gammerman A, Shafer G. (2005) Algorithmic Learning in a Random
World. Springer.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
X_ref <- matrix(rnorm(200 * 30), nrow = 200)
fit <- fit_conformal_anomaly(X_ref, k = 10L)
fit$n_cal
} # }
```
