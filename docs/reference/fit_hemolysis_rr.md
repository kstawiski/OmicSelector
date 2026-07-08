# Fit robust-regression hemolysis nuisance model from training controls

Fits per-feature Huber M-estimation robust regressions
([`MASS::rlm`](https://rdrr.io/pkg/MASS/man/rlm.html)) of log-expression
on a hemolysis proxy score (and optionally a platelet proxy). The
resulting intercept and slope coefficients are frozen at fit time. At
deployment, the predicted nuisance contribution is subtracted from
test-sample log-expressions via
[`apply_hemolysis_rr`](https://kstawiski.github.io/OmicSelector/reference/apply_hemolysis_rr.md).

Failed fits (convergence failures, NA coefficients) fall back to a
zero-slope model (intercept = feature mean) and are flagged in the
`stable` vector.

## Usage

``` r
fit_hemolysis_rr(
  X_train,
  hemo_score,
  platelet_score = NULL,
  method = c("huber", "rlm"),
  clip_max = 2
)
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features). Ideally a mix of clean
  and mildly contaminated training controls.

- hemo_score:

  Numeric vector of length `nrow(X_train)`. Hemolysis proxy score per
  training sample.

- platelet_score:

  Optional numeric vector of length `nrow(X_train)`. Platelet
  contamination proxy. Default `NULL` (no platelet correction).

- method:

  `"huber"` (default,
  [`MASS::rlm`](https://rdrr.io/pkg/MASS/man/rlm.html) with `psi.huber`)
  or `"rlm"` (alias).

- clip_max:

  Maximum absolute correction per feature in log-units. Default 2.
  Prevents over-correction from unstable slope estimates.

## Value

Object of class `hemolysis_rr_fit` with components: `betas` (features
\\\times\\ coefficients matrix), `stable`, `has_platelet`, `clip_max`,
`train_hemo_range`, `train_platelet_range`, `method`, `n_train`.

## References

Huber PJ. (1964) Robust Estimation of a Location Parameter. *The Annals
of Mathematical Statistics* 35(1): 73–101.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
X <- matrix(rlnorm(100 * 20, 5, 0.8), nrow = 100,
            dimnames = list(NULL, paste0("mir-", seq_len(20))))
hemo <- runif(100)
fit <- fit_hemolysis_rr(X, hemo)
fit$n_train
} # }
```
