# Fit frozen logistic-normal empirical-Bayes prior from training data

Estimates per-feature CLR mean and variance on a training cohort. At
apply time
([`apply_logistic_normal_eb`](https://kstawiski.github.io/OmicSelector/reference/apply_logistic_normal_eb.md)),
each test-sample CLR observation is shrunk toward the per-feature prior
using a precision-weighted Bayesian update: \$\$\hat{z}\_{i,f} = w_f
z\_{i,f} + (1 - w_f) \mu_f,\$\$ where \\w_f = \sigma^2_f / (\sigma^2_f +
\sigma^2\_{\mathrm{obs},i})\\. The observation variance
\\\sigma^2\_{\mathrm{obs},i}\\ scales inversely with sample total counts
so shallow-coverage samples shrink harder. Per-feature variance is truly
load-bearing (not a scalar heuristic).

## Usage

``` r
fit_logistic_normal_eb(x_train, pseudocount = 0.5)
```

## Arguments

- x_train:

  Numeric matrix (samples \\\times\\ features). Must have column names.
  All values must be non-negative.

- pseudocount:

  Numeric \\\> 0\\. Added before log transform. Default 0.5.

## Value

Object of class `logistic_normal_eb_fit` with components:
`feature_names`, `pseudocount`, `prior_mean`, `prior_var`,
`median_total_train`.

## References

Aitchison J, Shen SM. (1980) *Biometrika* 67: 261–272. Efron B. (2010)
Large-Scale Inference. Cambridge UP.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
x_train <- matrix(abs(rnorm(50 * 20, 100, 30)), nrow = 50,
                  dimnames = list(NULL, paste0("miR-", seq_len(20))))
fit <- fit_logistic_normal_eb(x_train)
length(fit$prior_var)
} # }
```
