# Fit a frozen monotone quantile calibrator from training data

Estimates per-feature empirical quantile grids on a training cohort and
freezes them for deployment. At apply time
([`apply_frozen_quantile`](https://kstawiski.github.io/OmicSelector/reference/apply_frozen_quantile.md)),
test-sample feature values are mapped to their empirical quantile under
the training distribution and then transformed back to the corresponding
training-scale value. Tail values beyond the training range are clamped
(intentional winsorisation, documented limitation).

## Usage

``` r
fit_frozen_quantile(
  x_train,
  n_quantiles = 1000L,
  na_action = c("median", "leave_na")
)
```

## Arguments

- x_train:

  Numeric matrix (samples \\\times\\ features). Must have column names.
  All values must be non-negative.

- n_quantiles:

  Integer — number of quantile grid levels (resolution). Default 1000.

- na_action:

  `"median"` (default, replaces NA test values with per-feature training
  median) or `"leave_na"`.

## Value

Object of class `frozen_quantile_fit` with components: `quantile_ref`
(matrix of quantile values), `quantile_levels`, `feature_names`,
`median_per_feat`, `n_quantiles`, `na_action`.

## References

Bolstad BM, et al. (2003) *Bioinformatics* 19: 185–193.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
x_train <- matrix(abs(rnorm(50 * 20, 100, 30)), nrow = 50,
                  dimnames = list(NULL, paste0("miR-", seq_len(20))))
fit <- fit_frozen_quantile(x_train, n_quantiles = 100L)
} # }
```
