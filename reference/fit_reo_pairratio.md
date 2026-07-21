# Fit REO-pairratio within-sample log-ratio discriminator

Learns a frozen penalized-logistic discriminator over selected pairwise
log-ratios. Features are pre-screened using only the training data by
the absolute two-sample t-statistic on `log(X_train + pseudocount)`. All
pairwise log-ratios among the screened features are then generated with
[`ws_logratio`](https://kstawiski.github.io/OmicSelector/reference/ws_logratio.md),
and `glmnet::cv.glmnet()` fits a binomial elastic-net model with
deterministic folds. No test data or test-time batch statistics are
used.

A retained ratio \\\log(x_a + c) - \log(x_b + c)\\ is computed from one
specimen's two features, so a frozen linear predictor over these ratios
is single-sample deployable. Larger returned linear predictors are more
case-like.

## Usage

``` r
fit_reo_pairratio(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be uniquely named feature ids.

- y_train:

  Integer/numeric 0/1 labels aligned to `X_train`; 1 is case/disease and
  0 is control.

- meta_train:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

- hp:

  List of frozen hyperparameters: `pseudocount` (default `0.5`,
  positive), `m_features` (default `min(40L, ncol(X_train))`, at least
  2), `alpha` (default `1`, in `[0, 1]`), `nfolds` (default `10L`, at
  least 3), `lambda_rule` (default `"1se"`, or `"min"`), and `seed`
  (default `1L`, used only to make fold assignment reproducible while
  restoring the global RNG state).

## Value

A plain list of class `reo_pairratio_model` containing `selected_pairs`,
nonzero `coefficients`, `intercept`, `pseudocount`, `screened_features`,
`feature_universe`, `lambda`, `lambda_rule`, `fit_status`, and resolved
`hp`. If every training log-ratio is constant, or if the chosen lambda
yields an intercept-only glmnet model, `selected_pairs` has zero rows
and scoring returns the intercept for every specimen. The former state
is recorded as `"intercept_only_zero_variance"` and has no fitted
lambda.

## References

Aitchison J. (1986) *The Statistical Analysis of Compositional Data*.
Chapman and Hall.

Friedman J, Hastie T, Tibshirani R. (2010) Regularization paths for
generalized linear models via coordinate descent. *Journal of
Statistical Software* 33: 1-22.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(60 * 20, shape = 2), nrow = 60,
            dimnames = list(NULL, paste0("miR-", seq_len(20))))
y <- rep(c(0, 1), each = 30)
X[y == 1, "miR-1"] <- X[y == 1, "miR-1"] + 10
X[y == 0, "miR-2"] <- X[y == 0, "miR-2"] + 10
model <- fit_reo_pairratio(X, y, hp = list(m_features = 8L, nfolds = 5L))
score <- score_reo_pairratio(model, X)
} # }
```
