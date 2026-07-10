# Score REO-pairratio within-sample log-ratio discriminator

Scores each row of `X` independently with the frozen pairwise log-ratio
coefficients learned by
[`fit_reo_pairratio`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_pairratio.md).
For one specimen and each stored pair `(a, b)`, if both features are
present in `X`, the term is
`log(X[i, a] + pseudocount) - log(X[i, b] + pseudocount)`. If either
feature is absent, that pair is dropped for that specimen. The returned
score is the frozen linear predictor, not
[`plogis()`](https://rdrr.io/r/stats/Logistic.html) of it.

If the model is intercept-only, or if none of the stored pairs is usable
for a scored specimen, the score is the intercept. Scoring uses no
random numbers and no scored-batch statistics.

## Usage

``` r
score_reo_pairratio(model, X, meta = NULL)
```

## Arguments

- model:

  A `reo_pairratio_model` object returned by
  [`fit_reo_pairratio`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_pairratio.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Numeric vector of one finite linear-predictor score per row of `X`.
Larger scores are more case-like.

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
score_reo_pairratio(model, X[1, , drop = FALSE])
} # }
```
