# Score Mondrian-conformal class-conditional LRT discriminator

Scores each row of `X` independently with a frozen
[`fit_conf_mondrian`](https://kstawiski.github.io/OmicSelector/reference/fit_conf_mondrian.md)
model. A specimen is transformed to within-sample rCLR coordinates over
the model features present in `X`; its class-specific nonconformity
scores are compared with the frozen class-conditional calibration
scores; and the returned score is \$\$\log p\_{case}(x) - \log
p\_{control}(x).\$\$

## Usage

``` r
score_conf_mondrian(model, X, meta = NULL)
```

## Arguments

- model:

  A `conf_mondrian_model` object returned by
  [`fit_conf_mondrian`](https://kstawiski.github.io/OmicSelector/reference/fit_conf_mondrian.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. Larger values are more
case-like. If fewer than `model$hp$min_features` model features are
present, the documented neutral score `0` is returned for every row.

## Details

At scoring time the present feature set is
`intersect(model$feature_universe, colnames(X))`. Class centroids,
diagonal variances, and calibration nonconformity scores are re-derived
over that present subset from the frozen raw anchors only. Each row's
score then uses only that row's own rCLR vector plus the frozen
model-derived representation, so no scored-batch statistic or cross-row
coupling is used.

## References

Vovk V, Gammerman A, Shafer G. (2005) *Algorithmic Learning in a Random
World*. Springer.

Shafer G, Vovk V. (2008) A tutorial on conformal prediction. *Journal of
Machine Learning Research* 9: 371-421.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(60 * 30, shape = 2), nrow = 60,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 30)
X[y == 1, 1:5] <- X[y == 1, 1:5] + 20
X[y == 0, 6:10] <- X[y == 0, 6:10] + 20
model <- fit_conf_mondrian(X, y)
score_conf_mondrian(model, X[1, , drop = FALSE])
} # }
```
