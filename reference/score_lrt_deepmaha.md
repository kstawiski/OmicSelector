# Score the deep-feature class-conditional Mahalanobis LRT discriminator

Scores each row of `X` independently with the frozen model from
[`fit_lrt_deepmaha`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_deepmaha.md),
in PURE R (no python). Each query is mapped to the per-sample robust CLR
over the frozen feature universe (absent universe features carry the
neutral rCLR value `0`), embedded by the exported MLP forward, and
scored by the class-conditional Mahalanobis log-likelihood ratio (case
minus control) with the frozen tied covariance. Larger = more case-like.

Queries with fewer than `model$hp$min_features` universe features
present (a column-overlap floor), or empty positive support over the
frozen universe on the (pre-rCLR) aligned abundances
(`!any(X_use[i, ] > 0)`), return the neutral score `0`. A FLAT
all-equal-positive composition maps to the rCLR origin but is a VALID
specimen and is scored normally. The score of a row depends only on that
row and the frozen model and is exactly invariant to per-specimen
positive scaling.

## Usage

``` r
score_lrt_deepmaha(model, X, meta = NULL)
```

## Arguments

- model:

  A `lrt_deepmaha_model` object from
  [`fit_lrt_deepmaha`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_deepmaha.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with named feature columns.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`; larger values are more
case-like.

## References

Lee K, Lee K, Lee H, Shin J. (2018) A Simple Unified Framework for
Detecting Out-of-Distribution Samples and Adversarial Attacks. *NeurIPS*
31. arXiv:1807.03888.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_lrt_deepmaha(X, y)
score_lrt_deepmaha(model, X[1, , drop = FALSE])
} # }
```
