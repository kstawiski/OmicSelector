# Score the VICReg self-supervised single-sample discriminator

Scores each row of `X` independently with the frozen model from
[`fit_ssl_vicreg`](https://kstawiski.github.io/OmicSelector/reference/fit_ssl_vicreg.md),
in PURE R (no python). Each query is mapped to the per-sample robust CLR
over the frozen feature universe (absent universe features carry the
neutral rCLR value `0`), embedded by the exported encoder forward, and
scored by the frozen linear-probe logit \\w^\top z + b\\. Larger = more
case-like. No augmentation and no RNG are used at scoring.

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
score_ssl_vicreg(model, X, meta = NULL)
```

## Arguments

- model:

  A `ssl_vicreg_model` object from
  [`fit_ssl_vicreg`](https://kstawiski.github.io/OmicSelector/reference/fit_ssl_vicreg.md).

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

Bardes A, Ponce J, LeCun Y. (2022) VICReg:
Variance-Invariance-Covariance Regularization for Self-Supervised
Learning. *ICLR*. arXiv:2105.04906.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_ssl_vicreg(X, y)
score_ssl_vicreg(model, X[1, , drop = FALSE])
} # }
```
