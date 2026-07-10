# Score the DANN domain-adversarial cross-cohort single-sample discriminator

Scores each row of `X` independently with the frozen model from
[`fit_dann`](https://kstawiski.github.io/OmicSelector/reference/fit_dann.md),
in PURE R (no python). Each query is mapped to the per-sample robust CLR
over the frozen feature universe (absent universe features carry the
neutral rCLR value `0`), embedded by the exported encoder forward, and
scored by the frozen linear LABEL logit \\w^\top z + b\\. Larger = more
case-like. The domain discriminator is NOT part of scoring (FIT-only
adversary); no cohort/kit input is used and no RNG is consumed at
scoring.

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
score_dann(model, X, meta = NULL)
```

## Arguments

- model:

  A `dann_model` object from
  [`fit_dann`](https://kstawiski.github.io/OmicSelector/reference/fit_dann.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with named feature columns.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method (the transfer aspect is in the fit, not the
  score).

## Value

Plain finite numeric vector of length `nrow(X)`; larger values are more
case-like.

## References

Ganin Y, et al. (2016) Domain-Adversarial Training of Neural Networks.
*JMLR* 17(59):1–35. arXiv:1505.07818.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_dann(X, y, meta_train = meta)
score_dann(model, X[1, , drop = FALSE])
} # }
```
