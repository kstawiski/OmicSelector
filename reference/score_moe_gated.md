# Score the self-gated mixture-of-experts single-sample discriminator

Scores each row of `X` independently with the frozen model from
[`fit_moe_gated`](https://kstawiski.github.io/OmicSelector/reference/fit_moe_gated.md),
in PURE R (no python). Each query is mapped to the per-sample robust CLR
over the frozen feature universe (absent universe features carry the
neutral rCLR value `0`), embedded by the exported LayerNorm encoder
forward to the pre-activation embedding \\z\\, and scored by the
self-gated mixture logit \\s(z) = \sum_e g_e(z)\\h_e(z)\\. Larger = more
case-like.

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
score_moe_gated(model, X, meta = NULL)
```

## Arguments

- model:

  A `moe_gated_model` object from
  [`fit_moe_gated`](https://kstawiski.github.io/OmicSelector/reference/fit_moe_gated.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with named feature columns.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method (the gate is self-driven).

## Value

Plain finite numeric vector of length `nrow(X)`; larger values are more
case-like.

## References

Jacobs RA, et al. (1991) Adaptive Mixtures of Local Experts. *Neural
Computation* 3(1):79-87. Shazeer N, et al. (2017) Outrageously Large
Neural Networks. *ICLR*. arXiv:1701.06538.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_moe_gated(X, y)
score_moe_gated(model, X[1, , drop = FALSE])
} # }
```
