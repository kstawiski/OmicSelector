# Score the spectral-normalized neural Gaussian process discriminator (SNGP)

Scores each row of `X` independently with the frozen model from
[`fit_unc_sngp`](https://kstawiski.github.io/OmicSelector/reference/fit_unc_sngp.md),
in PURE R (no python). Each query is mapped to the per-sample robust CLR
over the frozen feature universe (absent universe features carry the
neutral rCLR value `0`), embedded by the exported spectral-normalized
LayerNorm encoder forward, projected by the fixed random-Fourier-feature
map, and scored by the mean-field-adjusted RFF-GP posterior logit.
Larger = more case-like.

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
score_unc_sngp(model, X, meta = NULL)
```

## Arguments

- model:

  A `unc_sngp_model` object from
  [`fit_unc_sngp`](https://kstawiski.github.io/OmicSelector/reference/fit_unc_sngp.md).

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

Liu JZ, Lin Z, Padhy S, Tran D, Bedrax-Weiss T, Lakshminarayanan B.
(2020) Simple and Principled Uncertainty Estimation with Deterministic
Deep Learning via Distance Awareness. *NeurIPS* 33. arXiv:2006.10108.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_unc_sngp(X, y)
score_unc_sngp(model, X[1, , drop = FALSE])
} # }
```
