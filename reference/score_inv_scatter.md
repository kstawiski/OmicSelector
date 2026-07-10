# Score the 1D wavelet-scattering discriminator (forced row-by-row)

Scores each row of `X` independently with the frozen model from
[`fit_inv_scatter`](https://kstawiski.github.io/OmicSelector/reference/fit_inv_scatter.md).
The `ScatteringTorch1D` operator is REBUILT in python from the stored
config (no python pointer survives a fit), then each query is mapped to
the per-sample robust CLR over the frozen feature universe (absent
universe features carry the neutral rCLR value `0`), zero-padded to the
frozen `shape`, scattered in float64 ONE ROW AT A TIME, standardized
with the FROZEN training statistics, and passed through the FROZEN
ridge-logistic linear head. The score is \\\mathrm{intercept} + \sum_d
w_d \cdot \mathrm{std\\coeff}\_d(x)\\; larger = more case-like. Forcing
the \\n=1\\ scattering path uniformly (float64) makes a row's
coefficients – hence its score – bit-identical whether scored alone or
in any batch.

Queries with fewer than `model$hp$min_features` universe features
present (a column-overlap floor), or empty positive support over the
frozen universe on the (pre-rCLR) aligned ORIGINAL abundances
(`!any(X_use[i, ] > 0)`), return the neutral score `0`. A FLAT
all-equal-positive composition maps to the rCLR origin (all-zero signal)
but is a VALID specimen and is scored normally (NOT floored). The score
of a row depends only on that row and the frozen model, and is exactly
invariant to per-specimen positive scaling.

## Usage

``` r
score_inv_scatter(model, X, meta = NULL)
```

## Arguments

- model:

  An `inv_scatter_model` object from
  [`fit_inv_scatter`](https://kstawiski.github.io/OmicSelector/reference/fit_inv_scatter.md).

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

Anden J, Mallat S. (2014) Deep Scattering Spectrum. *IEEE TSP*
62(16):4114-4128.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_inv_scatter(X, y)
score_inv_scatter(model, X[1, , drop = FALSE])
} # }
```
