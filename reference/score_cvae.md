# Score the counterfactual class-conditional VAE single-sample discriminator

Scores each row of `X` independently with the frozen model from
[`fit_cvae`](https://kstawiski.github.io/OmicSelector/reference/fit_cvae.md),
in PURE R (no python). Each query is mapped to the per-sample robust CLR
over the frozen feature universe (absent universe features carry the
neutral rCLR value `0`), then scored by the counterfactual evidence
log-likelihood ratio \\s(x) = \mathrm{ELBO}(x,\mathrm{case}) -
\mathrm{ELBO}(x,\mathrm{control})\\ using the latent MEAN \\z = \mu\\
(DETERMINISTIC – no sampling, no RNG). Larger = more case-like.

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
score_cvae(model, X, meta = NULL)
```

## Arguments

- model:

  A `cvae_model` object from
  [`fit_cvae`](https://kstawiski.github.io/OmicSelector/reference/fit_cvae.md).

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

Kingma DP, Welling M. (2014) Auto-Encoding Variational Bayes. *ICLR*.
arXiv:1312.6114.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_cvae(X, y)
score_cvae(model, X[1, , drop = FALSE])
} # }
```
