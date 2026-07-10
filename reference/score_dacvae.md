# Score the DA-cVAE single-sample disease discriminator

Scores each row of `X` independently with the frozen
[`fit_dacvae`](https://kstawiski.github.io/OmicSelector/reference/fit_dacvae.md)
model, in PURE R (no python, no platform/cohort input). Each query is
mapped to its standardised Stage-1 structural-ILR balances, embedded by
the exported encoder to the deterministic latent \\z = \mu\\, and scored
by the frozen disease logit \\w^\top\mu + b_0\\. Larger = more
case-like. Specimens with empty positive support over the frozen
universe, fewer than `hp$min_features` universe features present, or a
degenerate (no-balance) model return the neutral score `0`. The score of
a row depends only on that row and the frozen model and is exactly
invariant to per-specimen positive scaling
([`singlesample_assert_row_equivariant`](https://kstawiski.github.io/OmicSelector/reference/singlesample_assert_row_equivariant.md)).

## Usage

``` r
score_dacvae(model, X, meta = NULL)
```

## Arguments

- model:

  A `dacvae_model` from
  [`fit_dacvae`](https://kstawiski.github.io/OmicSelector/reference/fit_dacvae.md).

- X:

  Numeric matrix (samples x features) of non-negative abundances; named
  cols.

- meta:

  Optional per-sample metadata; accepted for interface uniformity,
  ignored (the platform adversary is FIT-only).

## Value

Finite numeric vector of length `nrow(X)`; higher = more case-like.

## See also

[`fit_dacvae`](https://kstawiski.github.io/OmicSelector/reference/fit_dacvae.md),
[`score_dacvae_novelty`](https://kstawiski.github.io/OmicSelector/reference/score_dacvae_novelty.md)
