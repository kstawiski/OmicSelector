# Score the DA-cVAE single-sample novelty (out-of-distribution) capability

Returns, for each row of `X`, the mean Euclidean distance in the frozen
latent space from its encoded \\\mu\\ to its `model$novelty_k` nearest
TRAIN control (reference) latents — a kNN-distance novelty score (design
F9: a multi-reference score, not a single-centroid Mahalanobis). Larger
= more out-of-distribution relative to the training reference
population. This is a SECONDARY capability (it supports the
cross-platform / OOD use case) and is NOT the disease score. Floored /
degenerate specimens return `0`. Pure base-R, row-equivariant, exactly
scale-invariant.

## Usage

``` r
score_dacvae_novelty(model, X, meta = NULL)
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
  ignored.

## Value

Finite numeric vector of length `nrow(X)`; higher = more novel.

## See also

[`fit_dacvae`](https://kstawiski.github.io/OmicSelector/reference/fit_dacvae.md),
[`score_dacvae`](https://kstawiski.github.io/OmicSelector/reference/score_dacvae.md)
