# Convert qPCR Ct values to relative abundance with a fitted DA-cVAE's frozen reference

Applies a fitted `ct`-mode
[`fit_dacvae`](https://kstawiski.github.io/OmicSelector/reference/fit_dacvae.md)
model's FROZEN per-miRNA Ct reference (median train Ct), assay
efficiency \\E\\, and censored-Ct MZR floor to a matrix of qPCR
cycle-threshold values, returning the relative-abundance matrix
\\E^{-(Ct - \mathrm{ref})}\\ (censored / undetected Ct \\\to\\ the
floor). Exposed for the cross-platform harness to convert qPCR cohorts
to the abundance simplex with the SAME frozen reference the model was
fit under (single-sample-faithful). For an `abundance`-mode model this
is a no-op identity.

## Usage

``` r
dacvae_ct_to_abundance(model, X)
```

## Arguments

- model:

  A `dacvae_model` from
  [`fit_dacvae`](https://kstawiski.github.io/OmicSelector/reference/fit_dacvae.md).

- X:

  Numeric matrix of qPCR Ct values (samples x features); named feature
  columns.

## Value

Relative-abundance matrix with the same dimensions/dimnames as `X`.

## See also

[`fit_dacvae`](https://kstawiski.github.io/OmicSelector/reference/fit_dacvae.md),
[`score_dacvae`](https://kstawiski.github.io/OmicSelector/reference/score_dacvae.md)
