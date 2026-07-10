# Panel-internal dominance score (within-sample QC)

Reports the fraction of total panel signal carried by the top-K features
in a single sample, on the raw abundance scale. Dimensionless,
scale-invariant, in \[K/N, 1\]. Sera with extreme dominance (e.g. top-2
features carrying \> 60 (massive haemolysis, incomplete extraction).

This is a SCALAR QC metric, not a normalization peer. It is designed for
use as an input to claim-gating in Module D, not as a Module-A
normalization method comparable to rCLR / ILR / ALR / MAD-logratio.

## Usage

``` r
ws_dominance_score(x, top_k = 2L, pseudocount = 0.5)
```

## Arguments

- x:

  Samples x features matrix.

- top_k:

  Integer; number of top features to include in the dominance numerator.
  Default 2.

- pseudocount:

  Additive pseudocount. Default 0.5.

## Value

Numeric vector of length `nrow(x)` with values in \[0, 1\].
