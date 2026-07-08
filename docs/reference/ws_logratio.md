# Within-Sample Pairwise Log-Ratio Features

Computes all pairwise log-ratios between features within each sample.
For p features, produces p\*(p-1)/2 ratio features. Each ratio is
inherently self-normalizing — multiplicative batch effects cancel.

If the input is already on log scale (e.g., log2 expression), the
log-ratio is simply the difference: log(A/B) = log(A) - log(B).

## Usage

``` r
ws_logratio(x, feature_names = NULL)
```

## Arguments

- x:

  A numeric matrix (samples x features). Assumed log-scale.

- feature_names:

  Optional character vector of feature names

## Value

Matrix with p\*(p-1)/2 columns (pairwise ratios)
