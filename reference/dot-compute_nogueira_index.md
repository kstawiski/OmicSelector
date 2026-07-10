# Compute Nogueira Stability Index (Internal Wrapper)

Internal wrapper that calls the canonical implementation from
stability.R. This ensures consistency and avoids code duplication.

## Usage

``` r
.compute_nogueira_index(fold_features, p)
```

## Arguments

- fold_features:

  List of character vectors, each containing features selected in one
  fold

- p:

  Total number of features in the dataset (the candidate universe)

## Value

Numeric stability index between 0 and 1

## References

Nogueira, S., Sechidis, K., & Brown, G. (2018). On the Stability of
Feature Selection Algorithms. Journal of Machine Learning Research,
18(174), 1-54.
