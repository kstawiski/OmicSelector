# Compute Nogueira Stability Index

Computes the Nogueira Stability Index for feature selection consistency.
The index ranges from 0 (completely unstable) to 1 (perfectly stable).

## Usage

``` r
.compute_nogueira_index(fold_features, p)
```

## Arguments

- fold_features:

  List of character vectors, each containing features selected in one
  fold

- p:

  Total number of features in the dataset

## Value

Numeric stability index between 0 and 1

## Details

The Nogueira Stability Index is defined as: SI = 1 - (observed_variance
/ max_variance)

Where observed variance is the average pairwise disagreement between
feature sets, and max variance is what we'd expect from random
selection.

## References

Nogueira, S., Sechidis, K., & Brown, G. (2018). On the Stability of
Feature Selection Algorithms. Journal of Machine Learning Research,
18(174), 1-54.
