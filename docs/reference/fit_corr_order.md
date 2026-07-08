# Correlation-Ordered Grid Encoding

Reorders features by hierarchical clustering on feature-feature
correlation distance, then fills a grid in row-major order. Features
that co-vary end up as neighbors in the image, helping the CNN's local
receptive fields exploit correlation structure.

Must be FIT on training data only to avoid leakage.

## Usage

``` r
fit_corr_order(X_train, method = "average")
```

## Arguments

- X_train:

  Numeric matrix (n_train × p features) for fitting the order

- method:

  Linkage method for `hclust` (default: "average")

## Value

Integer vector of length p, the feature order to use for encoding

## References

Sharma A et al. (2019). DeepInsight: A methodology to transform a
non-image data to an image for convolution neural network architecture.
*Scientific Reports*, 9, 11399.

## Examples

``` r
X_train <- matrix(rnorm(30), nrow = 6, ncol = 5)
fit_corr_order(X_train)
#> [1] 3 1 2 4 5
```
