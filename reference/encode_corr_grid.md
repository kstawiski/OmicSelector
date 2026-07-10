# Encode with Correlation-Ordered Grid

Encode with Correlation-Ordered Grid

## Usage

``` r
encode_corr_grid(x, order, rows = NULL, cols = NULL)
```

## Arguments

- x:

  Numeric feature vector (length p)

- order:

  Integer vector from `fit_corr_order`

- rows:

  Number of rows in the output grid. Defaults to
  `ceiling(sqrt(length(x)))`.

- cols:

  Number of columns in the output grid. Defaults to
  `ceiling(length(x) / rows)`.

## Value

A rows × cols numeric matrix

## References

Sharma A et al. (2019). DeepInsight: A methodology to transform a
non-image data to an image for convolution neural network architecture.
*Scientific Reports*, 9, 11399.

## Examples

``` r
encode_corr_grid(1:6, order = c(6, 5, 4, 3, 2, 1), rows = 2, cols = 3)
#>      [,1] [,2] [,3]
#> [1,]    6    5    4
#> [2,]    3    2    1
```
