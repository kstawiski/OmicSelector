# Encode with DeepInsight Layout

Encode with DeepInsight Layout

## Usage

``` r
encode_deepinsight(x, layout)
```

## Arguments

- x:

  Numeric feature vector (length p)

- layout:

  Output of `fit_deepinsight`

## Value

A grid_size × grid_size numeric matrix

## References

Sharma A et al. (2019). DeepInsight: A methodology to transform a
non-image data to an image for convolution neural network architecture.
*Scientific Reports*, 9, 11399.

## Examples

``` r
layout <- list(positions = matrix(c(1, 1, 1, 2, 2, 1), ncol = 2, byrow = TRUE), grid_size = 2L)
encode_deepinsight(1:3, layout)
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    3    0
```
