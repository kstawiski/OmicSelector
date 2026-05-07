# Simple Row-Major Grid Encoding

Places feature values in a rectangular grid in row-major order. Simple
reshape with zero-padding if (rows × cols) \> p.

## Usage

``` r
encode_simple_grid(x, rows = NULL, cols = NULL)
```

## Arguments

- x:

  Numeric vector of feature values (length p)

- rows:

  Number of rows in the grid (default: ceiling(sqrt(p)))

- cols:

  Number of columns (default: ceiling(p / rows))

## Value

A rows × cols numeric matrix

## References

Sharma A et al. (2019). DeepInsight: A methodology to transform a
non-image data to an image for convolution neural network architecture.
*Scientific Reports*, 9, 11399.

## Examples

``` r
encode_simple_grid(1:6, rows = 2, cols = 3)
#>      [,1] [,2] [,3]
#> [1,]    1    2    3
#> [2,]    4    5    6
```
