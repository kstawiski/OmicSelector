# Create Pairwise Ratio Image from Expression Vector

Alias for \[encode_ratio_image()\] retained for backwards compatibility.

## Usage

``` r
make_ratio_image(x)
```

## Arguments

- x:

  Numeric vector of feature values on log scale.

## Value

A square numeric matrix with one row and one column per feature.

## References

Aitchison J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman and Hall.

## Examples

``` r
x <- c(5, 3, 1)
make_ratio_image(x)
#>      [,1] [,2] [,3]
#> [1,]    0    2    4
#> [2,]   -2    0    2
#> [3,]   -4   -2    0
```
