# Batch-Create Ratio Images from an Expression Matrix

Applies \[encode_ratio_image()\] row-wise to a log-scale expression
matrix.

## Usage

``` r
make_ratio_images(mat)
```

## Arguments

- mat:

  Numeric matrix with samples in rows and features in columns.

## Value

A 3D array with dimensions \`samples x features x features\`.

## References

Aitchison J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman and Hall.

## Examples

``` r
mat <- matrix(c(5, 3, 1, 4, 2, 0), nrow = 2, byrow = TRUE)
imgs <- make_ratio_images(mat)
dim(imgs)
#> [1] 2 3 3
```
