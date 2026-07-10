# Centered Log-Ratio Transformation

Applies the centered log-ratio transform row-wise. For log-scale inputs
this subtracts the within-sample mean from each feature. For raw
positive inputs the function applies \`log(x + pseudocount)\` before
centering.

## Usage

``` r
clr_transform(x, input_scale = c("log", "raw"), pseudocount = 1e-08)
```

## Arguments

- x:

  Numeric vector or matrix with samples in rows and features in columns.

- input_scale:

  Either \`"log"\` (default) or \`"raw"\`.

- pseudocount:

  Small positive constant added to raw inputs before logging. Ignored
  when \`input_scale = "log"\`.

## Value

A numeric vector or matrix of the same dimensions as \`x\`.

## References

Aitchison J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman and Hall.

Egozcue JJ, Pawlowsky-Glahn V, Mateu-Figueras G, Barcelo-Vidal C.
(2003). Isometric logratio transformations for compositional data
analysis. *Mathematical Geology*, 35(3), 279-300.

Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. (2017). Microbiome
datasets are compositional: and this is not optional. *Frontiers in
Microbiology*, 8, 2224.

## Examples

``` r
clr_transform(c(5, 3, 1))
#> [1]  2  0 -2
clr_transform(matrix(c(5, 3, 1, 4, 2, 0), nrow = 2, byrow = TRUE))
#>      [,1] [,2] [,3]
#> [1,]    2    0   -2
#> [2,]    2    0   -2
```
