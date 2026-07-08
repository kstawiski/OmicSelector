# Create Pairwise Log-Ratio Image

Converts a single sample's expression vector into a pairwise log-ratio
image. For log-scale inputs, each pixel is the within-sample difference
\\x_i - x_j\\.

## Usage

``` r
encode_ratio_image(x)
```

## Arguments

- x:

  Numeric vector of feature values on log scale.

## Value

A square numeric matrix with one row and one column per feature.

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
x <- c(5, 3, 1)
encode_ratio_image(x)
#>      [,1] [,2] [,3]
#> [1,]    0    2    4
#> [2,]   -2    0    2
#> [3,]   -4   -2    0
```
