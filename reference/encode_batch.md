# Batch-Encode an Expression Matrix with Any Encoding

Convenience wrapper that applies an encoding method to every row of an
expression matrix, returning a 3D array suitable as CNN input.

## Usage

``` r
encode_batch(
  X,
  method = c("simple", "corr", "deepinsight", "ratio"),
  fit_params = NULL,
  ...
)
```

## Arguments

- X:

  Numeric matrix (n × p)

- method:

  One of "simple", "corr", "deepinsight", "ratio"

- fit_params:

  For "corr" or "deepinsight": the fitted order/layout. For "simple" or
  "ratio": ignored.

- ...:

  Additional arguments passed to the encoding function

## Value

A 3D array (n × height × width)

## References

Sharma A et al. (2019). DeepInsight: A methodology to transform a
non-image data to an image for convolution neural network architecture.
*Scientific Reports*, 9, 11399.

Aitchison J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman and Hall.

Quinn TP, Erb I, Richardson MF, Crowley TM. (2020). Understanding
sequencing data as compositions: an outlook and review.
*Bioinformatics*, 36(16), 4424-4432.

## Examples

``` r
X <- matrix(rnorm(28), nrow = 4, ncol = 7)
encode_batch(X, method = "simple", rows = 1, cols = 7)
#> , , 1
#> 
#>            [,1]
#> [1,]  0.7423952
#> [2,] -0.8742927
#> [3,] -2.0828138
#> [4,]  0.0937679
#> 
#> , , 2
#> 
#>              [,1]
#> [1,] -0.001819812
#> [2,] -0.013100652
#> [3,]  0.667968075
#> [4,] -0.013167703
#> 
#> , , 3
#> 
#>            [,1]
#> [1,]  0.7760470
#> [2,] -2.0107353
#> [3,] -1.1281805
#> [4,]  0.3487996
#> 
#> , , 4
#> 
#>            [,1]
#> [1,] -0.3528985
#> [2,]  0.9447746
#> [3,] -1.0047199
#> [4,]  0.7239027
#> 
#> , , 5
#> 
#>             [,1]
#> [1,] -0.66883274
#> [2,] -1.11304019
#> [3,] -0.34280519
#> [4,]  0.04977932
#> 
#> , , 6
#> 
#>            [,1]
#> [1,] -1.2276820
#> [2,] -0.7640059
#> [3,] -1.2461821
#> [4,]  1.0167740
#> 
#> , , 7
#> 
#>            [,1]
#> [1,]  0.7236126
#> [2,] -1.0325271
#> [3,]  0.5573460
#> [4,] -0.2555814
#> 
```
