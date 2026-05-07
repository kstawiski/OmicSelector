# Create benchmark service from OmicPipeline

Convenience function to create a BenchmarkService

## Usage

``` r
omic_benchmark(pipeline, outer_folds = 5, inner_folds = 3, ...)
```

## Arguments

- pipeline:

  An OmicPipeline object

- outer_folds:

  Number of outer CV folds

- inner_folds:

  Number of inner CV folds

- ...:

  Additional arguments passed to BenchmarkService\$new()

## Value

A BenchmarkService object
