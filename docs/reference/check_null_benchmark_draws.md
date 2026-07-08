# Validate Random-Panel Null Benchmark

Checks whether the number of random draws is sufficient for a stable
percentile estimate. Warns when fewer than 500 draws are used.

## Usage

``` r
check_null_benchmark_draws(n_draws, percentile = 0.95)
```

## Arguments

- n_draws:

  Number of random panel draws performed

- percentile:

  The percentile being estimated (e.g., 0.95)

## Value

Invisible TRUE if adequate, FALSE with warning otherwise
