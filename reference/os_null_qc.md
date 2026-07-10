# Validate a Random-Panel Null Benchmark

Validate a Random-Panel Null Benchmark

## Usage

``` r
os_null_qc(null_result, min_draws = 1000L, expected_panel_size = NULL)
```

## Arguments

- null_result:

  Object from `os_panel_null_benchmark`.

- min_draws:

  Minimum accepted number of null draws.

- expected_panel_size:

  Expected locked-panel size.

## Value

A data frame with QC checks and pass/fail values.
