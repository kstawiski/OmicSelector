# Run Function in Parallel with Split-Aware Keys

Internal helper for parallel execution with proper handling of split
indices to prevent data leakage in cached computations.

## Usage

``` r
parallel_lapply(items, fun, split_id = NULL, ...)
```

## Arguments

- items:

  List of items to process

- fun:

  Function to apply to each item

- split_id:

  Current CV split identifier (for cache keys)

- ...:

  Additional arguments passed to fun

## Value

List of results
