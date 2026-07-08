# Cached Filter Computation

Wrapper for filter score computation with split-aware caching.
Automatically handles cache key generation and invalidation.

## Usage

``` r
cached_filter(filter_fun, task, filter_name, split_id, cache = NULL)
```

## Arguments

- filter_fun:

  Function that computes filter scores

- task:

  mlr3 Task object

- filter_name:

  Name of the filter (e.g., "anova", "mrmr")

- split_id:

  Current CV split identifier (REQUIRED)

- cache:

  Cache object from create_omic_cache()

## Value

Filter scores (from cache or freshly computed)
