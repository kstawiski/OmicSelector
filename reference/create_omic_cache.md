# Create Split-Aware Cache

Creates a cache instance configured for OmicSelector's zero-leakage
requirements. Uses an LRU (Least Recently Used) eviction policy.

## Usage

``` r
create_omic_cache(max_size = 200 * 1024^2, max_age = 3600, dir = NULL)
```

## Arguments

- max_size:

  Maximum cache size in bytes. Default: 200MB.

- max_age:

  Maximum age of cached items in seconds. Default: 3600 (1 hour).

- dir:

  Optional directory for disk-based caching. NULL for memory-only.

## Value

A cachem cache object
