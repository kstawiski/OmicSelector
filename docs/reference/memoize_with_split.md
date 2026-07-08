# Memoize Function with Split Context

Creates a memoized version of a function that includes split_id in the
cache key. This is the recommended way to add caching to custom
functions.

## Usage

``` r
memoize_with_split(f, cache = NULL)
```

## Arguments

- f:

  Function to memoize

- cache:

  Cache object (optional, creates new if NULL)

## Value

A memoized function that requires split_id as first argument

## Details

The returned function has signature: function(split_id, ...) The
split_id is used in the cache key but passed through to f.

## Examples

``` r
if (FALSE) { # \dontrun{
# Original function
expensive_computation <- function(data, param) {
  # ... expensive operation ...
}

# Memoized version
cached_computation <- memoize_with_split(expensive_computation)

# Usage in CV loop
for (i in 1:5) {
  split_id <- paste0("fold_", i)
  result <- cached_computation(split_id, data = fold_data, param = 0.5)
}
} # }
```
