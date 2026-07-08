# Generate Split-Aware Cache Key

Creates a cache key that includes split information to prevent leakage.
This is the CRITICAL function for safe caching in nested CV.

## Usage

``` r
generate_cache_key(operation, split_id, features = NULL, params = NULL)
```

## Arguments

- operation:

  Name of the operation being cached (e.g., "filter_anova")

- split_id:

  Identifier for the current CV split (REQUIRED)

- features:

  Character vector of features involved

- params:

  Named list of parameters affecting the result

## Value

A character string cache key

## Details

The key structure: "operation::split_id::features_hash::params_hash"

Example: "filter_anova::outer3_inner2::a1b2c3d4::x5y6z7"
