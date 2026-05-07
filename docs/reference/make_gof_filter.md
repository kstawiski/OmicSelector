# Create GOF Filter

Convenience function to create GOF filters.

## Usage

``` r
make_gof_filter(type = c("ks", "hurdle", "zero_prop"), ...)
```

## Arguments

- type:

  Filter type: "ks", "hurdle", or "zero_prop"

- ...:

  Additional parameters passed to the filter

## Value

A Filter object

## Examples

``` r
if (FALSE) { # \dontrun{
filt <- make_gof_filter("ks")
filt <- make_gof_filter("hurdle", min_nonzero_per_class = 5)
} # }
```
