# Apply a frozen specimen-local preprocessor

Apply a frozen specimen-local preprocessor

## Usage

``` r
apply_singlesample_preprocessor(preprocessor, X_raw)
```

## Arguments

- preprocessor:

  A `singlesample_preprocessor` returned by
  [`fit_singlesample_preprocessor`](https://kstawiski.github.io/OmicSelector/reference/fit_singlesample_preprocessor.md).

- X_raw:

  Raw specimens-by-features matrix/data frame or a single named vector.
  See
  [`fit_singlesample_preprocessor`](https://kstawiski.github.io/OmicSelector/reference/fit_singlesample_preprocessor.md)
  for conversion rules and error codes.

## Value

A list with numeric matrix `abundance` and logical matrix `detected`,
with matching dimensions and names. Columns follow the frozen order,
restricted to supplied features. If no frozen columns are supplied, both
matrices have zero columns.
