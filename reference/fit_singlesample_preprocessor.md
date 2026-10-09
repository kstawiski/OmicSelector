# Freeze specimen-local preprocessing for counts or Ct measurements

Records an ordered feature universe and conversion constants. No
training means, totals, quantiles or other estimated statistics are
retained.

## Usage

``` r
fit_singlesample_preprocessor(
  X_raw,
  type = c("counts", "ct"),
  ct_limit = 40,
  ct_base = 2,
  scale = c("cpm", "none"),
  features = colnames(X_raw)
)
```

## Arguments

- X_raw:

  Raw specimens-by-features matrix or data frame, with unique, non-empty
  column names. Counts must be numeric, finite and non-negative. Ct
  input may include character markers, which are converted to NA.

- type:

  Measurement type: `"counts"` or `"ct"`.

- ct_limit:

  Positive finite Ct detection limit (default 40).

- ct_base:

  Finite exponential base greater than 1 (default 2).

- scale:

  For counts, `"cpm"` divides each row by its own total over available
  frozen features and multiplies by one million; `"none"` leaves counts
  unchanged. Ignored for Ct conversion.

- features:

  Ordered character vector of unique feature names present in `X_raw`.
  Defaults to all training columns.

## Value

An S3 `singlesample_preprocessor` holding only `type`, `features`,
`ct_limit`, `ct_base`, and `scale`.

## Details

Counts of zero are non-detects. An all-zero row remains zero. Ct values
that are missing, non-finite, non-numeric, or at or above `ct_limit` are
non-detects. Their Ct is replaced by `ct_limit` before computing
`ct_base^(ct_limit - Ct)`, so they map to abundance 1. Negative finite
Ct values below the limit are converted by the same formula.

Columns outside the frozen universe are dropped. Missing frozen columns
remain absent, without imputation. For CPM, only available frozen
columns contribute to that specimen's total; missing columns can
therefore change its conversion. Scorer input guards determine whether
support is sufficient.

Invalid inputs raise `singlesample_input_error`. Codes include
`INVALID_INPUT`, `EMPTY_INPUT`, `INVALID_FEATURE_NAMES`,
`INVALID_PREPROCESSOR`, `INVALID_COUNT_VALUE` (negative counts),
`NONFINITE_REQUIRED_VALUE` (non-finite frozen count values), and
`NONFINITE_ABUNDANCE` (conversion or count-total overflow). Errors
involving measurement values report affected `rows`.

## See also

[`apply_singlesample_preprocessor`](https://kstawiski.github.io/OmicSelector/reference/apply_singlesample_preprocessor.md),
[`deploy_singlesample`](https://kstawiski.github.io/OmicSelector/reference/deploy_singlesample.md)

## Examples

``` r
counts <- rbind(c(a = 10, b = 30), c(a = 6, b = 2))
prep <- fit_singlesample_preprocessor(counts, type = "counts")
apply_singlesample_preprocessor(prep, counts[1, ])
#> $abundance
#>           a      b
#> [1,] 250000 750000
#> 
#> $detected
#>         a    b
#> [1,] TRUE TRUE
#> 
```
