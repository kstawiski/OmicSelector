# Limit-of-detection fallback imputation for qPCR non-detects

Replaces every NA in `expr_mat` with the supplied `lod_ct` value
(default 40). Used by the Paper-3 pipeline when Bayesian hierarchical
imputation via `nondetects` is unavailable or fails. Sets
`nondetect_method = "lod_fallback"` and `n_imputed` attributes.

## Usage

``` r
qpcr_nondetect_lod_fallback(expr_mat, lod_ct = 40, log_handler = message)
```

## Arguments

- expr_mat:

  Numeric matrix of Ct values; NAs represent non-detects.

- lod_ct:

  Numeric scalar; the Ct value used to impute NAs. Default 40.

- log_handler:

  Function used to emit progress messages. Default `message`.

## Value

Matrix of the same shape as `expr_mat`, with NAs filled by `lod_ct`.
Carries attributes `nondetect_method` and `n_imputed`.

## Examples

``` r
m <- matrix(c(20, NA, 25, 30, NA, NA), nrow = 2,
            dimnames = list(c("miR-1", "miR-2"), c("S1", "S2", "S3")))
qpcr_nondetect_lod_fallback(m, log_handler = function(...) invisible())
#> Error in qpcr_nondetect_lod_fallback(m, log_handler = function(...) invisible()): could not find function "qpcr_nondetect_lod_fallback"
```
