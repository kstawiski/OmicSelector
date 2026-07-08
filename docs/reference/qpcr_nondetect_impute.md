# Bayesian hierarchical imputation of qPCR non-detects

Wraps `nondetects::qpcrImpute` for Ct-scale qPCR expression matrices.
Falls back to LOD imputation
([`qpcr_nondetect_lod_fallback`](https://kstawiski.github.io/OmicSelector/reference/qpcr_nondetect_lod_fallback.md))
when either `nondetects` or `HTqPCR` is unavailable, or when the
Bayesian fit fails. Sets `nondetect_method` and `n_imputed` attributes
on the returned matrix.

## Usage

``` r
qpcr_nondetect_impute(expr_mat, log_handler = message)
```

## Arguments

- expr_mat:

  Numeric matrix (features-by-samples) of Ct values; NAs represent
  non-detects.

- log_handler:

  Function used to emit progress messages. Default `message`.

## Value

Matrix of the same shape as `expr_mat` with non-detects imputed. Carries
attributes `nondetect_method` (one of `"nondetects_qpcrImpute"`,
`"lod_fallback"`, `"none_needed"`) and `n_imputed` (integer count of
imputed cells).

## Examples

``` r
if (FALSE) { # \dontrun{
m <- matrix(c(20, NA, 25, 30, NA, NA), nrow = 2,
            dimnames = list(c("miR-1", "miR-2"), c("S1", "S2", "S3")))
qpcr_nondetect_impute(m)
} # }
```
