# Blondal hemolysis index (log miR-451a - log miR-23a-3p)

Computes the eLife-14 Blondal-style hemolysis index on a log-abundance
matrix as \\\log(\mathrm{miR\text{-}451a}) -
\log(\mathrm{miR\text{-}23a\text{-}3p})\\. Both anchor features must be
present (canonical names `hsa-miR-451a` / `hsa-miR-23a-3p`, or the
corresponding MIMAT accessions `MIMAT0001631` / `MIMAT0000078`);
otherwise `NULL` is returned and the caller should treat the cell as not
hemolysis-correctable. This is the convention the single-sample pipeline
uses to feed
[`fit_hemolysis_rr`](https://kstawiski.github.io/OmicSelector/reference/fit_hemolysis_rr.md)
/
[`apply_hemolysis_rr`](https://kstawiski.github.io/OmicSelector/reference/apply_hemolysis_rr.md).

## Usage

``` r
hemolysis_index_blondal(X_log)
```

## Arguments

- X_log:

  Numeric matrix (samples x features) on log-abundance scale. Column
  names must be canonical mature-miRNA names (post-resolver) or MIMAT
  accessions.

## Value

Numeric vector of length `nrow(X_log)`, or `NULL` if either anchor
feature is missing from `colnames(X_log)`.

## References

Blondal T, Jensby Nielsen S, Baker A, Andreasen D, Mouritzen P, Wrang
Teilum M, Dahlsveen IK. (2013) Assessing sample and miRNA profile
quality in serum and plasma or other biofluids. *Methods* 59(1): S1-S6.

## Examples

``` r
X_log <- matrix(c(8, 6, 4, 12, 5, 7, 3, 11), nrow = 2,
                dimnames = list(NULL,
                                c("hsa-miR-451a", "hsa-miR-23a-3p",
                                  "hsa-miR-1", "hsa-miR-2")))
hemolysis_index_blondal(X_log)
#> [1]  4 -6
```
