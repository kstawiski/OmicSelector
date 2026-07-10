# Adaptive Pseudocount log2 Transform

Compute `log2(x + pseudocount)` where the pseudocount is chosen
per-feature to avoid distorting the dynamic range. Critical when the
same analysis matrix mixes raw counts (typical scale 1-10000) with
pre-normalized ratios (typical scale 0.001-0.1) — adding a fixed +1 to a
0.003 ratio compresses its log2 to near zero and erases the signal.

## Usage

``` r
os_log_transform_adaptive(
  x,
  method = c("half_min", "fixed", "epsilon"),
  value = 1
)
```

## Arguments

- x:

  Numeric matrix or data.frame with samples in rows and features in
  columns. Must contain non-negative values; `NA` is propagated.

- method:

  One of `"half_min"`, `"fixed"`, or `"epsilon"`.

  - `"half_min"`: per-feature pseudocount = 0.5 \* smallest positive
    value in that column. Falls back to 1e-6 if a column has no positive
    values. Default.

  - `"fixed"`: same pseudocount across features, given by `value`
    (default 1.0). Backward-compatible with the `log2(x + 1)`
    convention.

  - `"epsilon"`: pseudocount = 1e-6 for all features.

- value:

  Numeric scalar; the pseudocount used when `method = "fixed"`. Ignored
  otherwise. Default 1.0.

## Value

A numeric matrix with the same dimensions and dimnames as `x`, carrying
an attribute `"pseudocount"` (length-`ncol(x)` vector) recording the
per-feature pseudocount used. The same pseudocount vector should be
reused on held-out data to avoid leakage.

## Examples

``` r
x <- matrix(c(0.001, 0.05, 0.13,
              2,     30,   1500), nrow = 2, byrow = TRUE)
colnames(x) <- c("ratio_a", "ratio_b", "count_c")
os_log_transform_adaptive(x, method = "half_min")
#>        ratio_a   ratio_b   count_c
#> [1,] -9.380822 -3.736966 -2.358454
#> [2,]  1.000361  4.908092 10.550809
#> attr(,"pseudocount")
#> ratio_a ratio_b count_c 
#>  0.0005  0.0250  0.0650 
os_log_transform_adaptive(x, method = "fixed", value = 1)
#>          ratio_a    ratio_b    count_c
#> [1,] 0.001441974 0.07038933  0.1763228
#> [2,] 1.584962501 4.95419631 10.5517083
#> attr(,"pseudocount")
#> ratio_a ratio_b count_c 
#>       1       1       1 
```
