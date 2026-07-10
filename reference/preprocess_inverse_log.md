# Inverse-log preprocessor for pre-log-transformed microarray deposits

Converts a matrix of log-transformed expression intensities back to a
positive abundance-like scale via \\base^x\\. This is required when a
public deposit has already applied a \\\log_2\\ (or \\\log\_{10}\\)
transform before archival (common in Toray 3D-Gene and some Agilent
deposits). The resulting matrix is on a positive scale comparable to raw
NGS counts and can be fed directly to within-sample Module A methods
(rCLR, ILR, ALR-pivot, MAD log-ratio, dominance score) which require
positive compositional input.

Non-finite values (Inf, -Inf, NaN) are replaced with the per-row median
of finite values, with a warning.

## Usage

``` r
preprocess_inverse_log(x, base = 2, pseudo = 0)
```

## Arguments

- x:

  Numeric matrix (features \\\times\\ samples, or samples \\\times\\
  features — the function operates element-wise). Both orientations are
  accepted; the transform is purely element-wise.

- base:

  Numeric. The log base used by the deposit. Default 2 (Toray 3D-Gene
  standard). Use 10 for log10 deposits.

- pseudo:

  Numeric. Small offset added to the output to keep strictly positive
  values. Default 0 (output of \\base^x\\ is always \\\> 0\\ for finite
  inputs).

## Value

Numeric matrix of the same dimensions as `x`, on a positive abundance
scale. The attribute `preprocessing` records the transform applied.

## References

Toray Industries, Inc. (3D-Gene microRNA chip) — data format
specification.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
log2_matrix <- matrix(rnorm(50 * 20), nrow = 50)
colnames(log2_matrix) <- paste0("S", seq_len(20))
rownames(log2_matrix) <- paste0("miR-", seq_len(50))
abundance <- preprocess_inverse_log(log2_matrix, base = 2)
stopifnot(all(abundance > 0))
} # }
```
