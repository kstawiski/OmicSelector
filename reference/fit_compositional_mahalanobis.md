# Fit robust Mahalanobis detector on compositional log-ratio coordinates

Fits a Minimum Covariance Determinant (MCD) robust Mahalanobis distance
model on ILR (default) or rCLR log-ratio coordinates of a training
compositional matrix. The MCD estimator is provided by the `robustbase`
package. The function fails closed if `robustbase` is unavailable and
`require_robust = TRUE` (the default).

## Usage

``` r
fit_compositional_mahalanobis(
  x_train,
  transform = c("ilr", "rclr"),
  alpha = 0.75,
  pseudocount = 0.5,
  require_robust = TRUE
)
```

## Arguments

- x_train:

  Numeric matrix (samples \\\times\\ features). All values must be
  non-negative.

- transform:

  `"ilr"` (default, \\D \to D-1\\ dimensional, full rank under healthy
  assumptions) or `"rclr"` (sum-to-zero constraint, effective rank
  \\D-1\\).

- alpha:

  Numeric in (0.5, 1). MCD coverage fraction. Default 0.75.

- pseudocount:

  Numeric \\\> 0\\. Added before log transform. Default 0.5.

- require_robust:

  Logical. If `TRUE` (default), stops when `robustbase` is unavailable
  rather than silently downgrading to classical covariance.

## Value

Object of class `compositional_mahalanobis_fit`.

## References

Rousseeuw PJ, Van Driessen K. (1999) *Technometrics* 41: 212-223.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
x_train <- matrix(abs(rnorm(80 * 20, 100, 30)), nrow = 80,
                  dimnames = list(NULL, paste0("miR-", seq_len(20))))
fit <- fit_compositional_mahalanobis(x_train, transform = "rclr",
                                     require_robust = FALSE)
} # }
```
