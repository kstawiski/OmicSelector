# Safe ROC Computation (Direction-Guarded)

Wrapper around [`pROC::roc()`](https://rdrr.io/pkg/pROC/man/roc.html)
that enforces `direction = "auto"` by default and warns loudly when
direction is hardcoded.

\*\*Common bug this prevents\*\*: Hardcoding `direction = "<"` when the
predictor's score direction depends on class weights, calibration, or
model type. With extreme class weighting (e.g., 17x), the positive class
can have LOWER predicted probabilities, inverting the AUC to 1 -
true_AUC.

## Usage

``` r
safe_roc(
  response,
  predictor,
  direction = "auto",
  levels = c(0, 1),
  quiet = TRUE,
  ...
)
```

## Arguments

- response:

  Binary outcome vector (0/1 or factor)

- predictor:

  Numeric prediction scores

- direction:

  ROC direction. Default `"auto"` (strongly recommended). If set to
  `"<"` or `">"`, a warning is issued.

- levels:

  Two-element vector specifying the levels of the response. Default
  `c(0, 1)`.

- quiet:

  Logical; suppress pROC messages. Default TRUE.

- ...:

  Additional arguments passed to
  [`pROC::roc()`](https://rdrr.io/pkg/pROC/man/roc.html)

## Value

A pROC roc object

## Examples

``` r
if (FALSE) { # \dontrun{
roc_obj <- safe_roc(y_true, y_pred)
auc_val <- as.numeric(pROC::auc(roc_obj))
} # }
```
