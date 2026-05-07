# Safe AUC with Confidence Interval

Computes AUC with DeLong confidence interval using
[`safe_roc()`](https://kstawiski.github.io/OmicSelector/reference/safe_roc.md).
Returns a named list with AUC, CI, and the ROC object.

## Usage

``` r
safe_auc(response, predictor, ci_method = "delong", ...)
```

## Arguments

- response:

  Binary outcome vector (0/1)

- predictor:

  Numeric prediction scores

- ci_method:

  CI method: "delong" (default) or "bootstrap"

- ...:

  Additional arguments passed to
  [`safe_roc()`](https://kstawiski.github.io/OmicSelector/reference/safe_roc.md)

## Value

A list with:

- auc:

  Numeric AUC value

- ci_lower:

  Lower bound of 95% CI

- ci_upper:

  Upper bound of 95% CI

- roc:

  The pROC roc object

- direction:

  Direction used by pROC
