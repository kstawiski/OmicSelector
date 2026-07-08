# Log-Contrast Lasso on CLR Features

Fits a lasso or elastic-net logistic model on centered log-ratio (CLR)
features, then projects the coefficient vector onto the zero-sum
hyperplane to obtain a valid log-contrast representation. Feature
importance is the absolute projected coefficient magnitude.

## Usage

``` r
codaFS_logcontrast_lasso(
  X,
  y,
  n_features = NULL,
  positive = NULL,
  input_scale = c("log", "raw"),
  alpha = 1,
  nfolds = 5L,
  nlambda = 100L,
  lambda_rule = c("lambda.1se", "lambda.min")
)
```

## Arguments

- X:

  Numeric matrix with samples in rows and features in columns.

- y:

  Binary outcome vector.

- n_features:

  Optional target panel size.

- positive:

  Optional positive-class label for non-numeric outcomes.

- input_scale:

  Either \`"log"\` (default) or \`"raw"\`.

- alpha:

  Elastic-net mixing parameter passed to \`glmnet\`.

- nfolds:

  Number of cross-validation folds used by \`cv.glmnet\`.

- nlambda:

  Number of lambda values used by \`cv.glmnet\`.

- lambda_rule:

  Either \`"lambda.1se"\` (default) or \`"lambda.min"\`.

## Value

A list with per-feature \`scores\`, the raw and zero-sum projected
coefficients, and \`selected_features\`.
