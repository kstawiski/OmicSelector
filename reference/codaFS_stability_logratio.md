# Stability Selection on Pairwise Log-Ratios

Applies subsampling-based stability selection to pairwise log-ratio
features using elastic-net logistic models. Stable ratios are aggregated
back to individual features by their incidence-weighted selection
frequency.

## Usage

``` r
codaFS_stability_logratio(
  X,
  y,
  n_features = NULL,
  positive = NULL,
  input_scale = c("log", "raw"),
  alpha = 0.9,
  subsample_frac = 0.7,
  n_subsamples = 20L,
  nfolds = 3L,
  nlambda = 60L,
  lambda_rule = c("lambda.1se", "lambda.min"),
  selection_threshold = 0.6,
  max_pair_features = 60L,
  prescreen_metric = c("clr_effect", "variance")
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

  Elastic-net mixing parameter.

- subsample_frac:

  Fraction of samples used in each stability iteration.

- n_subsamples:

  Number of stability iterations.

- nfolds:

  Number of cross-validation folds within each subsample.

- nlambda:

  Number of lambda values used by \`cv.glmnet\`.

- lambda_rule:

  Either \`"lambda.1se"\` (default) or \`"lambda.min"\`.

- selection_threshold:

  Frequency threshold used to define stable ratios.

- max_pair_features:

  Maximum number of prescreened features retained before constructing
  pairwise ratios.

- prescreen_metric:

  Either \`"clr_effect"\` or \`"variance"\`.

## Value

A list with per-feature \`scores\`, stable pairwise log-ratio
frequencies, and \`selected_features\`.
