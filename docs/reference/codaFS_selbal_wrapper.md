# Selbal-Style Forward Balance Feature Selection

Wraps the official \`selbal\` implementation when it is available. In
environments where \`selbal\` is not installed, or when \`backend =
"fallback"\` is requested, the function uses a native forward
balance-selection routine that greedily improves univariate logistic
deviance on training data only.

## Usage

``` r
codaFS_selbal_wrapper(
  X,
  y,
  n_features = NULL,
  positive = NULL,
  input_scale = c("log", "raw"),
  backend = c("auto", "official", "fallback"),
  max_terms = NULL,
  min_improvement = 1e-04,
  max_pair_features = 80L,
  prescreen_metric = c("clr_effect", "variance")
)
```

## Arguments

- X:

  Numeric matrix with samples in rows and features in columns.

- y:

  Binary outcome vector.

- n_features:

  Optional target panel size. Used as the default \`max_terms\`.

- positive:

  Optional positive-class label for non-numeric outcomes.

- input_scale:

  Either \`"log"\` (default) or \`"raw"\`.

- backend:

  One of \`"auto"\`, \`"official"\`, or \`"fallback"\`.

- max_terms:

  Maximum number of features in the selected balance. If \`NULL\`,
  defaults to \`n_features\`.

- min_improvement:

  Minimum deviance improvement required to add a new feature during the
  fallback forward search.

- max_pair_features:

  Maximum number of prescreened features considered by the forward
  search.

- prescreen_metric:

  Either \`"clr_effect"\` or \`"variance"\`.

## Value

A list with per-feature \`scores\`, \`selected_features\`, balance
membership, and backend details.
