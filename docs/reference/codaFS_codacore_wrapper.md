# CoDaCoRe-Derived Feature Selection

Selects features appearing in sparse CoDaCoRe balances. The official
\`codacore\` backend is used when available and functional; otherwise
the existing OmicSelector sparse-balance fallback is reused. Feature
scores are aggregated across balances using absolute slopes when
available.

## Usage

``` r
codaFS_codacore_wrapper(
  X,
  y,
  n_features = NULL,
  positive = NULL,
  input_scale = c("log", "raw"),
  backend = c("auto", "official", "fallback"),
  max_balances = NULL,
  top_k = 3L,
  class_weight = FALSE,
  verbose = FALSE,
  seed = 42L
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

- backend:

  One of \`"auto"\`, \`"official"\`, or \`"fallback"\`.

- max_balances:

  Maximum number of balances to extract.

- top_k:

  Initial sparse balance size for the fallback backend.

- class_weight:

  Logical; passed to the fallback backend.

- verbose:

  Logical; print backend messages.

- seed:

  Integer random seed.

## Value

A list with per-feature \`scores\`, \`selected_features\`, extracted
balances, and backend details.
