# Weighted Centered Log-Ratio Transformation

Computes a weighted within-sample CLR transform. The weighted center
uses down-weighted hemolysis-sensitive features, and the transformed
coordinates can optionally be shrunk by the same weights before
classification.

## Usage

``` r
weighted_clr_transform(
  x,
  weights = NULL,
  sensitivity = .paper1_hemolysis_coefficients(),
  input_scale = c("log", "raw"),
  pseudocount = 1e-08,
  shrink_features = TRUE,
  feature_names = NULL
)
```

## Arguments

- x:

  Numeric vector or matrix with samples in rows and features in columns.

- weights:

  Optional positive numeric vector of feature weights. If named, names
  are matched to \`feature_names\`.

- sensitivity:

  Optional named positive numeric vector of hemolysis coefficients. Used
  only when \`weights = NULL\`.

- input_scale:

  Either \`"log"\` (default) or \`"raw"\`.

- pseudocount:

  Small positive constant added before logging raw inputs.

- shrink_features:

  Logical; if \`TRUE\`, multiply the centered coordinates by their
  feature weights after centering.

- feature_names:

  Optional feature names.

## Value

Numeric vector or matrix of weighted CLR features.
