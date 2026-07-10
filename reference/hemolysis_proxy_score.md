# Hemolysis Proxy Score

Computes a log-scale hemolysis proxy score for panel robustness checks.
By default this is \`miR-92a-3p - miR-320d\`.

## Usage

``` r
hemolysis_proxy_score(
  x,
  numerator = "MIMAT0000092",
  denominator = "MIMAT0006764",
  feature_names = NULL
)
```

## Arguments

- x:

  Numeric vector or matrix with samples in rows and features in columns.

- numerator:

  Name of the numerator feature. Defaults to \`"MIMAT0000092"\`.

- denominator:

  Name of the denominator feature. Defaults to \`"MIMAT0006764"\`.

- feature_names:

  Optional feature names. If omitted, \`names(x)\` or \`colnames(x)\`
  are used.

## Value

Numeric vector of per-sample proxy scores, or a scalar for vector input.
