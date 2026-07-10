# Predict technology-residualized ALR scores

Predict technology-residualized ALR scores

## Usage

``` r
predict_tech_residualized_alr(
  fit,
  expr,
  meta,
  technology_col = fit$technology_col
)
```

## Arguments

- fit:

  Fit from \`fit_tech_residualized_alr()\`.

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`; disease labels are not used.

- technology_col:

  Column identifying technology class.

## Value

Numeric score vector with pivot and model audit attributes.
