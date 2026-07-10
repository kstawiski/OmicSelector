# Estimate a Provenance-Only Prediction Floor

Fits simple train-fold-only logistic provenance models to estimate
whether non-expression source variables alone can predict the outcome.
This is a lightweight package-level governance helper, not a replacement
for a full prespecified modeling script.

## Usage

``` r
os_provenance_floor_suite(
  data,
  outcome,
  provenance_predictors,
  group_id = NULL,
  positive = NULL,
  n_folds = 5L,
  seed = 1L
)
```

## Arguments

- data:

  Data frame containing outcome and provenance predictors.

- outcome:

  Column name for the binary outcome.

- provenance_predictors:

  Character vector of non-expression predictor columns.

- group_id:

  Optional column name or vector defining grouped folds.

- positive:

  Optional positive-class label.

- n_folds:

  Number of grouped folds.

- seed:

  Random seed.

## Value

An `os_provenance_floor` object.
