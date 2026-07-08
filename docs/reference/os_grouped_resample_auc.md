# Estimate AUC over Grouped Resampling Folds

Estimate AUC over Grouped Resampling Folds

## Usage

``` r
os_grouped_resample_auc(y, score, folds)
```

## Arguments

- y:

  Binary outcome vector.

- score:

  Numeric score vector where larger values indicate the positive class.

- folds:

  List of integer test-row indices.

## Value

A list with pooled and per-fold AUC values.
