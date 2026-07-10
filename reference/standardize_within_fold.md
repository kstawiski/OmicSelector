# Within-Fold Standardization (Leakage-Free)

Standardizes features using training-fold statistics only. Test-fold
features are centered and scaled using the same training-derived
parameters.

## Usage

``` r
standardize_within_fold(X_train, X_test)
```

## Arguments

- X_train:

  Training fold feature matrix

- X_test:

  Test fold feature matrix

## Value

A list with:

- X_train:

  Standardized training matrix

- X_test:

  Standardized test matrix (using training means/sds)

- means:

  Training-derived column means

- sds:

  Training-derived column SDs (zeros replaced with 1)
