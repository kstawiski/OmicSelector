# Detect Global Pre-Loop Imputation (Leakage Check)

Checks whether a feature matrix has already been imputed (no NAs
remaining) when NAs are expected. Call this at the top of a CV loop to
detect accidental global imputation.

## Usage

``` r
check_no_premature_imputation(X, expected_na_rate = 0, label = "X")
```

## Arguments

- X:

  Feature matrix to check

- expected_na_rate:

  Expected proportion of NAs (0 to 1). If the matrix has fewer NAs than
  expected, a warning is issued.

- label:

  Label for the matrix in warning messages (e.g., "panel_mat")

## Value

Invisible TRUE if check passes, FALSE with warning if suspicious
