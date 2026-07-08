# Fit a Top-k Oriented-Pair Panel Scorer

Fits a deterministic top-k oriented-pair scorer on a locked panel. Pair
orientation is selected from the training data so that larger vote
fractions indicate the positive class. This is a lightweight, auditable
k-TSP-style scorer for small biomarker panels; it is not a full
switchBox replacement.

## Usage

``` r
os_ktsp_fit(X, y, k = 11L)
```

## Arguments

- X:

  Numeric matrix or data frame with samples in rows and features in
  columns.

- y:

  Binary outcome vector. The second factor level is treated as positive.

- k:

  Number of oriented pairs to retain.

## Value

An object of class `"os_ktsp_model"`.
