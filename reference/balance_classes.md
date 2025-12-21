# Create Balanced Training Set

Creates a balanced training set using synthetic data augmentation.
Designed for use within cross-validation folds.

## Usage

``` r
balance_classes(task, method = c("smote", "noise", "undersample"), ratio = 1)
```

## Arguments

- task:

  An mlr3 classification task

- method:

  Augmentation method: "smote", "noise", "undersample"

- ratio:

  Target ratio of minority to majority (default: 1.0)

## Value

A balanced task
