# SMOTE Augmentation for Omics Data

Applies SMOTE to balance class distributions in omics data. Uses
k-nearest neighbors to generate synthetic samples.

## Usage

``` r
smote_augment(task, ratio = 1, k = 5L)
```

## Arguments

- task:

  An mlr3 classification task

- ratio:

  Target ratio of minority to majority (default: 1.0 = balanced)

- k:

  Number of nearest neighbors for SMOTE (default: 5)

## Value

A new task with augmented data

## Details

For omics data, SMOTE should be applied with caution: - Use within CV
folds only (to prevent data leakage) - Consider feature selection before
SMOTE (faster, better interpolation) - May create unrealistic expression
profiles in high-dimensional space

## Examples

``` r
if (FALSE) { # \dontrun{
# Create balanced task
task_balanced <- smote_augment(task, ratio = 1.0)
} # }
```
