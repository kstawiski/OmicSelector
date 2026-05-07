# Within-Fold Median Imputation (Leakage-Free)

Imputes missing values using training-fold medians only. This prevents
the data leakage that occurs when medians are computed on the full
dataset (including test-fold samples) before the CV loop.

\*\*Common bug this prevents\*\*: Computing \`median(X\[, j\])\` on the
full matrix before the train/test split, which leaks test-fold
information into the imputed training set.

## Usage

``` r
impute_within_fold(X_train, X_test)
```

## Arguments

- X_train:

  Training fold feature matrix (samples x features)

- X_test:

  Test fold feature matrix (samples x features)

## Value

A list with:

- X_train:

  Imputed training matrix (training-only medians)

- X_test:

  Imputed test matrix (same training medians)

- impute_medians:

  Named vector of medians used for imputation

- n_imputed_train:

  Number of values imputed in training set

- n_imputed_test:

  Number of values imputed in test set

## Examples

``` r
if (FALSE) { # \dontrun{
for (fold in folds) {
  X_train <- feature_mat[train_idx, ]
  X_test  <- feature_mat[test_idx, ]

  # CORRECT: within-fold imputation
  imp <- impute_within_fold(X_train, X_test)
  X_train <- imp$X_train
  X_test  <- imp$X_test

  # WRONG (leaks!): imputing before the loop
  # for (j in 1:ncol(mat)) mat[is.na(mat[,j]), j] <- median(mat[,j], na.rm=TRUE)
}
} # }
```
