# Apply Frozen ComBat Within Cross-Validation Folds

Properly applies FrozenComBat batch correction within each fold of
cross-validation. This is the CORRECT way to apply batch correction for
ML pipelines - it prevents data leakage by fitting parameters only on
training indices and applying to test.

## Usage

``` r
apply_frozen_combat_cv(
  data,
  batch,
  train_indices,
  test_indices = NULL,
  covariates = NULL,
  parametric = TRUE
)
```

## Arguments

- data:

  Full data matrix (samples x features)

- batch:

  Full batch vector

- train_indices:

  Row indices for training set

- test_indices:

  Row indices for test set (optional)

- covariates:

  Optional covariates data.frame (will be subset by indices)

- parametric:

  Use parametric empirical Bayes (default TRUE)

## Value

List with: - corrected_train: Batch-corrected training data -
corrected_test: Batch-corrected test data (if test_indices provided) -
frozen_combat: The fitted FrozenComBat object (for external validation)

## Details

\## IMPORTANT: Proper Usage in Nested CV

For nested cross-validation, you should use this function OR the PipeOp:

Option 1: Use a PipeOp in \`mlr3pipelines\`:
\`create_frozen_combat_pipeop(batch_col = "batch")\`.

Option 2: Fit within each custom CV fold, then transform the
corresponding test fold with the frozen parameters from that training
split.

\## WRONG: Do NOT do this! Do not apply ComBat to the full dataset
before splitting folds, because that leaks batch-adjustment information
from held-out samples into training.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
data <- matrix(rnorm(200), nrow = 40, ncol = 5)
batch <- rep(c("A", "B"), each = 20)

# 5-fold CV
folds <- split(1:40, rep(1:5, each = 8))

for (i in seq_along(folds)) {
  test_idx <- folds[[i]]
  train_idx <- setdiff(1:40, test_idx)

  result <- apply_frozen_combat_cv(
    data = data,
    batch = batch,
    train_indices = train_idx,
    test_indices = test_idx
  )

  # Train model on result$corrected_train
  # Evaluate on result$corrected_test
}
} # }
```
