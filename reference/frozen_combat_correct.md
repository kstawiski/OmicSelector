# Convenience function for frozen ComBat correction

Simple wrapper for common use case: fit on training, apply to training
and test.

## Usage

``` r
frozen_combat_correct(
  train_data,
  train_batch,
  test_data = NULL,
  test_batch = NULL,
  covariates = NULL,
  parametric = TRUE
)
```

## Arguments

- train_data:

  Training data matrix (samples x features)

- train_batch:

  Training batch vector

- test_data:

  Optional test data matrix

- test_batch:

  Optional test batch vector (required if test_data provided)

- covariates:

  Optional covariates data.frame

- parametric:

  Use parametric priors (default TRUE)

## Value

List with corrected_train and optionally corrected_test

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate example data
set.seed(42)
train <- matrix(rnorm(100), nrow = 20, ncol = 5)
test <- matrix(rnorm(50), nrow = 10, ncol = 5)
train_batch <- rep(c("A", "B"), each = 10)
test_batch <- rep(c("A", "B"), each = 5)

# Apply frozen ComBat
result <- frozen_combat_correct(
  train_data = train,
  train_batch = train_batch,
  test_data = test,
  test_batch = test_batch
)

# Use corrected data
corrected_train <- result$corrected_train
corrected_test <- result$corrected_test
} # }
```
