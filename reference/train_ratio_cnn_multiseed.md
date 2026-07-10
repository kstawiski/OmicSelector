# Train Ratio CNN Across Three Default Seeds

Repeats \[train_ratio_cnn()\] across the default reproducibility seeds
\`c(42, 7, 2026)\` and aggregates the predicted probabilities.

## Usage

``` r
train_ratio_cnn_multiseed(
  X_train,
  y_train,
  X_test,
  seeds = .paper1_default_seeds(),
  aggregate = c("mean", "median"),
  epochs = 30L,
  lr = 0.001,
  batch_size = 256L,
  class_weight = TRUE,
  device = NULL,
  verbose = TRUE,
  positive = NULL
)
```

## Arguments

- X_train:

  Numeric matrix with training samples in rows and features in columns.

- y_train:

  Binary outcome vector. Factors, characters, logical values, and
  integer \`0/1\` labels are supported.

- X_test:

  Numeric matrix with test samples in rows and features in columns.

- seeds:

  Integer vector of seeds. Defaults to the reproducibility triplet
  \`c(42, 7, 2026)\`.

- aggregate:

  How to aggregate per-seed predictions; either \`"mean"\` or
  \`"median"\`.

- epochs:

  Number of training epochs. Defaults to \`30\`.

- lr:

  Learning rate. Defaults to \`0.001\`.

- batch_size:

  Mini-batch size. Defaults to \`256\`.

- class_weight:

  Logical; if \`TRUE\`, use inverse-frequency positive-class weighting.
  Defaults to \`TRUE\`.

- device:

  One of \`"cpu"\`, \`"cuda"\`, or \`NULL\`/\`"auto"\` for automatic
  selection.

- verbose:

  Logical; print training progress. Defaults to \`TRUE\`.

- positive:

  Optional positive-class label for non-numeric outcomes.

## Value

A list with aggregated predictions, per-seed predictions, and the fitted
objects.

## References

Quinn TP, Erb I, Richardson MF, Crowley TM. (2020). Understanding
sequencing data as compositions: an outlook and review.
*Bioinformatics*, 36(16), 4424-4432.

## Examples

``` r
if (FALSE) { # \dontrun{
X_train <- matrix(rnorm(14 * 20), nrow = 20)
y_train <- factor(rep(c("control", "case"), each = 10))
X_test <- matrix(rnorm(14 * 5), nrow = 5)
fit <- train_ratio_cnn_multiseed(X_train, y_train, X_test, verbose = FALSE)
fit$predictions
} # }
```
