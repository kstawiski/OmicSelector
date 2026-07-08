# Train Ratio Image CNN

Fits a pairwise-ratio CNN on log-scale biomarker data and returns
predictions for new samples.

## Usage

``` r
train_ratio_cnn(
  X_train,
  y_train,
  X_test,
  epochs = 30L,
  lr = 0.001,
  batch_size = 256L,
  class_weight = TRUE,
  device = NULL,
  verbose = TRUE,
  seed = 42L,
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

- seed:

  Integer random seed used for reproducible torch initialization.
  Defaults to \`42\`.

- positive:

  Optional positive-class label for non-numeric outcomes.

## Value

A list with elements \`predictions\`, \`fit\`, \`model\`,
\`train_images\`, and \`image_stats\`.

## References

Aitchison J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman and Hall.

Egozcue JJ, Pawlowsky-Glahn V, Mateu-Figueras G, Barcelo-Vidal C.
(2003). Isometric logratio transformations for compositional data
analysis. *Mathematical Geology*, 35(3), 279-300.

Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. (2017). Microbiome
datasets are compositional: and this is not optional. *Frontiers in
Microbiology*, 8, 2224.

Quinn TP, Erb I, Richardson MF, Crowley TM. (2020). Understanding
sequencing data as compositions: an outlook and review.
*Bioinformatics*, 36(16), 4424-4432.

## Examples

``` r
if (FALSE) { # \dontrun{
X_train <- matrix(rnorm(14 * 20), nrow = 20)
y_train <- factor(rep(c("control", "case"), each = 10))
X_test <- matrix(rnorm(14 * 5), nrow = 5)
fit <- train_ratio_cnn(X_train, y_train, X_test, epochs = 5, verbose = FALSE)
fit$predictions
} # }
```
