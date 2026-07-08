# Train a CLR + MLP Classifier

Fits a dense CoDA baseline: a centered log-ratio transform followed by a
multilayer perceptron. The default backend uses torch when it is
available, and otherwise falls back to \`nnet::nnet()\` or logistic
regression for environments without torch.

## Usage

``` r
train_clr_mlp(
  X_train,
  y_train,
  epochs = 30L,
  lr = 0.001,
  batch_size = 256L,
  class_weight = TRUE,
  device = NULL,
  verbose = TRUE,
  seed = 42L,
  positive = NULL,
  input_scale = c("log", "raw"),
  backend = c("auto", "torch", "nnet", "glm"),
  hidden_units = c(152L, 136L),
  dropout = 0.3
)
```

## Arguments

- X_train:

  Numeric training matrix with samples in rows and features in columns.

- y_train:

  Binary outcome vector.

- epochs:

  Number of training epochs. Defaults to \`30\`.

- lr:

  Learning rate. Defaults to \`0.001\`.

- batch_size:

  Mini-batch size. Defaults to \`256\`.

- class_weight:

  Logical; if \`TRUE\`, use inverse-frequency weighting.

- device:

  One of \`"cpu"\`, \`"cuda"\`, or \`NULL\`/\`"auto"\`.

- verbose:

  Logical; print training progress. Defaults to \`TRUE\`.

- seed:

  Integer random seed. Defaults to \`42\`.

- positive:

  Optional positive-class label for non-numeric outcomes.

- input_scale:

  Either \`"log"\` (default) or \`"raw"\`.

- backend:

  One of \`"auto"\`, \`"torch"\`, \`"nnet"\`, or \`"glm"\`.

- hidden_units:

  Integer vector with the two hidden-layer sizes used by the torch
  backend. Defaults to \`c(152, 136)\`.

- dropout:

  Dropout rate used by the torch backend. Defaults to \`0.3\`.

## Value

An object of class \`"omicselector_clr_mlp"\` suitable for
\[predict_clr_mlp()\].

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
X_train <- matrix(rnorm(14 * 20), nrow = 20)
y_train <- factor(rep(c("control", "case"), each = 10))
fit <- train_clr_mlp(X_train, y_train, backend = "glm", verbose = FALSE)
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
fit$backend
#> [1] "glm"
```
