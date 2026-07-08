# Train a CoDaCoRe Classifier

Fits an official CoDaCoRe model when the \`codacore\` package is
available. If the package is unavailable or the TensorFlow backend
errors and \`backend = "auto"\`, the function falls back to a
deterministic sparse-balance logistic model built from ranked CLR
effects.

## Usage

``` r
train_codacore(
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
  backend = c("auto", "codacore", "fallback"),
  max_balances = 5L,
  top_k = 3L,
  cv_folds = 5L,
  fast = TRUE
)
```

## Arguments

- X_train:

  Numeric training matrix with samples in rows and features in columns.

- y_train:

  Binary outcome vector.

- epochs:

  Number of optimization epochs. Defaults to \`30\`.

- lr:

  Learning rate placeholder retained for API compatibility with the
  other Paper 1 learners. The fallback backend does not use it directly.

- batch_size:

  Mini-batch size passed to the official \`codacore\` optimizer.
  Defaults to \`256\`.

- class_weight:

  Logical; if \`TRUE\`, use inverse-frequency weighting in the fallback
  logistic model.

- device:

  Included for API consistency; currently ignored.

- verbose:

  Logical; print backend messages. Defaults to \`TRUE\`.

- seed:

  Integer random seed. Defaults to \`42\`.

- positive:

  Optional positive-class label for non-numeric outcomes.

- input_scale:

  Either \`"log"\` (default) or \`"raw"\`.

- backend:

  One of \`"auto"\`, \`"codacore"\`, or \`"fallback"\`.

- max_balances:

  Maximum number of sparse balances for the fallback model and
  \`maxBaseLearners\` for the official backend.

- top_k:

  Initial sparse balance size used by the fallback model.

- cv_folds:

  Number of internal cross-validation folds passed to the official
  backend.

- fast:

  Logical; passed through to the official \`codacore\` call.

## Value

An object of class \`"omicselector_codacore"\` suitable for
\[predict_codacore()\].

## References

Gordon-Rodriguez E, Susin A, McEwen JD, et al. (2022). CoDaCoRe:
learning sparse log-ratios for high-throughput sequencing data.
*Bioinformatics*, 38(12), 3179-3186.

Aitchison J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman and Hall.

Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. (2017). Microbiome
datasets are compositional: and this is not optional. *Frontiers in
Microbiology*, 8, 2224.

## Examples

``` r
X_train <- matrix(rnorm(14 * 20), nrow = 20)
y_train <- factor(rep(c("control", "case"), each = 10))
fit <- train_codacore(X_train, y_train, backend = "fallback", verbose = FALSE)
fit$backend
#> [1] "fallback"
```
