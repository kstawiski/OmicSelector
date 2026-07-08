# Convenience Wrapper for Cross-Platform Transfer

Fits a cross-platform adapter on the source data, adapts source and
target matrices, fits or reuses a prediction model, generates
predictions on the target platform, and evaluates performance when
target labels are available.

## Usage

``` r
cross_platform_transfer(
  model,
  source_data,
  source_labels,
  target_data,
  target_labels = NULL,
  strategy = "rank",
  reference_features = NULL,
  reference_summary = c("mean", "median"),
  ...
)
```

## Arguments

- model:

  A training function, an \`mlr3\` learner, or an already fitted model
  object accepted by \[stats::predict()\].

- source_data:

  Numeric source matrix with samples in rows and features in columns.

- source_labels:

  Source labels used for model fitting.

- target_data:

  Numeric target matrix with samples in rows and features in columns.

- target_labels:

  Optional target labels for evaluation.

- strategy:

  Adaptation strategy passed to \[CrossPlatformAdapter\].

- reference_features:

  Optional reference miRNA names for \`strategy = "reference"\`.

- reference_summary:

  Summary statistic for reference normalization.

- ...:

  Additional arguments forwarded to \[CrossPlatformAdapter\].

## Value

A list with the fitted adapter, adapted matrices, fitted model,
predictions, optional evaluation metrics, and domain-shift summaries
before and after adaptation.

## Examples

``` r
set.seed(42)
src <- matrix(rnorm(100 * 5), 100, 5)
tgt <- matrix(rnorm(50 * 5, mean = 1.5), 50, 5)
y_src <- rbinom(100, 1, 0.4)
y_tgt <- rbinom(50, 1, 0.4)

fit_glm <- function(x, y) {
  stats::glm(
    y ~ .,
    data = data.frame(y = as.integer(y), x),
    family = stats::binomial()
  )
}

transfer <- cross_platform_transfer(
  model = fit_glm,
  source_data = src,
  source_labels = y_src,
  target_data = tgt,
  target_labels = y_tgt,
  strategy = "rank"
)

transfer$evaluation$auc
#> [1] 0.5909091
```
