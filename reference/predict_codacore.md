# Predict with a CoDaCoRe Classifier

Generates positive-class probabilities for new samples from a fit
created by \[train_codacore()\].

## Usage

``` r
predict_codacore(fit, X_new)
```

## Arguments

- fit:

  A fitted \`"omicselector_codacore"\` object.

- X_new:

  Numeric matrix of new samples.

## Value

A numeric vector of positive-class probabilities.

## References

Gordon-Rodriguez E, Susin A, McEwen JD, et al. (2022). CoDaCoRe:
learning sparse log-ratios for high-throughput sequencing data.
*Bioinformatics*, 38(12), 3179-3186.

## Examples

``` r
X_train <- matrix(rnorm(14 * 20), nrow = 20)
y_train <- factor(rep(c("control", "case"), each = 10))
fit <- train_codacore(X_train, y_train, backend = "fallback", verbose = FALSE)
predict_codacore(fit, X_train[1:3, , drop = FALSE])
#> [1] 0.4790714 0.7788656 0.2962273
```
