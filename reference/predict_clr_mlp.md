# Predict with a CLR + MLP Classifier

Generates positive-class probabilities for new samples from a fit
created by \[train_clr_mlp()\].

## Usage

``` r
predict_clr_mlp(fit, X_new)
```

## Arguments

- fit:

  A fitted \`"omicselector_clr_mlp"\` object.

- X_new:

  Numeric matrix of new samples.

## Value

A numeric vector of positive-class probabilities.

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
predict_clr_mlp(fit, X_train[1:3, , drop = FALSE])
#> [1] 1 1 1
```
