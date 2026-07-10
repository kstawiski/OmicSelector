# Predict with a Weighted CLR + MLP Classifier

Generates positive-class probabilities for new samples from a fit
created by \[train_weighted_clr_mlp()\].

## Usage

``` r
predict_weighted_clr_mlp(fit, X_new)
```

## Arguments

- fit:

  A fitted \`"omicselector_weighted_clr_mlp"\` object.

- X_new:

  Numeric matrix of new samples.

## Value

Numeric vector of positive-class probabilities.
