# Score REO-kTSP within-sample pair-order votes

Scores each specimen using the oriented pairs learned by
[`fit_reo_ktsp`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_ktsp.md).
With the default `strict = FALSE`, only pairs with both features present
contribute. The score is the mean of `feature_a < feature_b` votes over
these available pairs. Missing pairs reduce the denominator: this is a
reduced-support score, not the complete fitted model's score. Its
discrimination and calibration on a reduced panel are not established by
the complete-panel fit.

Empty or constant fitted support, absence of all required features,
all-nondetect rows, and panels with no complete pair raise a
`singlesample_input_error`; see
[`score_specimen`](https://kstawiski.github.io/OmicSelector/reference/score_specimen.md).
Scoring uses only each specimen and the frozen model.

## Usage

``` r
score_reo_ktsp(model, X, meta = NULL, strict = FALSE)
```

## Arguments

- model:

  A `reo_ktsp_model` object returned by
  [`fit_reo_ktsp`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_ktsp.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

- strict:

  Logical scalar. The default `FALSE` preserves scoring over available
  pairs for backward compatibility. Set `TRUE` to reject any missing
  feature from a fitted pair. Unused training features are not required.
  Existing invalid-input and unusable-support checks apply in both
  modes.

## Value

Numeric vector of one finite vote-fraction score per row. With all
fitted-pair features present, numeric values equal
`predict(model$ktsp, X)`. The `support` attribute is a data frame with
one row per score: `n_available` counts fitted pairs with both named
columns present, `n_fitted` counts all fitted pairs, and `full_support`
reports whether these counts agree. Counts describe column availability,
not detection or assay quality. Save `attr(score, "support")` before
subsetting or numeric coercion, which can discard attributes.

## References

Geman D, d'Avignon C, Naiman DQ, Winslow RL. (2004) Classifying gene
expression profiles from pairwise mRNA comparisons. *Statistical
Applications in Genetics and Molecular Biology* 3: Article19.

Tan AC, Naiman DQ, Xu L, Winslow RL, Geman D. (2005) Simple decision
rules for classifying human cancers from gene expression profiles.
*Bioinformatics* 21: 3896-3904.

## Examples

``` r
# An illustrative frozen two-pair model.
pairs <- data.frame(feature_a = c("a", "c"), feature_b = c("b", "d"))
model <- structure(list(ktsp = list(pairs = pairs)), class = "reo_ktsp_model")
X <- rbind(c(a = 2, b = 1, c = 1, d = 3), c(a = 1, b = 2, c = 4, d = 2))
full <- score_reo_ktsp(model, X, strict = TRUE)
partial <- score_reo_ktsp(model, X[, c("a", "b"), drop = FALSE])
attr(partial, "support")  # 1 available / 2 fitted pairs; full_support = FALSE
#>   n_available n_fitted full_support
#> 1           1        2        FALSE
#> 2           1        2        FALSE
stopifnot(identical(as.numeric(partial), c(0, 1)))
tryCatch(score_reo_ktsp(model, X[, c("a", "b"), drop = FALSE], strict = TRUE),
         singlesample_input_error = function(e) e$code)
#> [1] "MISSING_REQUIRED_FEATURES"
```
