# Score REO-pairratio within-sample log-ratio discriminator

Scores each specimen with the pairwise log-ratio coefficients learned by
[`fit_reo_pairratio`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_pairratio.md).
Each available pair `(a, b)` contributes its frozen coefficient times
`log(X[i, a] + pseudocount) - log(X[i, b] + pseudocount)`. With the
default `strict = FALSE`, pairs missing either feature are omitted. The
frozen intercept and remaining coefficients are unchanged; there is no
rescaling. This is a reduced-support linear predictor, not the complete
fitted model's score. Its discrimination and calibration on a reduced
panel are not established by the complete-panel fit. The returned score
is the linear predictor, not
[`plogis()`](https://rdrr.io/r/stats/Logistic.html) of it.

Empty or constant fitted support, absence of all required features,
all-nondetect rows, and panels without a complete pair raise a
`singlesample_input_error`; see
[`score_specimen`](https://kstawiski.github.io/OmicSelector/reference/score_specimen.md).
Scoring uses no random numbers and no scored-batch statistics.

## Usage

``` r
score_reo_pairratio(model, X, meta = NULL, strict = FALSE)
```

## Arguments

- model:

  A `reo_pairratio_model` object returned by
  [`fit_reo_pairratio`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_pairratio.md).

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

Numeric vector of one finite linear-predictor score per row. Larger
scores are more case-like. The `support` attribute is a data frame with
one row per score: `n_available` counts fitted pairs with both named
columns present, `n_fitted` counts all fitted pairs, and `full_support`
reports whether these counts agree. Counts describe column availability,
not detection or assay quality. Save `attr(score, "support")` before
subsetting or numeric coercion, which can discard attributes.

## References

Aitchison J. (1986) *The Statistical Analysis of Compositional Data*.
Chapman and Hall.

Friedman J, Hastie T, Tibshirani R. (2010) Regularization paths for
generalized linear models via coordinate descent. *Journal of
Statistical Software* 33: 1-22.

## Examples

``` r
# An illustrative frozen two-pair model.
pairs <- data.frame(feature_a = c("a", "c"), feature_b = c("b", "d"))
model <- structure(list(selected_pairs = pairs, coefficients = c(2, -1),
                        intercept = 0.25, pseudocount = 0.5),
                   class = "reo_pairratio_model")
X <- rbind(c(a = 2, b = 1, c = 1, d = 3), c(a = 1, b = 2, c = 4, d = 2))
full <- score_reo_pairratio(model, X, strict = TRUE)
partial <- score_reo_pairratio(model, X[, c("a", "b"), drop = FALSE])
attr(partial, "support")  # 1 available / 2 fitted pairs; full_support = FALSE
#>   n_available n_fitted full_support
#> 1           1        2        FALSE
#> 2           1        2        FALSE
expected <- 0.25 + 2 * (log(X[, "a"] + 0.5) - log(X[, "b"] + 0.5))
stopifnot(isTRUE(all.equal(as.numeric(partial), expected)))
tryCatch(score_reo_pairratio(model, X[, c("a", "b"), drop = FALSE], strict = TRUE),
         singlesample_input_error = function(e) e$code)
#> [1] "MISSING_REQUIRED_FEATURES"
```
