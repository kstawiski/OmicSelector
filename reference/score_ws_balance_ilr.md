# Score within-cohort ILR balance discriminator

Scores each row of `X` independently with the frozen model learned by
[`fit_ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/fit_ws_balance_ilr.md).
The scorer first restricts the specimen panel to features from the
frozen training universe that are present in `X`. If fewer than
`model$min_sbp_features` canonical SBP features are present, the
specimen/panel receives the documented neutral score `0`. Otherwise,
[`ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/ws_balance_ilr.md)
computes the ILR balance matrix and columns with any non-`NA` value are
summed per row with `na.rm = TRUE`; rows with no valid balance therefore
also receive `0`. The stored orientation sign is then applied so higher
scores are more case-like.

A specimen with fewer than the minimum canonical SBP features has a
degenerate low-information score. This is the single-sample analogue of
the benchmark's panel-coverage ineligibility, surfaced as finite `0`
instead of `NA` so foundation validators can require finite outputs.

## Usage

``` r
score_ws_balance_ilr(model, X, meta = NULL)
```

## Arguments

- model:

  A `ws_balance_ilr_model` object returned by
  [`fit_ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/fit_ws_balance_ilr.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. Scores use only each
row's own values plus the frozen model, with no scored-batch centering,
quantiles, or renormalization.

## References

Egozcue J. J., Pawlowsky-Glahn V. (2005) Groups of parts and their
balances in compositional data analysis. *Mathematical Geology* 37(7):
795-828.

## Examples

``` r
if (FALSE) { # \dontrun{
sbp <- ws_default_sbp()
features <- unique(unlist(lapply(sbp, function(b) {
  c(b$numerator, b$denominator)
}), use.names = FALSE))
features <- setdiff(features, c("ALL_NON_RBC", "ALL_NON_PLT",
                                "TOP_K_BY_ABUNDANCE", "TAIL_BY_ABUNDANCE"))
X <- matrix(stats::rgamma(40 * length(features), shape = 2), nrow = 40,
            dimnames = list(NULL, features))
y <- rep(c(0, 1), each = 20)
model <- fit_ws_balance_ilr(X, y)
score_ws_balance_ilr(model, X[1, , drop = FALSE])
} # }
```
