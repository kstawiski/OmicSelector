# Score within-cohort ILR balance discriminator

Scores each row of `X` independently with the frozen model learned by
[`fit_ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/fit_ws_balance_ilr.md).
The scorer first restricts the specimen panel to features from the
frozen training universe that are present in `X`. Empty or constant
fitted support, absent required features, all-nondetect rows,
insufficient canonical SBP features, and rows with no valid balance
raise a `singlesample_input_error`; see
[`score_specimen`](https://kstawiski.github.io/OmicSelector/reference/score_specimen.md).
Valid balances are summed per row and the frozen orientation sign is
applied so higher scores are more case-like.

## Usage

``` r
score_ws_balance_ilr(model, X, meta = NULL, strict = FALSE)
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

- strict:

  Logical scalar; defaults to `FALSE`, preserving partial balance
  scoring. If `TRUE`, every distinct balance member present at fitting
  must be supplied, or `MISSING_REQUIRED_FEATURES` is raised. Complement
  and abundance-ranked balances use the entire frozen training universe
  as their member set. Other balances count only their named members
  present at fitting. Missing universe metadata fails closed with
  `MISSING_FITTED_UNIVERSE`. Existing input guards apply in both modes.

## Value

Finite numeric vector of length `nrow(X)`. Scores use only each row's
own values plus the frozen model, with no scored-batch centering,
quantiles, or renormalization. The `support` data-frame attribute has
one row per score: `n_available` counts supplied distinct frozen balance
members, `n_fitted` counts all such members, and `full_support` reports
whether these counts agree. This records named-column availability, not
detection, assay quality or calibration. Other training columns can
still affect the default specimen pseudocount. Save
`attr(score, "support")` before subsetting or numeric coercion, which
can discard attributes.

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
