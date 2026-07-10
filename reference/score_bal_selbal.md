# Score selbal-selected single-balance discriminator

Scores each row of `X` independently with the frozen balance learned by
[`fit_bal_selbal`](https://kstawiski.github.io/OmicSelector/reference/fit_bal_selbal.md).
For one specimen, the self-contained balance is computed over the frozen
numerator/denominator groups using only the features that are PRESENT in
`X` and strictly POSITIVE in that specimen: \$\$B = \sqrt{r s / (r + s)}
\\\bigl(\overline{\log v_N} - \overline{\log v_D}\bigr),\$\$ with \\r\\,
\\s\\ the per-specimen counts of present-positive numerator and
denominator features. The returned score is the frozen logistic linear
predictor \\intercept + slope \cdot B\\; larger is more case-like.

If a specimen has fewer than `model$min_group_coverage` present-positive
features in EITHER group (so the balance is undefined/degenerate), it
receives the documented neutral score `0`. If selbal selected an
empty/degenerate balance at fit time (`model$degenerate`), every row
receives `0`. The balance is exactly invariant to per-sample positive
scaling, and scoring uses no random numbers and no scored-batch
statistics.

## Usage

``` r
score_bal_selbal(model, X, meta = NULL)
```

## Arguments

- model:

  A `bal_selbal_model` object returned by
  [`fit_bal_selbal`](https://kstawiski.github.io/OmicSelector/reference/fit_bal_selbal.md).

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

Rivera-Pinto J, Egozcue JJ, Pawlowsky-Glahn V, Paredes R, Noguera-Julian
M, Calle ML. (2018) Balances: a new perspective for microbiome analysis.
*mSystems* 3: e00053-18.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(120 * 15, shape = 30, rate = 2), nrow = 120,
            dimnames = list(NULL, paste0("f", sprintf("%02d", seq_len(15)))))
y <- rep(c(0, 1), each = 60)
X[y == 1, c("f01", "f02", "f03")] <- X[y == 1, c("f01", "f02", "f03")] * 1.8
X[y == 1, c("f13", "f14", "f15")] <- X[y == 1, c("f13", "f14", "f15")] / 1.8
model <- fit_bal_selbal(X, y, hp = list(logit_acc = "Dev"))
score_bal_selbal(model, X[1, , drop = FALSE])
} # }
```
