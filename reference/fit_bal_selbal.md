# Fit selbal-selected single-balance discriminator

Learns the two feature groups (numerator vs denominator) of a SINGLE
isometric-log-ratio (ILR) balance from the training data with the selbal
forward-selection algorithm, freezes those groups, and stores a frozen
univariate logistic calibration of the self-contained per-sample
balance. The selbal selection is the only part that touches the training
labels and is therefore leakage-safe: it runs once at fit time, and the
frozen numerator/denominator groups are then evaluated on each specimen
independently. selbal differs from
[`fit_ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/fit_ws_balance_ilr.md)
in exactly one way: instead of a prespecified sequential binary
partition, the two sides of the balance are LEARNED from the training
data.

At inference (see
[`score_bal_selbal`](https://kstawiski.github.io/OmicSelector/reference/score_bal_selbal.md))
the balance over the frozen groups is a self-contained per-sample
log-contrast: \$\$B = \sqrt{r s / (r + s)} \\\bigl(\overline{\log v_N} -
\overline{\log v_D}\bigr),\$\$ where the means are over the features of
each frozen group that are both PRESENT in the scored panel and strictly
POSITIVE in that specimen, and \\r\\, \\s\\ are the per-specimen counts
of those present-positive numerator and denominator features. Because
the balance is a difference of means of logs, it is EXACTLY invariant to
per-sample positive scaling (the per-sample scale cancels), and it uses
only that specimen's own values, so the method is single-sample
deployable and row-equivariant. The calibration recomputes the training
balance through the SAME self-contained helper used at score time (not
selbal's internal `$balance.values`, which is computed over the full
zero-replaced closure) so fit and score share one balance definition and
the stored intercept/slope match exactly what
[`score_bal_selbal`](https://kstawiski.github.io/OmicSelector/reference/score_bal_selbal.md)
computes.

The entire selbal call (whose internal accuracy evaluation may touch the
RNG) is wrapped in a global-`.Random.seed` save/restore guard with
`set.seed(hp$seed)` inside it, so fitting is deterministic and leaves
the caller's RNG state untouched.

## Usage

``` r
fit_bal_selbal(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be uniquely named feature ids.

- y_train:

  Integer/numeric 0/1 labels aligned to `X_train`; 1 is case/disease and
  0 is control.

- meta_train:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

- hp:

  List of frozen hyperparameters: `max_vars` (selbal `maxV`, the cap on
  the total number of features in the balance, default `10L`, integer
  \\\ge 2\\); `logit_acc` (selbal accuracy measure driving forward
  selection, default `"Dev"`, one of `c("AUC", "Dev", "Rsq", "Tjur")`;
  deviance is a continuous objective that keeps adding informative
  features after class separation, recovering a more robust
  multi-feature balance than rank-saturating `"AUC"`, which tends to
  stop early at a tiny balance); `zero_rep` (selbal zero replacement,
  default `"bayes"`, one of `c("bayes", "one")`); `min_group_coverage`
  (minimum present-positive features required in EACH frozen group for a
  non-neutral score, default `1L`, integer \\\ge 1\\); and `seed`
  (default `1L`, used only to make the selbal selection reproducible
  while restoring the global RNG state).

## Value

A plain list of class `bal_selbal_model` containing the frozen
`numerator` and `denominator` feature-name vectors, the logistic
`intercept` and `slope`, an `orientation_sign` fallback (used only when
the logistic slope is degenerate/non-finite), `min_group_coverage`,
`feature_universe`, a `degenerate` flag, and the resolved `hp`. If
selbal returns an empty/degenerate balance (a group with no features, or
overlapping groups), `degenerate` is `TRUE` and
[`score_bal_selbal`](https://kstawiski.github.io/OmicSelector/reference/score_bal_selbal.md)
returns the neutral score `0` for every specimen.

## References

Rivera-Pinto J, Egozcue JJ, Pawlowsky-Glahn V, Paredes R, Noguera-Julian
M, Calle ML. (2018) Balances: a new perspective for microbiome analysis.
*mSystems* 3: e00053-18.

Egozcue JJ, Pawlowsky-Glahn V. (2005) Groups of parts and their balances
in compositional data analysis. *Mathematical Geology* 37(7): 795-828.

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
score <- score_bal_selbal(model, X)
} # }
```
