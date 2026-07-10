# Fit the Group-DRO kit x biofluid robust discriminator

Canonical-contract fit wrapper for the Group-DRO kit x biofluid scorer
primitive. The training matrix, labels, and (optional) kit x biofluid
metadata are forwarded to
[`fit_group_dro_scorer`](https://kstawiski.github.io/OmicSelector/reference/fit_group_dro_scorer.md),
which performs all train-only CLR panel selection, optional nested-CV
hyperparameter tuning, and exponentiated-gradient Group-DRO logistic
optimisation. When `meta_train` is `NULL` or lacks the configured kit /
biofluid columns the fit degenerates gracefully to a single constant
group (ERM); the realised number of groups is recorded as an auditable
engagement signal.

## Usage

``` r
fit_dro_group(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with named feature columns.

- y_train:

  Numeric / integer 0/1 labels (1 = case), length `nrow(X_train)`, with
  both classes present.

- meta_train:

  Optional per-sample metadata carrying the kit and biofluid columns
  named by `hp$kit_col` / `hp$biofluid_col`. When `NULL` or missing
  those columns, the fit runs as a single-group (ERM) Group-DRO model.

- hp:

  Optional list of hyperparameters (strict allow-list). Fields, with
  primitive-matching defaults: `kit_col` (`"kit_family"`),
  `biofluid_col` (`"biofluid"`), `panel_size` (`20L`), `panel_features`
  (`NULL`; train-only CLR panel selected when unset), `eta_grid`
  (`c(0.01, 0.05, 0.10)`), `l2_grid` (`c(0.001, 0.01, 0.10)`),
  `inner_folds` (`3L`), `learning_rate` (`0.05`), `epochs` (`200L`),
  `tune` (`TRUE`), `seed` (`42L`).

## Value

Object of class `dro_group_model`: a list with the frozen primitive
`fit` (a `group_dro_scorer` object), `n_groups` (the number of distinct
kit x biofluid groups engaged), `groups_engaged` (`n_groups > 1`),
`train_groups`, the resolved `kit_col` / `biofluid_col`, and `hp`.

## References

Sagawa S, Koh PW, Hashimoto TB, Liang P. (2020) Distributionally Robust
Neural Networks for Group Shifts. *ICLR*. arXiv:1911.08731.

## See also

[`score_dro_group`](https://kstawiski.github.io/OmicSelector/reference/score_dro_group.md),
[`fit_group_dro_scorer`](https://kstawiski.github.io/OmicSelector/reference/fit_group_dro_scorer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 80; p <- 25; k <- 6
L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
            dimnames = list(NULL, paste0("hsa-miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
X <- exp(L)
meta <- data.frame(kit_family = rep(c("A", "B"), length.out = n),
                   biofluid = rep(c("plasma", "serum"), each = n / 2))
model <- fit_dro_group(X, y, meta)
score_dro_group(model, X[1, , drop = FALSE])
} # }
```
