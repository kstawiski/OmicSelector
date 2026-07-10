# Paired DeLong SE for single-sample AUC lift

Computes the empirical DeLong variance of a paired AUC difference
(method score minus baseline score). When \`fold\` is supplied, the
per-fold DeLong variances are aggregated as the variance of the mean
fold-level AUC difference, matching the \`auc_obs_cv =
mean(auc_obs_folds)\` mean-of-folds estimand.

## Usage

``` r
singlesample_paired_auc_diff_se(
  y,
  score_method,
  score_baseline,
  fold = NULL,
  orient = c("median", "auc")
)
```

## Arguments

- y:

  Binary outcome vector, coded 0/1 or a two-level factor.

- score_method:

  Numeric scores from the candidate method.

- score_baseline:

  Numeric scores from the paired baseline method.

- fold:

  Optional fold identifier. If \`NULL\`, all observations are treated as
  one fold.

- orient:

  Orientation rule. \`"median"\` mirrors pROC's within-cohort
  auto-orientation convention; \`"auc"\` reports \`max(AUC, 1 - AUC)\`
  and is used by the transfer-family technology-aware engines.

## Value

A list with \`lift\`, \`se\`, \`variance\`, \`n_folds\`, \`auc_method\`,
\`auc_baseline\`, and a per-fold \`fold_results\` data frame.

## References

DeLong ER, DeLong DM, Clarke-Pearson DL. (1988) Comparing the areas
under two or more correlated receiver operating characteristic curves: a
nonparametric approach. *Biometrics* 44:837-845.
