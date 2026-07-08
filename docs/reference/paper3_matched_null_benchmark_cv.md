# 5-fold nested-CV matched-null benchmark

Nested cross-validation variant of
[`paper3_matched_null_benchmark`](https://kstawiski.github.io/OmicSelector/reference/paper3_matched_null_benchmark.md).
Outer folds are stratified by outcome; when `group_id` is supplied
(e.g., matched-set patient IDs from a provenance audit), grouped folds
guarantee no patient group is split across train/test. For each outer
fold a panel of size `panel_size` is selected on the TRAINING fold by
univariate AUC, and the matched-null strata are computed on the training
fold; the held-out test fold is then used to compute both the observed
AUC and `K` matched-null AUCs.

This is the engine that produced the manuscript's per-cell numbers in
Supplementary Table S1 / Figure 3 / Table 2.

## Usage

``` r
paper3_matched_null_benchmark_cv(
  X,
  y,
  panel_size = 20L,
  scoring_fn = NULL,
  K = 10000L,
  outer_k = 5L,
  group_id = NULL,
  feature_pool = NULL,
  exclude_features = character(0),
  post_hemolysis_corrected = FALSE,
  seed = 42L,
  min_valid_folds = 4L,
  min_pos_per_fold = 5L,
  min_neg_per_fold = 5L,
  null_parallel_cores = 1L,
  auc_ci_method = c("hanley_mcneil", "delong", "none")
)
```

## Arguments

- X:

  Numeric matrix (samples \\\times\\ features) with column names.

- y:

  Binary outcome (0/1 or factor with two levels).

- panel_size:

  Integer. Panel size selected per outer fold on training data. Default
  20.

- scoring_fn:

  Function `function(X_panel)` returning per-sample scores; if `NULL`
  the default `rowSums(log(X + 0.5))` is used. The closure must NEVER
  look at held-out data.

- K:

  Integer. Null replicates per fold. Default 10000.

- outer_k:

  Integer. Outer CV folds. Default 5.

- group_id:

  Vector of length `nrow(X)`, or `NULL`. When provided, grouped K-fold
  is used and groups are never split.

- feature_pool:

  Character. Full candidate pool. `NULL` uses `colnames(X)`.

- exclude_features:

  Character. Features excluded from the null candidate pool.

- post_hemolysis_corrected:

  Logical. If `TRUE` prepends the canonical hemolysis-marker list to
  `exclude_features`.

- seed:

  Integer in \[0, 21474\]. Fold-assignment and per-draw null seeds are
  derived deterministically from this value.

- min_valid_folds:

  Integer. Default 4 (manuscript / plan v0.6 sec 4.1).

- min_pos_per_fold:

  Integer. Default 5.

- min_neg_per_fold:

  Integer. Default 5.

- null_parallel_cores:

  Integer. Number of forked workers per fold. Default 1 (serial).

- auc_ci_method:

  Character; one of `"hanley_mcneil"` (default), `"delong"` (uses pooled
  fold-level scores via `pROC::ci.auc(method = "delong")`), or `"none"`.
  Returned as `auc_obs_ci_lo` / `auc_obs_ci_hi` on the cohort-level
  pooled-fold prediction.

## Value

Named list including `auc_obs_cv`, `auc_obs_ci_lo`, `auc_obs_ci_hi`,
`p_emp_cv`, `n_valid_folds`, `eligible`, `fold_panels`, and a per-fold
audit. See the manuscript Methods §"Statistical evaluation".

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(99)
X <- matrix(rlnorm(150 * 120, 5, 1), nrow = 150,
            dimnames = list(NULL, paste0("f-", seq_len(120))))
y <- rbinom(150, 1, 0.5)
X[y == 1, paste0("f-", 1:10)] <- X[y == 1, paste0("f-", 1:10)] * 3
res <- paper3_matched_null_benchmark_cv(X, y, panel_size = 10L,
                                         K = 100L, outer_k = 5L)
res$p_emp_cv
} # }
```
