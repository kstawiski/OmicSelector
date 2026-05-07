# Paper 3 matched-null benchmark for within-sample miRNA panel scoring

For an observed scoring method `scoring_fn` applied to a panel of size
`k`, generates `K` random panels of the same size drawn from the cohort
feature pool, jointly stratified on per-feature detection-rate and
log-mean-abundance quartile bins (16 strata total). The empirical
p-value measures whether the observed panel AUC exceeds the matched-null
distribution.

Name note: the simpler (non-stratified) package predecessor is
[`os_panel_null_benchmark`](https://kstawiski.github.io/OmicSelector/reference/os_panel_null_benchmark.md)
(`panel-gates.R`). This function uses the name
`paper3_matched_null_benchmark` to avoid collision.

## Usage

``` r
paper3_matched_null_benchmark(
  X,
  y,
  panel,
  scoring_fn,
  K = 10000L,
  feature_pool = NULL,
  exclude_features = character(0),
  post_hemolysis_corrected = FALSE,
  seed = 42L,
  auc_ci_method = c("hanley_mcneil", "delong", "none")
)
```

## Arguments

- X:

  Matrix (samples \\\times\\ features). Must have column names. Should
  be the full feature pool, not just the observed panel.

- y:

  Binary outcome (0/1 or factor with 2 levels).

- panel:

  Character vector of feature names in the observed panel.

- scoring_fn:

  Function `function(X_subset)` returning a per-sample numeric score.

- K:

  Integer — number of null panels. Default 10,000 (per plan v0.3.3).

- feature_pool:

  Character vector — full feature pool. Defaults to `colnames(X)`.

- exclude_features:

  Character — feature names to exclude from the candidate pool (e.g.,
  hemolysis markers). Default `character(0)`.

- post_hemolysis_corrected:

  Logical. If `TRUE`, prepends the canonical hemolysis-marker list to
  `exclude_features`. Default `FALSE`.

- seed:

  Integer — RNG seed. Default 42.

- auc_ci_method:

  Character; one of `"hanley_mcneil"` (default; closed-form
  Hanley-McNeil 1982 standard error), `"delong"`
  (`pROC::ci.auc(method="delong")`; requires `pROC` \>= 1.16), or
  `"none"` (skip CI). The 95% CI is reported alongside the observed AUC
  as `auc_obs_ci_lo` / `auc_obs_ci_hi`.

## Value

A list with components:

- `auc_obs`:

  Observed panel AUC.

- `auc_obs_ci_lo`, `auc_obs_ci_hi`:

  95% confidence interval for the observed AUC; NA when
  `auc_ci_method = "none"` or computation fails.

- `auc_null`:

  Length-`K` numeric vector of null AUCs (NA for degenerate draws).

- `p_emp`:

  Empirical p-value, or `NA_real_` if not eligible.

- `panel_size`:

  Number of features in the observed panel.

- `K`, `seed`:

  Parameters as supplied.

- `null_panels`:

  List of length `K` of random feature name sets.

- `excluded_features`:

  Final exclusion list applied.

- `eligible`:

  Logical: whether matched-null p-value is valid.

- `fallback_tier`:

  1, 2, or 3 (3 = not eligible).

- `fallback_count`:

  Number of draws that fell back to tier 3.

## References

Stawiski K. (in preparation) Paper 3 of the OmicSelector programme.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
X <- matrix(rlnorm(100 * 50, meanlog = 5, sdlog = 1), nrow = 100,
            dimnames = list(NULL, paste0("mir-", seq_len(50))))
y <- rbinom(100, 1, 0.5)
panel <- paste0("mir-", 1:10)
res <- paper3_matched_null_benchmark(X, y, panel,
         scoring_fn = function(M) rowSums(log(M + 1)),
         K = 200L)
res$p_emp
} # }
```
