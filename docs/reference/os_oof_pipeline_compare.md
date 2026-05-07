# Out-of-Fold Pipeline Comparison

Run multiple preprocessing-and-scoring pipelines on the same cohort
under grouped cross-validation, pool out-of-fold predictions per seed,
and report AUCs with clustered bootstrap CIs plus paired DeLong tests vs
a reference pipeline.

## Usage

``` r
os_oof_pipeline_compare(
  y,
  groups,
  pipelines,
  X = NULL,
  batch = NULL,
  plate = NULL,
  n_folds = 5L,
  seeds = c(7L, 42L, 2026L),
  reference = NULL,
  verbose = FALSE
)
```

## Arguments

- y:

  Integer or logical vector of labels (0/1), length `N`.

- groups:

  Vector of grouping IDs, length `N`; samples within the same group
  never appear on opposite sides of a fold split.

- pipelines:

  Named list of pipeline functions. Each function must have signature
  `function(X_train, y_train, X_test, batch_train, batch_test, plate_train, plate_test) -> list(scores_train, scores_test)`.
  The function is responsible for any preprocessing, fold-frozen
  parameter estimation, and learner fitting; it returns numeric
  prediction scores for the train rows (informational) and test rows
  (used).

- X:

  Either a single feature matrix shared across all pipelines, or a named
  list of matrices keyed by pipeline name. When a list is supplied its
  names must match `names(pipelines)` exactly.

- batch:

  Optional character/factor vector, length `N`, passed to each pipeline.
  Use for ComBat-style corrections.

- plate:

  Optional character/factor vector, length `N`, passed to each pipeline.
  Use for plate-median corrections.

- n_folds:

  Integer; number of CV folds (default 5).

- seeds:

  Integer vector of seeds (default `c(7L, 42L, 2026L)`).

- reference:

  Character; name of the pipeline used as comparator in pairwise DeLong
  tests. Set to `NULL` to skip pairwise tests.

- verbose:

  Logical; print fold-level progress.

## Value

A list with:

- `per_seed_auc`:

  data.frame: pipeline, seed, auc, ci_lo, ci_hi.

- `summary`:

  data.frame: pipeline, mean_auc, mean_ci_lo, mean_ci_hi (across seeds).

- `oof_predictions`:

  Named list of length `length(seeds)` — each element is a named list of
  OOF score vectors keyed by pipeline.

- `pairwise_delong`:

  data.frame of paired DeLong results vs `reference` per seed (NULL if
  `reference` is NULL).

- `pairwise_summary`:

  data.frame of meta-combined deltas and Stouffer-combined p-values
  across seeds (NULL if reference NULL).
