# Create Bayesian-Optimized AutoTuner for XGBoost

Creates an AutoTuner with XGBoost learner and omics-optimized search
space.

## Usage

``` r
make_autotuner_xgboost(
  task,
  n_evals = NULL,
  inner_folds = 3,
  measure = "classif.auc",
  early_stopping = TRUE,
  nthread = 4L
)
```

## Arguments

- task:

  An mlr3 Task (used to set data-aware defaults)

- n_evals:

  Number of Bayesian optimization iterations (default: auto based on n)

- inner_folds:

  Number of inner CV folds for tuning (default: 3)

- measure:

  Performance measure (default: classif.auc)

- early_stopping:

  Use early stopping with nrounds (default: TRUE)

- nthread:

  Number of threads for XGBoost (default: 4)

## Value

An AutoTuner object ready for nested CV

## Examples

``` r
if (FALSE) { # \dontrun{
at <- make_autotuner_xgboost(task, n_evals = 30)
rr <- resample(task, at, rsmp("cv", folds = 5))
} # }
```
