# Create Bayesian-Optimized AutoTuner for Random Forest

Creates an AutoTuner with ranger and omics-optimized search space.

## Usage

``` r
make_autotuner_ranger(
  task,
  n_evals = NULL,
  inner_folds = 3,
  measure = "classif.auc",
  num_trees = 500L
)
```

## Arguments

- task:

  An mlr3 Task

- n_evals:

  Number of Bayesian optimization iterations

- inner_folds:

  Number of inner CV folds

- measure:

  Performance measure

- num_trees:

  Number of trees (fixed, usually not worth tuning)

## Value

An AutoTuner object
