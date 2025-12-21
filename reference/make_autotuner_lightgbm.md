# Create Bayesian-Optimized AutoTuner for LightGBM

Creates an AutoTuner with LightGBM learner and omics-optimized search
space.

## Usage

``` r
make_autotuner_lightgbm(
  task,
  n_evals = NULL,
  inner_folds = 3,
  measure = "classif.auc",
  nthread = 4L
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

- nthread:

  Number of threads

## Value

An AutoTuner object
