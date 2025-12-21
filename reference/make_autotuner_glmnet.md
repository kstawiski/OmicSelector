# Create Bayesian-Optimized AutoTuner for glmnet

Creates an AutoTuner with glmnet (elastic net) and search space for
alpha.

## Usage

``` r
make_autotuner_glmnet(
  task,
  n_evals = NULL,
  inner_folds = 3,
  measure = "classif.auc"
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

## Value

An AutoTuner object
