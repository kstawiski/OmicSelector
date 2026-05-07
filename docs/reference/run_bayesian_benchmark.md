# Run Bayesian Optimization Benchmark

Convenience function to run Bayesian-optimized learners in nested CV and
compare performance.

## Usage

``` r
run_bayesian_benchmark(
  task,
  learners = c("xgboost", "ranger", "glmnet"),
  outer_folds = 5,
  inner_folds = 3,
  n_evals = NULL,
  measure = "classif.auc",
  seed = NULL
)
```

## Arguments

- task:

  An mlr3 Task

- learners:

  Character vector of learner names: "xgboost", "lightgbm", "ranger",
  "glmnet"

- outer_folds:

  Number of outer CV folds (default: 5)

- inner_folds:

  Number of inner CV folds (default: 3)

- n_evals:

  Number of MBO iterations per learner (default: auto)

- measure:

  Performance measure (default: classif.auc)

- seed:

  Random seed

## Value

A benchmark result with tuned learners

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_bayesian_benchmark(
  task,
  learners = c("xgboost", "ranger", "glmnet"),
  outer_folds = 5,
  n_evals = 30
)
result$aggregate(msr("classif.auc"))
} # }
```
