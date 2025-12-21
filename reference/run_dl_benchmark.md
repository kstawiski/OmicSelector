# Deep Learning Benchmark

Benchmarks deep learning models against traditional ML on an omics task.

## Usage

``` r
run_dl_benchmark(
  task,
  models = c("mlp", "ranger", "glmnet"),
  outer_folds = 5L,
  epochs = 50L,
  seed = NULL
)
```

## Arguments

- task:

  An mlr3 classification task

- models:

  Character vector of models to include: "mlp", "tabtransformer",
  "ranger", "xgboost", "glmnet"

- outer_folds:

  Number of CV folds (default: 5)

- epochs:

  Number of training epochs for DL models (default: 50)

- seed:

  Random seed

## Value

A benchmark result
