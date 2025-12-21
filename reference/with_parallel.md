# With Parallel Scope

Execute code within a parallelization scope, automatically cleaning up.

## Usage

``` r
with_parallel(expr, workers = NULL, plan = "multisession")
```

## Arguments

- expr:

  Expression to evaluate

- workers:

  Number of workers

- plan:

  Parallelization plan

## Value

Result of expression
