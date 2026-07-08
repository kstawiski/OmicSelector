# Create AutoTuner with MBO

Internal helper to create an AutoTuner with mlr3mbo.

## Usage

``` r
create_autotuner(learner, search_space, n_evals, inner_folds, measure)
```

## Arguments

- learner:

  An mlr3 Learner

- search_space:

  A paradox ParamSet

- n_evals:

  Number of evaluations

- inner_folds:

  Inner CV folds

- measure:

  Measure name or object

## Value

An AutoTuner object
