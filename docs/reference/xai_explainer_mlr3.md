# Create DALEX Explainer from mlr3 Learner

Converts a trained mlr3 learner into a DALEX explainer for
interpretability.

## Usage

``` r
xai_explainer_mlr3(learner, task, data = NULL, label = NULL, features = NULL)
```

## Arguments

- learner:

  A trained mlr3 Learner (must have predict_type = "prob")

- task:

  The mlr3 task used for training

- data:

  Optional: explicit data for explainer (if NULL, uses task data)

- label:

  Label for the explainer (default: learner ID)

- features:

  Optional: subset of features to include (speeds up explanations)

## Value

A DALEX explainer object

## Examples

``` r
if (FALSE) { # \dontrun{
learner <- lrn("classif.ranger", predict_type = "prob")
learner$train(task)
explainer <- xai_explainer_mlr3(learner, task)
} # }
```
