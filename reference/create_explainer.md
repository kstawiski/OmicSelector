# Create Model Explainer

Creates a DALEX-style explainer for any mlr3 learner. The explainer
wraps the model for uniform interpretability access.

## Usage

``` r
create_explainer(learner, data, y, label = NULL, predict_function = NULL)
```

## Arguments

- learner:

  A trained mlr3 Learner or GraphLearner

- data:

  Data frame or matrix used for explanations

- y:

  Response variable (factor for classification)

- label:

  Optional label for the model

- predict_function:

  Custom predict function (optional)

## Value

An OmicExplainer object

## Examples

``` r
if (FALSE) { # \dontrun{
# Create explainer from trained learner
explainer <- create_explainer(
  learner = trained_learner,
  data = train_data,
  y = train_labels
)

# Compute feature importance
importance <- feature_importance(explainer)
} # }
```
