# Extract Selected Features from Trained GraphLearner

Extracts the selected feature names from a trained GraphLearner by
inspecting the filter PipeOp's state.

## Usage

``` r
extract_selected_features(learner, filter_id = "filter")
```

## Arguments

- learner:

  A trained GraphLearner object

- filter_id:

  The ID of the filter PipeOp in the graph (default: "filter")

## Value

Character vector of selected feature names, or NULL if extraction fails
