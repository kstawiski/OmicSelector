# Extract Selected Features from a Trained Learner

Internal function to extract selected feature names from various learner
types. Handles AutoFSelector, GraphLearner with filter PipeOps, and
other patterns.

## Usage

``` r
.extract_features_from_learner(lrn)
```

## Arguments

- lrn:

  A trained learner object

## Value

Character vector of selected feature names, or NULL if not extractable
