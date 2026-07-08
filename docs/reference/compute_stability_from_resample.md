# Compute Stability from ResampleResult

Convenience function to compute Nogueira Stability directly from an mlr3
ResampleResult.

## Usage

``` r
compute_stability_from_resample(
  resample_result,
  all_features,
  filter_id = "filter"
)
```

## Arguments

- resample_result:

  An mlr3 ResampleResult object

- all_features:

  Character vector of all candidate features

- filter_id:

  The ID of the filter PipeOp (default: "filter")

## Value

A NogueiraStability object
