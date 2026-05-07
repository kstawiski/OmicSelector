# Extract Features from All Folds in ResampleResult

Extracts selected features from each fold of a ResampleResult.

## Usage

``` r
extract_features_from_resample(resample_result, filter_id = "filter")
```

## Arguments

- resample_result:

  An mlr3 ResampleResult object

- filter_id:

  The ID of the filter PipeOp (default: "filter")

## Value

A list of character vectors, one per fold
