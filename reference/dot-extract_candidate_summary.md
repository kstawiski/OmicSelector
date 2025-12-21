# Extract Candidate Summary from Nested CV Results

Internal function to extract per-candidate statistics across outer
folds.

## Usage

``` r
.extract_candidate_summary(nested_result, metric, metric_higher_better)
```

## Arguments

- nested_result:

  NestedCVResult object

- metric:

  Performance metric column name

- metric_higher_better:

  Whether higher metric is better

## Value

data.table with candidate-level summary
