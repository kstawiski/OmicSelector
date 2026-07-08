# Get Consensus Features from Best Signature

Extract the consensus feature set from the selected best signature(s).
Aggregates features selected across outer folds.

## Usage

``` r
get_consensus_features(nested_result, best_signature, min_frequency = 0.5)
```

## Arguments

- nested_result:

  A NestedCVResult object

- best_signature:

  Result from select_best_signature() or learner_id string

- min_frequency:

  Minimum proportion of folds where feature must be selected (default:
  0.5)

## Value

A data.table with columns: feature, frequency, selected (TRUE if \>=
min_frequency)

## Examples

``` r
if (FALSE) { # \dontrun{
best <- select_best_signature(result)
consensus <- get_consensus_features(result, best, min_frequency = 0.8)
final_features <- consensus[selected == TRUE]$feature
} # }
```
