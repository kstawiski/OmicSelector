# Get Selected Features Per Fold

Extract the raw list of selected features for each fold from a nested CV
result. Useful for debugging feature selection and examining
fold-to-fold variability.

## Usage

``` r
get_selected_features_per_fold(nested_result, learner_id = NULL)
```

## Arguments

- nested_result:

  A NestedCVResult object from BenchmarkService\$run()

- learner_id:

  Optional learner ID to filter (default: all learners)

## Value

A data.table with columns: learner_id, fold, n_features, features (list
column)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- benchmark_service$run()
features_df <- get_selected_features_per_fold(result)
print(features_df)

# View features for a specific learner/fold
features_df[learner_id == "my_learner" & fold == 1]$features[[1]]
} # }
```
