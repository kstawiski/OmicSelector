# Quick pipeline creation from data

Convenience function to create an OmicPipeline from a data.frame

## Usage

``` r
omic_pipeline(data, target, positive = NULL, ...)
```

## Arguments

- data:

  A data.frame containing features and target

- target:

  Name of the target column

- positive:

  Positive class label (for binary classification)

- ...:

  Additional arguments passed to OmicPipeline\$new()

## Value

An OmicPipeline object
