# Create a Frozen ComBat PipeOp for mlr3pipelines

Creates an mlr3pipelines PipeOp that wraps FrozenComBat for use in
leakage-free ML pipelines. The batch correction parameters are learned
during training and frozen for prediction.

## Usage

``` r
create_frozen_combat_pipeop(
  batch_col = "batch",
  parametric = TRUE,
  mean_only = FALSE
)
```

## Arguments

- batch_col:

  Name of the batch column in the task (default: "batch")

- parametric:

  Use parametric empirical Bayes (default: TRUE)

- mean_only:

  Only correct means, not variances (default: FALSE)

## Value

A PipeOp object

## Examples

``` r
if (FALSE) { # \dontrun{
# Create pipeline with frozen ComBat
library(mlr3pipelines)

po_combat <- create_frozen_combat_pipeop(batch_col = "batch")
po_scale <- po("scale")
po_learner <- po("learner", lrn("classif.ranger"))

graph <- po_combat %>>% po_scale %>>% po_learner
} # }
```
