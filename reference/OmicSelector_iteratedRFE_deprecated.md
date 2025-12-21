# Wrapper for deprecated OmicSelector_iteratedRFE

This function is DEPRECATED due to data leakage when testSet is
provided. The test set labels are used during feature selection, which
inflates performance estimates by 10-30

## Usage

``` r
OmicSelector_iteratedRFE_deprecated(trainSet, testSet = NULL, ...)
```

## Arguments

- trainSet:

  Training data

- testSet:

  Test data (CAUSES LEAKAGE if provided)

- ...:

  Other arguments

## Value

Original function result with warning

## Details

\## Why is this deprecated?

The original implementation uses \`testSet\` during RFE iteration: “\`
rfeIter(testX = testSet\[, initFeatures\], testY = testSet\[,
classLab\], ...) “\`

This means the model can "see" test data labels during training, leading
to overly optimistic accuracy estimates.

\## Replacement

Use \`OmicPipeline\$create_graph_learner()\` with nested
cross-validation: “\` pipeline \<- OmicPipeline\$new(data, target =
"outcome") learner \<- pipeline\$create_graph_learner(filter = "anova",
model = "ranger") service \<- BenchmarkService\$new(pipeline,
outer_folds = 5, inner_folds = 3) service\$add_learner(learner) result
\<- service\$run() “\`

## See also

\[OmicPipeline\], \[BenchmarkService\]
