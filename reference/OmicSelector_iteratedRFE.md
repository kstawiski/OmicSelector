# OmicSelector_iteratedRFE

Helper for Marcin Kaszkowiak propritary method. Very longs.. needs
optimizing.

## Usage

``` r
OmicSelector_iteratedRFE(
  trainSet,
  testSet = NULL,
  initFeatures = colnames(trainSet),
  classLab,
  checkNFeatures = 25,
  votingIterations = 1000,
  useCV = F,
  nfolds = 10,
  initRandomState = 42
)
```
