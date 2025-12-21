# Create Robust Fallback Learner

Creates the best available fallback learner when mlr3torch is not
installed. Follows the fallback chain: xgboost -\> ranger -\> glmnet

## Usage

``` r
.create_robust_fallback()
```

## Value

An mlr3 Learner object
