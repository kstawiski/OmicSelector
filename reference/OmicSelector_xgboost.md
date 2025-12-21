# OmicSelector_xgboost

Train xgboost model with Bayesian optimalization. Code based on
http://www.mysmu.edu/faculty/jwwang/post/hyperparameters-tuning-for-xgboost-using-bayesian-optimization/

## Usage

``` r
OmicSelector_xgboost(
  features = "all",
  train = OmicSelector_load_datamix(use_smote_not_rose = T)[[1]],
  test = OmicSelector_load_datamix(use_smote_not_rose = T)[[2]],
  valid = OmicSelector_load_datamix(use_smote_not_rose = T)[[2]],
  eta = c(0, 1),
  gamma = c(0, 100),
  max_depth = c(2L, 10L),
  min_child_weight = c(1, 25),
  subsample = c(0.25, 1),
  nfold = c(3L, 10L),
  initPoints = 8,
  iters.n = 10
)
```

## Arguments

- features:

  Vector of features to be used. If "all", all features starting with
  'hsa' will be used.

- train:

  Training dataset with column Class ('Case' vs. 'Control') and features
  starting with 'hsa'.

- test:

  Testing dataset with column Class ('Case' vs. 'Control') and features
  starting with 'hsa'.

- valid:

  Testing dataset with column Class ('Case' vs. 'Control') and features
  starting with 'hsa'.

- eta:

  Bonderies of 'eta' parameter in XGBoost training, must be a vector of
  2.

- gamma:

  Bonderies of 'gamma' parameter in XGBoost training, must be a vector
  of 2.

- max_depth:

  Bonderies of 'max_depth' parameter in XGBoost training, must be a
  vector of 2.

- min_child_weight:

  Bonderies of 'min_child_weight' parameter in XGBoost training, must be
  a vector of 2.

- subsample:

  Bonderies of 'subsample' parameter in XGBoost training, must be a
  vector of 2.

- nfold:

  Bonderies of 'nfold' parameter in XGBoost training, must be a vector
  of 2.

## Value

Xgboost model
