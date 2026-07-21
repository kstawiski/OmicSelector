# Fit boosted rank trees external competitor

Fits the boosting algorithm from \`ranktreeEnsemble\` 0.24 using the
same fold-local, collision-free feature mapping and deployable
\`datrank=FALSE\` contract as \[fit_reo_rankforest()\]. The backend's
exposed \`weights\` argument is not forwarded by version 0.24 and is
therefore not represented as active in this adapter.

## Usage

``` r
fit_reo_rankboost(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric non-negative matrix of training samples by features.

- y_train:

  Numeric/integer 0/1 training labels; 1 is the case class.

- meta_train:

  Optional training metadata, checked for row alignment and otherwise
  ignored.

- hp:

  Optional fixed hyperparameters: \`ntree\` (500), \`seed\` (1),
  \`dimreduce\` (\`TRUE\`), \`nodedepth\` (3), \`nodesize\` (5),
  \`shrinkage\` (0.05), requested \`bag_fraction\` (0.5), and
  \`max_screen_features\` (128). The effective bag fraction applies a
  conservative training-size-only one-observation margin above the
  pinned \`gbm\` node-size boundary; both the requested and effective
  values are retained in the fitted state. A training set with no
  variable pair ordering returns an exact intercept-only state.

## Value

A frozen \`reo_rankboost_model\` suitable for singleton scoring.
