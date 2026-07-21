# Fit a random rank forest external competitor

Fits the random-forest algorithm from \`ranktreeEnsemble\` 0.24 entirely
on training-fold data. The audited adapter deliberately fixes
\`datrank=FALSE\`: version 0.24 otherwise converts training columns to
across-sample ranks but predicts from raw within-row pairs, giving
different train and deployment representations. Feature names are
replaced by a frozen collision-free map before the dependency constructs
pair indicators. The preliminary \`randomForestSRC\` screen is
unweighted because version 0.24 exposes no screen-weight hook; the
downstream pair forest uses inverse-frequency case weights, with outcome
level \`case\` (numeric target 1) as the positive class. Screening is
confined to \`rf.cores=1\` because version 0.24 does not forward its
seed to that fit; the caller's prior option is restored.

## Usage

``` r
fit_reo_rankforest(X_train, y_train, meta_train = NULL, hp = list())
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
  \`dimreduce\` (\`TRUE\`), \`class_balance\` (\`TRUE\`), and
  \`max_screen_features\` (128). The last value preserves the backend's
  top-quartile importance screen unless it would be expected to
  materialize more than 128 preliminary features before pair expansion.
  The realized count can differ when preliminary importance values are
  tied. If no pair-order indicator varies in the training data, fitting
  returns an exact intercept-only state instead of invoking a backend
  that cannot split it.

## Value

A frozen \`reo_rankforest_model\` suitable for singleton scoring.

## References

Lu M, Yin R, Chen XS. Ensemble methods of rank-based trees for single
sample classification with gene expression profiles. Journal of
Translational Medicine. 2024;22:140. doi:10.1186/s12967-024-04940-2.
