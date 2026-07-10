# Create a CatBoost Learner

Create an \`mlr3\` classification learner backed by the \`catboost\` R
package.

## Usage

``` r
make_catboost_learner(
  iterations = 500L,
  learning_rate = 0.03,
  depth = 6L,
  l2_leaf_reg = 5,
  border_count = 64L,
  patience = 40L,
  validation_fraction = 0.15,
  class_weights = "balanced",
  device = c("auto", "cpu", "cuda"),
  seed = 20260331L
)
```

## Arguments

- iterations:

  Maximum boosting iterations. Default: \`500\`.

- learning_rate:

  Learning rate. Default: \`0.03\`.

- depth:

  Tree depth. Default: \`6\`.

- l2_leaf_reg:

  L2 regularization on leaf values. Default: \`5\`.

- border_count:

  Number of numerical split candidates. Default: \`64\`.

- patience:

  Overfitting detector patience for the internal validation split.
  Default: \`40\`.

- validation_fraction:

  Fraction of the training fold reserved for early-stopping. Default:
  \`0.15\`.

- class_weights:

  Class weighting scheme. Use \`"balanced"\` (default), \`"none"\`, or a
  numeric vector with one weight per class.

- device:

  \`"auto"\` (default), \`"cpu"\`, or \`"cuda"\`.

- seed:

  Random seed used inside the backend. Default: \`20260331\`.

## Value

An \`mlr3::LearnerClassif\` with \`predict_type = "prob"\`.

## References

Prokhorenkova, L., et al. (2018). CatBoost: unbiased boosting with
categorical features.
