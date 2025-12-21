# Stability-Based Ensemble Selection

Builds an ensemble of models using stability-based feature selection.
The approach: 1. Run N bootstrap resamples with feature selection 2.
Compute feature selection frequencies 3. Create nested overlapping
subsets based on stability tiers 4. Train K base models on different
stability tiers 5. Combine predictions via stability-weighted voting

This approach addresses the "reproducibility crisis" in biomarker
discovery by prioritizing stable, reproducible feature sets over
single-run selections.

## References

Meinshausen & Buhlmann (2010). Stability Selection. JRSS-B. Nogueira et
al. (2018). On the Stability of Feature Selection Algorithms.

## Public fields

- `n_bootstrap`:

  Number of bootstrap resamples

- `filter`:

  Feature selection filter method

- `n_features`:

  Number of features to select per bootstrap

- `base_learner`:

  Base learner for ensemble members

- `threshold_mode`:

  "adaptive" (target set size) or "fixed" (percentage)

- `target_set_size`:

  Target number of stable features (for adaptive mode)

- `fixed_thresholds`:

  Thresholds for fixed mode (creates nested subsets)

- `aggregation`:

  "weighted" (stability-weighted) or "average" (simple)

- `verbose`:

  Print progress messages

- `frequencies`:

  Feature selection frequencies from last fit

- `stable_features`:

  List of stable feature sets by tier

- `tier_weights`:

  Weights for each tier in ensemble

- `fitted_models`:

  List of fitted models per tier

## Methods

### Public methods

- [`StabilityEnsemble$new()`](#method-StabilityEnsemble-new)

- [`StabilityEnsemble$compute_frequencies()`](#method-StabilityEnsemble-compute_frequencies)

- [`StabilityEnsemble$identify_stable_subsets()`](#method-StabilityEnsemble-identify_stable_subsets)

- [`StabilityEnsemble$create_tier_learner()`](#method-StabilityEnsemble-create_tier_learner)

- [`StabilityEnsemble$fit()`](#method-StabilityEnsemble-fit)

- [`StabilityEnsemble$predict()`](#method-StabilityEnsemble-predict)

- [`StabilityEnsemble$get_feature_importance()`](#method-StabilityEnsemble-get_feature_importance)

- [`StabilityEnsemble$get_summary()`](#method-StabilityEnsemble-get_summary)

- [`StabilityEnsemble$print()`](#method-StabilityEnsemble-print)

- [`StabilityEnsemble$clone()`](#method-StabilityEnsemble-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new StabilityEnsemble

#### Usage

    StabilityEnsemble$new(
      n_bootstrap = 100,
      filter = "anova",
      n_features = 50,
      base_learner = "glmnet",
      threshold_mode = "adaptive",
      target_set_size = c(10, 50),
      fixed_thresholds = c(0.9, 0.7, 0.5),
      aggregation = "weighted",
      verbose = TRUE
    )

#### Arguments

- `n_bootstrap`:

  Number of bootstrap resamples (default: 100)

- `filter`:

  Feature selection filter: "anova", "mrmr", "importance"

- `n_features`:

  Features to select per bootstrap (default: 50)

- `base_learner`:

  Learner for ensemble: "glmnet", "ranger", "xgboost"

- `threshold_mode`:

  "adaptive" or "fixed"

- `target_set_size`:

  Target stable set size range for adaptive mode

- `fixed_thresholds`:

  Thresholds for fixed mode

- `aggregation`:

  "weighted" or "average"

- `verbose`:

  Print progress

#### Returns

A new StabilityEnsemble object

------------------------------------------------------------------------

### Method `compute_frequencies()`

Compute feature selection frequencies from bootstrap resamples

#### Usage

    StabilityEnsemble$compute_frequencies(task, seed = NULL)

#### Arguments

- `task`:

  An mlr3 Task

- `seed`:

  Random seed for reproducibility

#### Returns

Named numeric vector of selection frequencies

------------------------------------------------------------------------

### Method `identify_stable_subsets()`

Identify stable feature subsets based on frequency thresholds

#### Usage

    StabilityEnsemble$identify_stable_subsets()

#### Returns

List of feature vectors for each stability tier

------------------------------------------------------------------------

### Method `create_tier_learner()`

Create mlr3 learner for a specific feature subset

#### Usage

    StabilityEnsemble$create_tier_learner(features)

#### Arguments

- `features`:

  Character vector of feature names

#### Returns

An mlr3 Learner

------------------------------------------------------------------------

### Method `fit()`

Fit the stability ensemble on a task

#### Usage

    StabilityEnsemble$fit(task, seed = NULL)

#### Arguments

- `task`:

  An mlr3 Task

- `seed`:

  Random seed

#### Returns

The fitted StabilityEnsemble (invisibly)

------------------------------------------------------------------------

### Method [`predict()`](https://rdrr.io/r/stats/predict.html)

Predict on new data using the ensemble

#### Usage

    StabilityEnsemble$predict(task)

#### Arguments

- `task`:

  An mlr3 Task (or data.frame with same features)

#### Returns

A data.table with predicted probabilities

------------------------------------------------------------------------

### Method `get_feature_importance()`

Get feature importance based on selection frequencies

#### Usage

    StabilityEnsemble$get_feature_importance(top_n = NULL)

#### Arguments

- `top_n`:

  Number of top features to return (default: all)

#### Returns

Named vector of feature frequencies, sorted descending

------------------------------------------------------------------------

### Method `get_summary()`

Get summary of the ensemble configuration and results

#### Usage

    StabilityEnsemble$get_summary()

#### Returns

A list with ensemble summary

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print summary

#### Usage

    StabilityEnsemble$print()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    StabilityEnsemble$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
