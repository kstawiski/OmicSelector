# Simple Weighted Averaging Ensemble

Alternative to stacking: weighted average of base model predictions.
Simpler and more interpretable, often performs comparably.

## Public fields

- `modalities`:

  List of OmicModalitySpec objects

- `weights`:

  Named numeric vector of modality weights

- `resampling`:

  Resampling strategy for weight estimation

- `task_type`:

  "classif" or "regr"

## Methods

### Public methods

- [`OmicWeightedEnsemble$new()`](#method-OmicWeightedEnsemble-new)

- [`OmicWeightedEnsemble$train()`](#method-OmicWeightedEnsemble-train)

- [`OmicWeightedEnsemble$predict()`](#method-OmicWeightedEnsemble-predict)

- [`OmicWeightedEnsemble$print()`](#method-OmicWeightedEnsemble-print)

- [`OmicWeightedEnsemble$clone()`](#method-OmicWeightedEnsemble-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new weighted ensemble

#### Usage

    OmicWeightedEnsemble$new(
      modalities,
      weights = NULL,
      resampling = NULL,
      task_type = "classif"
    )

#### Arguments

- `modalities`:

  List of OmicModalitySpec objects

- `weights`:

  Optional named weights. If NULL, estimated from OOF performance.

- `resampling`:

  Resampling strategy

- `task_type`:

  "classif" or "regr"

#### Returns

An OmicWeightedEnsemble object

------------------------------------------------------------------------

### Method `train()`

Train the weighted ensemble

#### Usage

    OmicWeightedEnsemble$train(
      y,
      target_name = "target",
      weight_method = c("inverse_logloss", "auc", "uniform")
    )

#### Arguments

- `y`:

  Target vector

- `target_name`:

  Name for target column

- `weight_method`:

  How to estimate weights: "inverse_logloss", "auc", "uniform"

#### Returns

Self (invisibly)

------------------------------------------------------------------------

### Method [`predict()`](https://rdrr.io/r/stats/predict.html)

Predict using weighted averaging

#### Usage

    OmicWeightedEnsemble$predict(newdata_list)

#### Arguments

- `newdata_list`:

  Named list of data matrices

#### Returns

Data frame with weighted average predictions

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print method

#### Usage

    OmicWeightedEnsemble$print()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    OmicWeightedEnsemble$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
