# OmicStackedEnsemble R6 Class

Orchestrates multi-omics stacking with leakage-free OOF prediction
generation. Trains base learners per modality and combines using a
meta-learner.

Based on Super Learner methodology (van der Laan et al., 2007).

## References

van der Laan, M.J., Polley, E.C., Hubbard, A.E. (2007). Super Learner.
Statistical Applications in Genetics and Molecular Biology, 6(1).

## Public fields

- `modalities`:

  List of OmicModalitySpec objects

- `meta_learner`:

  Learner for combining base predictions

- `resampling`:

  Resampling strategy for OOF predictions

- `prediction_type`:

  "prob" for probabilities, "response" for class labels

- `task_type`:

  "classif" or "regr"

## Methods

### Public methods

- [`OmicStackedEnsemble$new()`](#method-OmicStackedEnsemble-new)

- [`OmicStackedEnsemble$train()`](#method-OmicStackedEnsemble-train)

- [`OmicStackedEnsemble$predict()`](#method-OmicStackedEnsemble-predict)

- [`OmicStackedEnsemble$is_fitted()`](#method-OmicStackedEnsemble-is_fitted)

- [`OmicStackedEnsemble$get_oof_stack()`](#method-OmicStackedEnsemble-get_oof_stack)

- [`OmicStackedEnsemble$get_modality_weights()`](#method-OmicStackedEnsemble-get_modality_weights)

- [`OmicStackedEnsemble$print()`](#method-OmicStackedEnsemble-print)

- [`OmicStackedEnsemble$clone()`](#method-OmicStackedEnsemble-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new stacked ensemble

#### Usage

    OmicStackedEnsemble$new(
      modalities,
      meta_learner = NULL,
      resampling = NULL,
      prediction_type = "prob",
      task_type = "classif"
    )

#### Arguments

- `modalities`:

  List of OmicModalitySpec objects

- `meta_learner`:

  Learner for stacking (default: logistic regression)

- `resampling`:

  Resampling strategy (default: 5-fold CV)

- `prediction_type`:

  "prob" or "response"

- `task_type`:

  "classif" or "regr"

#### Returns

An OmicStackedEnsemble object

------------------------------------------------------------------------

### Method `train()`

Train the stacked ensemble

#### Usage

    OmicStackedEnsemble$train(y, target_name = "target")

#### Arguments

- `y`:

  Target vector (named with sample IDs)

- `target_name`:

  Name for target column

#### Returns

Self (invisibly), for method chaining

------------------------------------------------------------------------

### Method [`predict()`](https://rdrr.io/r/stats/predict.html)

Predict using the stacked ensemble

#### Usage

    OmicStackedEnsemble$predict(newdata_list)

#### Arguments

- `newdata_list`:

  Named list of data matrices, one per modality. Names must match
  modality IDs. Missing modalities can be NULL.

#### Returns

Prediction object

------------------------------------------------------------------------

### Method `is_fitted()`

Check if ensemble is fitted

#### Usage

    OmicStackedEnsemble$is_fitted()

#### Returns

Logical

------------------------------------------------------------------------

### Method `get_oof_stack()`

Get OOF stacking features (for auditing)

#### Usage

    OmicStackedEnsemble$get_oof_stack()

#### Returns

Data frame with OOF predictions used for meta-learner training

------------------------------------------------------------------------

### Method `get_modality_weights()`

Get modality contribution (meta-learner coefficients if linear)

#### Usage

    OmicStackedEnsemble$get_modality_weights()

#### Returns

Named vector or NULL if meta-learner is not linear

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print method

#### Usage

    OmicStackedEnsemble$print()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    OmicStackedEnsemble$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
