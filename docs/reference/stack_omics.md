# Create Multi-Omics Stacked Ensemble (Convenience Function)

High-level interface for creating and optionally training a stacked
ensemble.

## Usage

``` r
stack_omics(
  data_list,
  learner_list = NULL,
  y = NULL,
  meta_learner = NULL,
  resampling = NULL,
  task_type = "classif"
)
```

## Arguments

- data_list:

  Named list of data matrices (one per modality)

- learner_list:

  Optional named list of learners (or single learner for all)

- y:

  Optional target vector (if provided, trains immediately)

- meta_learner:

  Learner for stacking

- resampling:

  Resampling strategy

- task_type:

  "classif" or "regr"

## Value

OmicStackedEnsemble object

## Examples

``` r
if (FALSE) { # \dontrun{
# Create and train ensemble
ensemble <- stack_omics(
  data_list = list(mRNA = rna_matrix, miRNA = mirna_matrix),
  learner_list = list(mRNA = "classif.ranger", miRNA = "classif.glmnet"),
  y = y_vector,
  meta_learner = "classif.log_reg"
)

# Predict
preds <- ensemble$predict(list(mRNA = new_rna, miRNA = new_mirna))
} # }
```
