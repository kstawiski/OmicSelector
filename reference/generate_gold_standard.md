# Generate Gold Standard Synthetic Dataset

Creates a synthetic dataset with known ground truth for validating
feature selection pipelines against data leakage.

## Usage

``` r
generate_gold_standard(
  n_patients = 100,
  samples_per_patient = 3,
  n_features = 1000,
  n_causal = 20,
  n_trap = 30,
  n_batches = 4,
  seed = 42,
  test_ratio = 0.3
)
```

## Arguments

- n_patients:

  Number of patients (default: 100)

- samples_per_patient:

  Number of samples per patient (default: 3)

- n_features:

  Total number of features (default: 1000)

- n_causal:

  Number of true causal features (default: 20)

- n_trap:

  Number of trap/leakage detector features (default: 30)

- n_batches:

  Number of batch effects (default: 4)

- seed:

  Random seed for reproducibility (default: 42)

- test_ratio:

  Ratio of patients for test set (default: 0.3)

## Value

A list of class \`OmicSelectorGoldStandard\` containing:

- data:

  Full dataset (data.frame)

- train_data:

  Training split (data.frame)

- test_data:

  Test split (data.frame)

- train_indices:

  Row indices for training

- test_indices:

  Row indices for testing

- ground_truth:

  Named list of feature indices and names

- train_patient_ids, test_patient_ids:

  Patient IDs per split

- metadata:

  Dataset generation parameters

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate a small dataset for testing
gs <- generate_gold_standard(
  n_patients = 50,
  samples_per_patient = 2,
  n_features = 200,
  seed = 42
)
print(gs)
} # }
```
