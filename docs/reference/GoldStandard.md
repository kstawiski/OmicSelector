# Gold Standard Synthetic Dataset for Leakage Detection

Generates and evaluates synthetic datasets designed to detect data
leakage in feature selection and cross-validation pipelines.

## Details

The Gold Standard dataset includes: - \*\*Causal features\*\*: Known
ground truth signal correlated with outcome - \*\*Trap features\*\*:
Leakage detectors that appear informative only when test data improperly
leaks into training - \*\*Batch features\*\*: Nuisance features
confounded with batch effects - \*\*Pathway blocks\*\*: Correlated
feature groups mimicking biological pathways

A correctly implemented cross-validation pipeline should: 1. Select
mostly CAUSAL features 2. Select ZERO TRAP features (selecting any
indicates leakage)
