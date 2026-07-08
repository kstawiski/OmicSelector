# Bayesian Hyperparameter Optimization for Omics

Provides Bayesian optimization for hyperparameter tuning using mlr3mbo.
Optimized for high-dimensional omics data with sensible default search
spaces for XGBoost, LightGBM, Random Forest, and glmnet.

## Details

Key features: - Omics-optimized search spaces (strong regularization,
shallow trees) - Automatic budget scaling based on sample size - RF
surrogate model (robust for mixed integer/continuous params) - Expected
Improvement (EI) acquisition function

## References

Bischl et al. (2023). mlr3mbo: Bayesian Optimization in mlr3.
