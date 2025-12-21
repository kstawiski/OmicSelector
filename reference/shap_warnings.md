# Correlation-Aware SHAP Interpretation

Functions to compute SHAP values with warnings about correlated
features. High correlations between features can make SHAP values
misleading because: 1. Importance can be arbitrarily split between
correlated features 2. SHAP values may not reflect true causal
relationships 3. Removing one correlated feature could dramatically
change others' SHAP values

## References

Hooker et al. (2021). Unrestricted Permutation Forces Extrapolation.
Lundberg & Lee (2017). SHAP: A Unified Approach to Interpreting Model
Predictions.
