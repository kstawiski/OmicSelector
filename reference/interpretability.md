# Model Interpretability for OmicSelector 2.0

Provides model-agnostic interpretability via DALEX and SHAP-like
methods. Essential for understanding which features drive predictions
and for generating hypotheses in biomarker discovery.

## Details

\## Why Interpretability Matters

For biomarker discovery, understanding WHY a model makes predictions is
as important as accuracy. Key use cases: - Identify which genes/miRNAs
drive classification - Detect potential confounders or batch effects -
Generate hypotheses for biological validation - Build trust for clinical
adoption

\## Methods Provided

\- \*\*Feature Importance\*\*: Permutation-based importance
(model-agnostic) - \*\*Partial Dependence\*\*: How features affect
predictions marginally - \*\*SHAP Values\*\*: Additive feature
attributions (via iml or approximation) - \*\*Correlation Warnings\*\*:
Flag high-correlation features that may mislead

\## Fold-Aware Design

All interpretability methods respect the cross-validation structure: -
Explanations are computed per-fold using training data only - Aggregated
explanations show stability across folds - Unstable features are flagged
