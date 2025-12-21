# Calibration Metrics for OmicSelector 2.0

Functions for assessing and improving probability calibration of
classification models. Well-calibrated probabilities are essential for
clinical decision-making.

## Details

\## Why Calibration Matters

A model with high AUC can still produce poorly calibrated probabilities.
For clinical use, we need: - When model says 80 - Brier score decomposes
into discrimination + calibration + uncertainty - ECE (Expected
Calibration Error) measures mean deviation from perfect calibration

\## Calibration Repair

Post-hoc calibration can improve probability estimates: - \*\*Platt
scaling\*\*: Fits logistic regression on predictions - \*\*Isotonic
regression\*\*: Non-parametric monotonic recalibration - \*\*Temperature
scaling\*\*: Single parameter softmax recalibration

## References

Niculescu-Mizil, A., & Caruana, R. (2005). Predicting good probabilities
with supervised learning. ICML.
