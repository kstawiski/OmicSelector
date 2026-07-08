# Safe Evaluation Utilities for Biomarker Validation

Functions that enforce correct ROC computation and non-inferiority
testing. These utilities prevent the most common evaluation bugs:
hardcoded ROC direction (inverting AUC), wrong non-inferiority
references, and unstable percentile estimates from insufficient random
draws.
