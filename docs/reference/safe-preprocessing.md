# Safe Preprocessing Utilities for Cross-Validation

Functions that enforce leakage-free preprocessing within
cross-validation folds. These utilities prevent the most common data
leakage bugs in biomarker validation: global median imputation, global
scaling, and test-fold self-imputation.
