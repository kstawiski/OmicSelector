# Nogueira Stability Index for Feature Selection

Functions for computing the Nogueira Stability Index, which measures the
consistency of feature selection across cross-validation folds.

## Details

The Nogueira Stability Index is preferred over Jaccard/Kuncheva indices
because it is invariant to the number of selected features. A high
stability (close to 1) indicates that the same features are consistently
selected across different training subsets, suggesting robust signal
rather than overfitting to specific data splits.

Formula: \$\$S = 1 - \frac{V}{V_0}\$\$

Where: - V = (1/p) \* sum(pi_j \* (1 - pi_j)) \[observed variance\] -
V_0 = (k_bar/p) \* (1 - k_bar/p) \[expected variance under random
selection\] - pi_j = selection frequency for feature j - k_bar = average
number of selected features - p = total number of candidate features
