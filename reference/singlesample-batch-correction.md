# Single-sample batch-correction methods (Module C)

Frozen batch-correction primitives introduced in the single-sample
scoring bank (Module C). All methods follow a fit-then-freeze paradigm:
factor loadings or basis vectors are estimated once on training data,
then applied to test samples without re-estimation, making them suitable
for clinical deployment where a population-matched training cohort may
not be available at prediction time.

Methods provided:

- [`fit_frozen_ruv`](https://kstawiski.github.io/OmicSelector/reference/fit_frozen_ruv.md)
  /
  [`apply_frozen_ruv`](https://kstawiski.github.io/OmicSelector/reference/apply_frozen_ruv.md):
  frozen RUV (Remove Unwanted Variation) using empirical
  negative-control features. Reference: Risso et al. (Nat Biotechnol
  2014).

- [`fit_robust_pca_residual`](https://kstawiski.github.io/OmicSelector/reference/fit_robust_pca_residual.md)
  /
  [`apply_robust_pca_residual`](https://kstawiski.github.io/OmicSelector/reference/apply_robust_pca_residual.md):
  MAD-scaled SVD residual filter. NOT the Candès RPCA (ALM
  sparse-low-rank) decomposition.

## References

Risso D, Ngai J, Speed TP, Dudoit S. (2014) Normalization of RNA-seq
data using factor analysis of control genes or samples. *Nature
Biotechnology* 32(9): 896–902.

Hubert M, Rousseeuw PJ, Vanden Branden K. (2005) ROBPCA: A New Approach
to Robust Principal Component Analysis. *Technometrics* 47(1): 64–79.
