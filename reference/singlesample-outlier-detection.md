# Paper 3 outlier detection and conformal claim-gating (Module D)

Outlier detection and conformal anomaly-scoring primitives introduced in
Paper 3 (Module D; Stawiski et al., in preparation). These methods flag
out-of-distribution samples before a biomarker claim is reported,
providing either distribution-free FPR guarantees (conformal approach)
or robust compositional distance metrics (Mahalanobis / isolation
forest).

Methods provided:

- [`fit_compositional_mahalanobis`](https://kstawiski.github.io/OmicSelector/reference/fit_compositional_mahalanobis.md)
  /
  [`apply_compositional_mahalanobis`](https://kstawiski.github.io/OmicSelector/reference/apply_compositional_mahalanobis.md):
  MCD-robust Mahalanobis distance on ILR or rCLR log-ratio coordinates.

- [`fit_conformal_anomaly`](https://kstawiski.github.io/OmicSelector/reference/fit_conformal_anomaly.md)
  /
  [`os_conformal_anomaly`](https://kstawiski.github.io/OmicSelector/reference/os_conformal_anomaly.md):
  conformal p-value via k-NN distance against a held-out calibration
  partition.

- [`fit_isolation_forest_logratio`](https://kstawiski.github.io/OmicSelector/reference/fit_isolation_forest_logratio.md)
  /
  [`apply_isolation_forest_logratio`](https://kstawiski.github.io/OmicSelector/reference/apply_isolation_forest_logratio.md):
  pure-R isolation forest on rCLR-transformed inputs.

## References

Filzmoser P, Hron K, Reimann C. (2009) Principal component analysis for
compositional data with outliers. *Environmetrics* 20: 621-632.

Rousseeuw PJ, Van Driessen K. (1999) A Fast Algorithm for the Minimum
Covariance Determinant Estimator. *Technometrics* 41(3): 212-223.

Vovk V, Gammerman A, Shafer G. (2005) Algorithmic Learning in a Random
World. Springer.

Liu FT, Ting KM, Zhou Z-H. (2008) Isolation Forest. *ICDM*.

Stawiski K. (in preparation) Paper 3 of the OmicSelector programme.
