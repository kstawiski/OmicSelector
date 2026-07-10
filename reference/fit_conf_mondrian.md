# Fit Mondrian-conformal class-conditional LRT discriminator

Learns a frozen two-class discriminator from training data using
class-conditional conformal nonconformity scores in per-sample robust
centered log-ratio (rCLR) coordinates. Each class stores raw training
anchors, a centroid in within-sample rCLR space, and a calibration set
of class-specific nonconformity scores. No test data or test-time batch
statistics are used.

## Usage

``` r
fit_conf_mondrian(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be uniquely named feature ids.

- y_train:

  Integer/numeric 0/1 labels aligned to `X_train`; 1 is case/disease and
  0 is control.

- meta_train:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

- hp:

  List of frozen hyperparameters: `nonconformity`, either `"euclid"`
  (default) or `"mahalanobis"` for diagonal class-variance weighting;
  `max_anchors_per_class` positive integer anchor cap per class (default
  `200L`); `min_features` positive integer feature-overlap floor at
  scoring (default `3L`); `eps` positive finite variance floor used by
  the Mahalanobis path (default `1e-6`); and `seed` non-negative integer
  used only for deterministic anchor subsampling while restoring the
  global RNG state (default `1L`).

## Value

A plain list of class `conf_mondrian_model` containing
`feature_universe`, raw `case_anchors` and `control_anchors`,
full-universe `repr` for inspection, and the resolved `hp`.

## Details

The conformal p-value is used here as a discriminative atypicality
contrast, not as an exact coverage-calibrated guarantee: the same
training anchors define the centroids and the calibration nonconformity
scores, so exchangeability is approximate.

## References

Vovk V, Gammerman A, Shafer G. (2005) *Algorithmic Learning in a Random
World*. Springer.

Shafer G, Vovk V. (2008) A tutorial on conformal prediction. *Journal of
Machine Learning Research* 9: 371-421.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(60 * 30, shape = 2), nrow = 60,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 30)
X[y == 1, 1:5] <- X[y == 1, 1:5] + 20
X[y == 0, 6:10] <- X[y == 0, 6:10] + 20
model <- fit_conf_mondrian(X, y)
score <- score_conf_mondrian(model, X)
} # }
```
