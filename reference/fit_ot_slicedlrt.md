# Fit sliced random-projection Gaussian LRT discriminator

Learns a frozen two-class discriminator from training data by projecting
each specimen's within-sample robust centered log-ratio (rCLR) vector
onto deterministic random one-dimensional slices. Along every slice, the
case and control training anchors are summarized by frozen Gaussian
moments, and scoring uses the average Gaussian log-likelihood ratio.
Larger scores are more case-like.

The rCLR transform is computed separately within each specimen:
`z[pos] = log(v[pos]) - mean(log(v[pos]))` for positive abundances and
`z[!pos] = 0`. This makes the representation invariant to multiplying a
specimen by a positive scalar, without estimating any statistic from the
scoring batch.

## Usage

``` r
fit_ot_slicedlrt(X_train, y_train, meta_train = NULL, hp = list())
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

  List of frozen hyperparameters: `n_slices` positive integer number of
  random projection directions (default `100L`); `max_anchors_per_class`
  positive integer anchor cap per class (default `200L`); `min_features`
  positive integer feature-overlap floor for fitting and scoring
  (default `3L`); `eps` positive finite variance floor (default `1e-6`);
  and `seed` non-negative integer used for deterministic anchor
  subsampling and random directions while restoring the global RNG state
  (default `1L`).

## Value

A plain list of class `ot_slicedlrt_model` containing the frozen
`feature_universe`, raw `case_anchors`, raw `control_anchors`, the
full-universe sliced Gaussian representation `repr`, and resolved `hp`.

## Details

Fitting stores raw training anchors rather than only fitted slice
moments so that scoring can re-derive the same representation over the
feature subset present at inference. The helper `.ot_slicedlrt_repr()`
is used both at fit time over the full training feature universe and at
score time over `intersect(model$feature_universe, colnames(X))`. Its
only inputs are frozen anchors, frozen hyperparameters, and the
present-feature dimension. Consequently, no test-batch centering,
quantile, variance, bandwidth, or other cross-row statistic is estimated
from scored specimens. If the present overlap contains fewer than
`hp$min_features` features,
[`score_ot_slicedlrt`](https://kstawiski.github.io/OmicSelector/reference/score_ot_slicedlrt.md)
returns the documented neutral score `0` for every row.

## References

Bonneel N, Rabin J, Peyre G, Pfister H. (2015) Sliced and Radon
Wasserstein barycenters of measures. *Journal of Mathematical Imaging
and Vision* 51: 22-45.

Neyman J, Pearson ES. (1933) On the problem of the most efficient tests
of statistical hypotheses. *Philosophical Transactions of the Royal
Society of London. Series A* 231: 289-337.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(80 * 30, shape = 2), nrow = 80,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 40)
X[y == 1, 1:6] <- X[y == 1, 1:6] + 25
X[y == 0, 7:12] <- X[y == 0, 7:12] + 25
model <- fit_ot_slicedlrt(X, y)
score <- score_ot_slicedlrt(model, X)
} # }
```
