# Fit kernel mean-embedding witness-at-a-point discriminator

Learns a frozen two-class kernel mean-embedding witness function from
training data. Each specimen is represented by a within-sample
compositional transform: it is first closed to proportions \\p = v /
\sum v\\ and then mapped to pseudocount centred log-ratio coordinates
\\\phi(v) = \log(p + pc) - \mathrm{mean}(\log(p + pc))\\. Closing to
proportions first makes \\\phi\\ *exactly* invariant to per-sample
scaling (library size / input amount) for any `pc`, because \\p(c v) =
p(v)\\ for \\c \> 0\\; a raw-abundance pseudocount log-ratio \\\log(v +
pc) - \mathrm{mean}(\log(v + pc))\\ is only approximately
scale-invariant and so would not deliver library-size invariance at the
default pseudocount. The RBF bandwidth is either supplied as `hp$gamma`
or learned by the median heuristic on the stored training anchors:
\$\$gamma = 1 / (2 \* median(\|\|phi(x_i) - phi(x_j)\|\|^2)),\$\$ using
all anchor pairs; if the median squared distance is zero or non-finite,
`gamma = 1`. This bandwidth is panel-frozen: it is learned once over the
full training feature universe and applied unchanged at scoring,
including under partial feature overlap, because re-tuning `gamma` to
the scored overlap would require a scored-batch statistic and break
single-sample deployability. Fitting stores raw case/control anchor
rows, the frozen training feature universe, `pc`, `gamma`, and resolved
hyperparameters. No test data or test-time batch statistic is used.

## Usage

``` r
fit_kme_witness(X_train, y_train, meta_train = NULL, hp = list())
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

  List of frozen hyperparameters: `pc` positive pseudocount on the
  closed proportions (default `1e-3`; small so it bounds the log-ratio
  for absent/near-absent features without flattening the composition),
  `gamma` positive RBF bandwidth (default `NULL`, median heuristic),
  `max_anchors_per_class` positive integer anchor cap per class (default
  `200L`), `min_features` positive integer feature-overlap floor at
  scoring (default `2L`), and `seed` non-negative integer used only for
  deterministic anchor subsampling while restoring the global RNG state
  (default `1L`).

## Value

A plain list of class `kme_witness_model` containing `case_anchors`,
`control_anchors`, `feature_universe`, `gamma`, `pc`, and resolved `hp`.

## References

Gretton A, Borgwardt KM, Rasch MJ, Scholkopf B, Smola A. (2012) A kernel
two-sample test. *Journal of Machine Learning Research* 13: 723-773.

Sutherland DJ, Tung HY, Strathmann H, De S, Ramdas A, Smola A, Gretton
A. (2017) Generative models and model criticism via optimized maximum
mean discrepancy. *International Conference on Learning
Representations*.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(60 * 30, shape = 2), nrow = 60,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 30)
X[y == 1, 1:5] <- X[y == 1, 1:5] + 20
X[y == 0, 6:10] <- X[y == 0, 6:10] + 20
model <- fit_kme_witness(X, y)
score <- score_kme_witness(model, X)
} # }
```
