# Fit Fisher-Rao geodesic class-conditional LRT discriminator

Learns a frozen two-class discriminator from training data using the
information geometry of compositions. Each specimen is closed to
proportions \\p = v / \sum v\\ (optionally additively smoothed toward
the uniform composition by the pseudocount, \\(p + pc) / \sum (p +
pc)\\) and mapped to the square-root embedding \\u = \sqrt{p}\\. Because
\\\sum_i p_i = 1\\ we have \\\\u\\ = 1\\, so every composition lands on
the positive orthant of the unit sphere \\S^{D-1}\\, the manifold on
which the Fisher-Rao (multinomial information) metric is the round
metric. The Fisher-Rao geodesic distance between two compositions is the
great-circle distance of their embeddings, \$\$d(p, q) = 2 \\
\arccos\\\Big(\sum_i \sqrt{p_i \\ q_i}\Big),\$\$ the term inside being
the Bhattacharyya coefficient \\= \langle u_p, u_q \rangle\\. Closing to
proportions first – and applying the pseudocount to the already-closed
proportions – makes the embedding (and hence every distance and the
final score) *exactly* invariant to per-sample scaling (library size /
input amount) for any `pc`, because \\p(c v) = p(v)\\ for \\c \> 0\\ and
the smoothing is a fixed function of the scale-invariant \\p\\.

Each class \\c \in \\case, control\\\\ is modelled as a
Riemannian-Gaussian on the sphere: a class mean \\\mu_c\\ and a scalar
Fisher-Rao dispersion \\\sigma_c^2\\. The mean is the extrinsic
(projected) Frechet mean \\\mu_c = \bar u_c / \\\bar u_c\\\\ with \\\bar
u_c\\ the average of the class embeddings; compositional data are
confined to one orthant of the sphere (any two embeddings have
non-negative inner product, so the round-sphere angle between them,
\\\arccos\langle u, w \rangle\\, is \\\le \pi/2\\ – half the injectivity
radius \\\pi\\; the factor-2 Fisher-Rao distance \\d = 2\arccos\\ is
then \\\le \pi\\), where the extrinsic mean is a consistent,
deterministic, closed-form proxy for the iterative Karcher mean. The
dispersion is the mean squared geodesic distance of the class embeddings
to \\\mu_c\\ (the Frechet variance), floored at `hp$eps` so a
near-degenerate class cannot produce an infinite weight.

Fitting stores the raw case/control anchor rows (so the class
representation can be re-derived consistently over any present-feature
subset at scoring), the frozen training feature universe, the resolved
hyperparameters, and – for inspection and tests – the full-universe
class means and dispersions. No test data and no test-time batch
statistic are used.

## Usage

``` r
fit_ig_fisherrao(X_train, y_train, meta_train = NULL, hp = list())
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

  List of frozen hyperparameters: `pc` non-negative pseudocount applied
  as additive smoothing of the closed proportions toward the uniform
  composition (default `0`, the pure Fisher-Rao geometry; a positive
  value moves boundary/zero compositions off the simplex boundary while,
  by acting on the already-closed proportions, preserving exact
  scale-invariance for any `pc`), `shared_dispersion` logical (default
  `FALSE`; `TRUE` pools a single within-class dispersion, reducing the
  score to a Fisher-Rao nearest-centroid contrast), `eps` positive
  dispersion floor (default `1e-6`), `max_anchors_per_class` positive
  integer anchor cap per class (default `200L`), `min_features` positive
  integer feature-overlap floor at scoring (default `2L`), and `seed`
  non-negative integer used only for deterministic anchor subsampling
  while restoring the global RNG state (default `1L`).

## Value

A plain list of class `ig_fisherrao_model` containing `case_anchors`,
`control_anchors`, `feature_universe`, `mu_case`, `mu_control`,
`sigma2_case`, `sigma2_control`, `sigma2_pooled` (all over the full
training universe), and the resolved `hp`.

## Details

The single-specimen score (see
[`score_ig_fisherrao`](https://kstawiski.github.io/OmicSelector/reference/score_ig_fisherrao.md))
is the class-conditional Riemannian-Gaussian log-likelihood ratio
\$\$S(x) = \frac{d(x, \mu\_{ctrl})^2}{2 \sigma\_{ctrl}^2} - \frac{d(x,
\mu\_{case})^2}{2 \sigma\_{case}^2},\$\$ i.e. \\\log p\_{case}(x) - \log
p\_{ctrl}(x)\\ up to an additive constant (the class log-normalizers and
the log-dispersion-ratio term) that does not depend on \\x\\ and so
leaves ranking, AUC, and any monotone threshold calibration unchanged.
Because each dispersion is estimated in the same squared-distance units
as the numerator, the arbitrary Fisher-Rao distance constant (the factor
2 above) cancels in every \\d^2/\sigma^2\\ ratio whenever the dispersion
is the estimated Frechet variance; the fixed `hp$eps` floor is a
numerical regularizer that does not rescale with the convention, so for
a degenerate (floor-bound) class that convention-independence is only
approximate. Larger \\S\\ means geometrically closer
(dispersion-weighted) to the case class.

## References

Rao CR. (1945) Information and accuracy attainable in the estimation of
statistical parameters. *Bulletin of the Calcutta Mathematical Society*
37: 81-91.

Amari S, Nagaoka H. (2000) *Methods of Information Geometry*. American
Mathematical Society / Oxford University Press.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(60 * 30, shape = 2), nrow = 60,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 30)
X[y == 1, 1:5] <- X[y == 1, 1:5] + 20
X[y == 0, 6:10] <- X[y == 0, 6:10] + 20
model <- fit_ig_fisherrao(X, y)
score <- score_ig_fisherrao(model, X)
} # }
```
