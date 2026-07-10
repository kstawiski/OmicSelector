# Fit class-conditional Gaussian-copula LRT discriminator

Learns a frozen two-class discriminator from training data using a
*Gaussian copula* for each class. Each specimen is first mapped to a
self-contained, per-sample robust centred-log-ratio (rCLR)
representation \\z\\, where \\z_j = \log v_j -
\mathrm{mean}\_{k:\\v_k\>0}\log v_k\\ on the nonzero parts of the
specimen and \\z_j = 0\\ on its zero parts. Because the centring uses
only that specimen's own nonzero values, \\z(c\\v) = z(v)\\ for any
\\c\>0\\: the representation – and therefore every downstream score – is
exactly invariant to per-sample scaling (library size / input amount).

Each class \\c \in \\\mathrm{case}, \mathrm{control}\\\\ is modelled by
a Gaussian copula. The copula is built from two frozen ingredients
estimated on *training only*:

1.  a per-feature empirical marginal used as a probability-integral
    transform (PIT). The marginal is the *linearly interpolated*
    plotting-position ECDF through the training order statistics with
    Weibull positions \\u\_{(i)} = i / (n_c + 1)\\ (constant
    extrapolation past the training min/max), clamped to \\\[\delta,
    1-\delta\]\\, \\\delta = \\ `hp$pit_clamp`. It equals the step
    plotting-position ECDF \\\\\\z^{\mathrm{train}} \le x\\ / (n_c +
    1)\\ at the data points but is continuous in \\x\\, so the score is
    exactly (to numerical tolerance) invariant to per-sample scaling
    instead of jumping when the irreducible floating-point noise of
    rescaling pushes \\z\\ across a breakpoint; and

2.  a per-class correlation matrix \\R_c\\ of the Gaussian scores \\q =
    \Phi^{-1}(u)\\ of the training anchors, regularized by shrinkage
    toward the identity (\\R_c \leftarrow (1-s)R_c + sI\\, \\s = \\
    `hp$shrink`) and a ridge that is increased until \\R_c\\ is positive
    definite.

The marginals are matched out by the PIT, so the copula captures only
the class-specific *dependence* structure between features. The single
specimen score (see
[`score_lrt_copula`](https://kstawiski.github.io/OmicSelector/reference/score_lrt_copula.md))
is the Gaussian-copula log-density ratio of case vs control; larger
means more case-like.

Fitting stores the raw case/control anchor rows (so the class copula can
be re-derived consistently over any present-feature subset at scoring),
the frozen training feature universe, the resolved hyperparameters, and
– for inspection and tests – the full-universe class representation. No
test data and no test-time batch statistic are used.

## Usage

``` r
fit_lrt_copula(X_train, y_train, meta_train = NULL, hp = list())
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

  List of frozen hyperparameters: `shrink` correlation shrinkage toward
  the identity in \\\[0, 1)\\ (default `0.05`); `pit_clamp` PIT clamp
  \\\delta\\ in \\(0, 0.5)\\ keeping \\\Phi^{-1}(u)\\ finite (default
  `1e-4`); `max_anchors_per_class` positive-integer anchor cap per class
  (default `200L`); `min_features` positive-integer feature-overlap
  floor at scoring (default `3L`); `eps` positive ridge increment /
  dispersion floor (default `1e-6`); and `seed` non-negative integer
  used only for deterministic anchor subsampling while restoring the
  global RNG state (default `1L`).

## Value

A plain list of class `lrt_copula_model` containing `case_anchors`,
`control_anchors`, `feature_universe`, `repr` (the full-universe frozen
class representation), and the resolved `hp`.

## Details

The score is the Gaussian-copula log-density ratio \$\$S(z) =
\ell\_{\mathrm{case}}(z) - \ell\_{\mathrm{control}}(z),\qquad \ell_c(z)
= -\tfrac{1}{2} q_c^\top (R_c^{-1} - I) q_c - \tfrac{1}{2}\log\det
R_c,\$\$ where \\q_c = \Phi^{-1}(\mathrm{PIT}\_c(z))\\ are the
class-\\c\\ Gaussian scores of the specimen. The \\-I\\ term removes the
standard-normal marginal contribution, so \\\ell_c\\ is the copula
(dependence-only) log-density; with \\R_c = I\\ the copula log-density
is identically \\0\\, i.e. a specimen is scored purely by how well its
feature *dependence* matches each class.

Degenerate features are handled without any caller-side special-casing.
A feature that is constant across a class's training anchors (e.g. zero
in all of them) collapses to a single-knot marginal, so its PIT maps
directly to that knot's plotting position (Gaussian score \\\approx 0\\)
with no interpolation; and its zero-variance Gaussian-score column –
which makes [`stats::cor`](https://rdrr.io/r/stats/cor.html) return `NA`
– is forced to the independent (identity) row before the correlation is
shrunk and ridge-regularized to positive definite. The fit and every
score therefore remain finite. `hp$pit_clamp` is validated to the open
interval \\(0, 0.5)\\ so \\\Phi^{-1}\\ stays finite.

## References

Sklar A. (1959) Fonctions de repartition a n dimensions et leurs marges.
*Publications de l'Institut de Statistique de l'Universite de Paris* 8:
229-231.

Joe H. (2014) *Dependence Modeling with Copulas*. Chapman and Hall/CRC.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 120; p <- 30
L <- matrix(stats::rnorm(n * p, 4, 0.6), nrow = n,
            dimnames = list(NULL, paste0("miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
f <- stats::rnorm(sum(y == 1))
L[y == 1, 1:8] <- L[y == 1, 1:8] + outer(f, rep(1.4, 8))  # case-only block
X <- exp(L)
model <- fit_lrt_copula(X, y)
score <- score_lrt_copula(model, X)
} # }
```
