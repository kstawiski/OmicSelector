# Fit class-conditional Student-t-copula LRT discriminator

Learns a frozen two-class discriminator from training data using a
*Student-t copula* for each class. It is the heavy-tailed generalization
of the Gaussian-copula discriminator
([`fit_lrt_copula`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_copula.md)):
the per-sample representation, the empirical marginals, and the
shrink/ridge correlation are identical, but the copula scores are
Student-t quantiles \\q = t\_\nu^{-1}(u)\\ (shared degrees of freedom
\\\nu = \\ `hp$df`) and the dependence is scored by the Student-t-copula
log-density rather than the Gaussian one. Because \\t\_\nu^{-1} \to
\Phi^{-1}\\ as \\\nu \to \infty\\, this method reduces to `lrt-copula`
in the Gaussian limit.

Each specimen is first mapped to a self-contained, per-sample robust
centred-log-ratio (rCLR) representation \\z\\, where \\z_j = \log v_j -
\mathrm{mean}\_{k:\\v_k\>0}\log v_k\\ on the nonzero parts of the
specimen and \\z_j = 0\\ on its zero parts. Because the centring uses
only that specimen's own nonzero values, \\z(c\\v) = z(v)\\ for any
\\c\>0\\: the representation – and therefore every downstream score – is
exactly invariant to per-sample scaling (library size / input amount).

Each class \\c \in \\\mathrm{case}, \mathrm{control}\\\\ is modelled by
a Student-t copula built from frozen ingredients estimated on *training
only*:

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

2.  a per-class correlation matrix \\R_c\\ of the Student-t scores \\q =
    t\_\nu^{-1}(u)\\ of the training anchors, regularized by shrinkage
    toward the identity (\\R_c \leftarrow (1-s)R_c + sI\\, \\s = \\
    `hp$shrink`) and a ridge that is increased until \\R_c\\ is positive
    definite.

The marginals are matched out by the PIT, so the copula captures only
the class-specific *dependence* structure between features. The single
specimen score (see
[`score_lrt_tcopula`](https://kstawiski.github.io/OmicSelector/reference/score_lrt_tcopula.md))
is the Student-t-copula log-density ratio of case vs control; larger
means more case-like.

Fitting stores the raw case/control anchor rows (so the class copula can
be re-derived consistently over any present-feature subset at scoring),
the frozen training feature universe, the resolved hyperparameters, and
– for inspection and tests – the full-universe class representation. No
test data and no test-time batch statistic are used.

## Usage

``` r
fit_lrt_tcopula(X_train, y_train, meta_train = NULL, hp = list())
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
  \\\delta\\ in \\(0, 0.5)\\ keeping \\t\_\nu^{-1}(u)\\ finite (default
  `1e-4`); `df` Student-t degrees of freedom \\\nu\\, a single finite
  number \\\ge 1\\ (default `5`; \\\nu=1\\ is the Cauchy copula, \\\nu
  \to \infty\\ the Gaussian copula); `max_anchors_per_class`
  positive-integer anchor cap per class (default `200L`); `min_features`
  positive-integer feature-overlap floor at scoring (default `3L`);
  `eps` positive ridge increment / dispersion floor (default `1e-6`);
  and `seed` non-negative integer used only for deterministic anchor
  subsampling while restoring the global RNG state (default `1L`).

## Value

A plain list of class `lrt_tcopula_model` containing `case_anchors`,
`control_anchors`, `feature_universe`, `repr` (the full-universe frozen
class representation), and the resolved `hp`.

## Details

The score is the Student-t-copula log-density ratio \$\$S(z) =
\ell\_{\mathrm{case}}(q\_{\mathrm{case}}) -
\ell\_{\mathrm{control}}(q\_{\mathrm{control}}),\$\$ where \\q_c =
t\_\nu^{-1}(\mathrm{PIT}\_c(z))\\ are the class-\\c\\ Student-t scores
of the specimen over its \\d\\ present features and \$\$\ell_c(q) =
\log\Gamma\\\Big(\tfrac{\nu+d}{2}\Big) +
(d-1)\log\Gamma\\\Big(\tfrac{\nu}{2}\Big) -
d\\\log\Gamma\\\Big(\tfrac{\nu+1}{2}\Big) - \tfrac{1}{2}\log\det R_c -
\tfrac{\nu+d}{2}\\\log\\\Big(1 + \tfrac{q^\top R_c^{-1} q}{\nu}\Big) +
\tfrac{\nu+1}{2}\sum_j \log\\\Big(1 + \tfrac{q_j^2}{\nu}\Big).\$\$ This
is the multivariate-t log-density minus the sum of the univariate-t
marginal log-densities, i.e. the copula (dependence-only) log-density
\\\log f\_{\mathrm{MVt}}(q; R_c, \nu) - \sum_j \log f_t(q_j; \nu)\\; the
\\\tfrac{d}{2}\log(\nu\pi)\\ terms cancel between the joint and the
marginals. The \\\nu/d\\-only \\\log\Gamma\\ constant is identical for
the two classes (shared \\\nu\\, same present \\d\\) and cancels in the
difference; it is kept in each \\\ell_c\\ for clarity. As \\\nu \to
\infty\\, \\-\tfrac{\nu+d}{2}\log(1+\tfrac{q^\top R_c^{-1} q}{\nu}) \to
-\tfrac{1}{2} q^\top R_c^{-1} q\\ and \\\tfrac{\nu+1}{2}\sum_j
\log(1+\tfrac{q_j^2}{\nu}) \to \tfrac{1}{2} q^\top q\\, so \\\ell_c \to
-\tfrac{1}{2} q^\top (R_c^{-1} - I) q - \tfrac{1}{2}\log\det R_c\\, the
Gaussian-copula log-density.

This method is a deliberately *tail-aware* dependence discriminator and
is the intended complement to
[`fit_lrt_copula`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_copula.md)
(the Gaussian-copula limit). At finite \\\nu\\ the Student-t copula with
\\R_c = I\\ is *not* the independence copula: the marginal-correction
term \\\tfrac{\nu+1}{2}\sum_j \log(1+q_j^2/\nu)\\ does not cancel the
joint term exactly (it does only as \\\nu \to \infty\\), so \\\ell_c\\
assigns high density to *jointly tail-extreme* score vectors. The
density is sign-invariant in each coordinate and permutation-invariant
(it depends on the \\q_j\\ only through \\q_j^2\\), but it is *not*
purely radial: only the joint-t term depends on \\q^\top R_c^{-1} q\\
alone, whereas the coordinatewise marginal term is concave in each
\\q_j^2\\, so at a fixed radius it *rewards multi-coordinate tail
extremeness* (mass spread across several features) over a single
dominant coordinate. For example, at \\R_c = I\\, \\\nu = 5\\, the
vectors \\(3, 0)\\ and \\(\sqrt{4.5}, \sqrt{4.5})\\ share the same
quadratic form \\q^\top R_c^{-1} q = 9\\ (Euclidean radius \\\lVert q
\rVert = 3\\) but score \\-0.415\\ and \\+0.347\\ respectively. The LRT
therefore responds to class-specific *joint-tail* behaviour as well as
class-specific correlation, which is the entire point of the
heavy-tailed generalization. A practical consequence is that, unlike the
Gaussian copula (whose exact \\-I\\ term removes the marginal
contribution), this score is *not* designed to capture a pure
per-feature mean/location shift between classes: such a shift is matched
out by the per-class PIT, and a specimen that is marginally extreme
under the opposite class's marginals can be rewarded by that class's
tail. The score is best read as a tail-aware dependence contrast; its
empirical discrimination on any given cohort is established by the
benchmark, not assumed. The orientation is the prespecified
likelihood-ratio contrast (\\\ell\_{\mathrm{case}} -
\ell\_{\mathrm{control}}\\, larger = more case-like) and is never
flipped post hoc.

Degenerate features are handled without any caller-side special-casing.
A feature that is constant across a class's training anchors (e.g. zero
in all of them) collapses to a single-knot marginal, so its PIT maps
directly to that knot's plotting position (Student-t score \\\approx
0\\) with no interpolation; and its zero-variance score column – which
makes [`stats::cor`](https://rdrr.io/r/stats/cor.html) return `NA` – is
forced to the independent (identity) row before the correlation is
shrunk and ridge-regularized to positive definite. The fit and every
score therefore remain finite. `hp$pit_clamp` is validated to the open
interval \\(0, 0.5)\\ so \\t\_\nu^{-1}\\ stays finite, and `hp$df` is
validated to \\\ge 1\\ for tail stability.

## References

Sklar A. (1959) Fonctions de repartition a n dimensions et leurs marges.
*Publications de l'Institut de Statistique de l'Universite de Paris* 8:
229-231.

Demarta S., McNeil A.J. (2005) The t copula and related copulas.
*International Statistical Review* 73(1): 111-129.

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
L[y == 1, 1:8] <- L[y == 1, 1:8] + outer(f, rep(c(1.4, -1.4), 4))  # case block
X <- exp(L)
model <- fit_lrt_tcopula(X, y, hp = list(df = 5))
score <- score_lrt_tcopula(model, X)
} # }
```
