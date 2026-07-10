# Fit class-conditional vine-copula LRT discriminator

Learns a frozen two-class discriminator from training data using a *vine
copula* (pair-copula construction) for each class. It is the flexible
generalization of the Gaussian-copula
([`fit_lrt_copula`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_copula.md))
and Student-t-copula
([`fit_lrt_tcopula`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_tcopula.md))
discriminators: instead of a single elliptical copula whose dependence
is one correlation matrix, each class is modelled by a vine copula – a
cascade of bivariate pair-copulas, each free to take its own parametric
family (Gaussian, Clayton, Gumbel, Frank) or collapse to independence –
so the class dependence structure can be asymmetric, tail-dependent, and
pair-specific rather than purely elliptical.

Each specimen is first mapped to a self-contained, per-sample robust
centred-log-ratio (rCLR) representation \\z\\, where \\z_j = \log v_j -
\mathrm{mean}\_{k:\\v_k\>0}\log v_k\\ on the nonzero parts of the
specimen and \\z_j = 0\\ on its zero parts. Because the centring uses
only that specimen's own nonzero values, \\z(c\\v) = z(v)\\ for any
\\c\>0\\: the representation – and therefore every downstream score – is
exactly invariant to per-sample scaling (library size / input amount).

Each class \\c \in \\\mathrm{case}, \mathrm{control}\\\\ is modelled by
a vine copula built from two frozen ingredients estimated on *training
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

2.  a per-class vine copula fitted to the pseudo-observations \\u =
    \mathrm{PIT}\_c(z)\\ of the training anchors via
    `rvinecopulib::vinecop`, with a truncation level (`hp$trunc_lvl`)
    that keeps only the strongest dependence trees and a restricted,
    parametric `hp$family_set` (including `"indep"`) so weak pairs
    collapse to independence. The marginals are matched out by the PIT,
    so the vine captures only the class-specific *dependence* structure
    between features.

The single specimen score (see
[`score_lrt_vinecopula`](https://kstawiski.github.io/OmicSelector/reference/score_lrt_vinecopula.md))
is the vine-copula log-density ratio of case vs control; larger means
more case-like.

**Design decisions (flagged for the gate).**

- *Frozen top-\\d\\ vine feature universe.* A vine copula has \\O(d^2)\\
  pair-copulas and is intractable / heavily overfit on the full
  hundreds-of-feature panel with only \\\sim 200\\ anchors. The model is
  therefore restricted to the \\d = \\ `hp$n_vine_features` features
  with the largest *training column mean* (ties broken by column index).
  These highest-abundance miRNAs are also the most platform-portable,
  which suits the single-sample, batch-robust contract. This frozen
  feature set is the model's `feature_universe`; the vine is never
  evaluated outside it.

- *Present-subset re-derivation.* As in the sibling copula scorers, the
  per-class marginals and vine are re-derived from the frozen raw
  anchors over exactly the feature subset present at scoring time
  (`intersect(feature_universe, colnames(X))`) – using only frozen
  training data, never a scored-batch statistic – so partial-overlap
  scores stay consistent and equal the fit-time representation at full
  overlap.

- *Truncated, restricted, parametric vine.* Truncation (`trunc_lvl`,
  default `3L`) and a parametric `family_set` with an explicit
  independence family control overfitting and keep the log-density
  closed-form and deterministic.

Fitting stores the raw case/control anchor rows restricted to the frozen
\\d\\-feature universe (so the class copula can be re-derived
consistently over any present-feature subset at scoring), the frozen
feature universe, the resolved hyperparameters, and – for inspection and
tests – the full-universe class representation. No test data and no
test-time batch statistic are used.

## Usage

``` r
fit_lrt_vinecopula(X_train, y_train, meta_train = NULL, hp = list())
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

  List of frozen hyperparameters: `n_vine_features` positive integer
  \\\ge 2\\ size of the frozen top-mean feature universe the vine is
  built on (default `10L`); `trunc_lvl` positive-integer vine truncation
  level (default `3L`); `family_set` non-empty character vector of
  `rvinecopulib` pair-copula family names (default
  `c("gaussian", "clayton", "gumbel", "frank", "indep")`); `pit_clamp`
  PIT clamp \\\delta\\ in \\(0, 0.5)\\ keeping the pseudo-observations
  in the open unit interval (default `1e-4`); `max_anchors_per_class`
  positive-integer anchor cap per class (default `200L`); `min_features`
  positive-integer feature-overlap floor at scoring (default `3L`);
  `eps` positive numerical floor reserved for interface parity with the
  sibling scorers (default `1e-6`); and `seed` non-negative integer used
  only for deterministic anchor subsampling and as a guard around the
  (deterministic) vine fit while restoring the global RNG state (default
  `1L`).

## Value

A plain list of class `lrt_vinecopula_model` containing `case_anchors`,
`control_anchors` (raw anchor rows restricted to the frozen
\\d\\-feature universe), `feature_universe` (the frozen top-\\d\\
feature ids, in their training-mean order), `repr` (the full-universe
frozen class representation), and the resolved `hp`.

## Details

The score is the vine-copula log-density ratio \$\$S(z) = \log
c\_{\mathrm{case}}(u\_{\mathrm{case}}) - \log
c\_{\mathrm{control}}(u\_{\mathrm{control}}),\$\$ where \\u_c =
\mathrm{PIT}\_c(z)\\ are the class-\\c\\ pseudo-observations of the
specimen over its present features and \\c_c\\ is the class-\\c\\ vine
copula density (`rvinecopulib::dvinecop`). The vine density is by
construction a density on uniform margins, i.e. it is the *copula*
(dependence-only) density: the per-feature marginals are matched out by
the PIT and contribute no term to the score. This matches the sibling
copula scorers, which are likewise copula-only – `lrt-copula` via the
\\-I\\ term that removes the standard-normal marginal contribution,
`lrt-tcopula` by explicitly subtracting the univariate-t marginal
log-densities. A specimen is therefore scored purely by how well its
feature *dependence* matches each class. The orientation is the
prespecified likelihood-ratio contrast (\\\log c\_{\mathrm{case}} - \log
c\_{\mathrm{control}}\\, larger = more case-like) and is never flipped
post hoc.

Degenerate features are handled without any caller-side special-casing.
A feature that is constant across a class's training anchors (e.g. zero
in all of them) collapses to a single-knot marginal, so its PIT maps
directly to that knot's plotting position with no interpolation; its
near-constant pseudo-observation column carries no usable dependence and
is selected to the independence family by `vinecop`. The fit and every
score therefore remain finite. `hp$pit_clamp` is validated to the open
interval \\(0, 0.5)\\ so the pseudo-observations stay strictly inside
the unit interval.

## References

Sklar A. (1959) Fonctions de repartition a n dimensions et leurs marges.
*Publications de l'Institut de Statistique de l'Universite de Paris* 8:
229-231.

Aas K., Czado C., Frigessi A., Bakken H. (2009) Pair-copula
constructions of multiple dependence. *Insurance: Mathematics and
Economics* 44(2): 182-198.

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
model <- fit_lrt_vinecopula(X, y)
score <- score_lrt_vinecopula(model, X)
} # }
```
