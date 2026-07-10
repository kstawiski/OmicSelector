# Fit Invariant Causal Prediction discriminator (transfer at training, single-sample at inference)

Learns a frozen logistic discriminator over the subset of robust-CLR
(rCLR) features that are both PREDICTIVE of the label and INVARIANT
(cohort-stable) across the training environments (cohorts) named by
`meta_train[[cohort_col]]`, following the Invariant Causal Prediction
(ICP) principle of Peters, Buhlmann & Meinshausen (2016): a predictor
whose conditional relationship with the label is INVARIANT across
environments is a plausible (batch-robust, transferable) causal signal,
whereas a predictor whose effect VARIES by cohort signals confounding /
non-causal, batch-driven association and is dropped. This makes
`fit_icp` a transfer-estimand method at TRAINING (selection uses the
cross-environment structure), but its deployment is fully SINGLE-SAMPLE:
the frozen model is a linear logistic score over the selected features
evaluated on one specimen's own rCLR vector. The transfer aspect is
entirely in the SELECTION.

icp is the structural twin of
[`fit_sel_stablemate`](https://kstawiski.github.io/OmicSelector/reference/fit_sel_stablemate.md)
with ONE change: the feature-SELECTION step. Where sel-stablemate calls
StableMate's predictivity+stability bootstrap ensembles, icp selects
features by an explicit, tractable per-feature ICP invariance screen:

1.  PREDICTIVITY screen: a pooled marginal logistic
    `glm(y ~ rclr_j, binomial)` over ALL training rows; feature \\j\\ is
    PREDICTIVE if the Wald p-value of its slope \\\< \alpha\_{pred}\\
    (default `0.05`).

2.  INVARIANCE screen (the ICP core): a per-cohort marginal logistic
    `glm(y ~ rclr_j, binomial)` WITHIN each VALID cohort (\\\ge 1\\ case
    AND \\\ge 1\\ control AND \\\ge\\ `min_env_n` rows), giving
    per-cohort slopes \\b_j^{(e)}\\ with standard errors \\SE_e\\. The
    cross-cohort HOMOGENEITY of \\\\b_j^{(e)}\\\\ is tested with a
    Cochran's fixed-effect heterogeneity statistic \\Q = \sum_e w_e
    (b_e - \bar{b})^2\\, \\w_e = 1/SE_e^2\\, \\\bar{b} = \sum_e w_e b_e
    / \sum_e w_e\\, with \\Q \sim \chi^2\_{E-1}\\ under homogeneity.
    Feature \\j\\ is INVARIANT if the Q p-value \\\> \alpha\_{inv}\\
    (default `0.10`) — i.e. we FAIL to reject homogeneity, so the \\y
    \mid x_j\\ relationship is stable across cohorts.

3.  SELECT feature \\j\\ if PREDICTIVE AND INVARIANT.

A frozen `glm(y ~ ., binomial)` is then fit over the rCLR of the
SELECTED features, pooled over all training rows, and the intercept +
per-feature coefficients are stored. The per-cohort glm sweep can touch
the RNG only through deterministic IRLS (it does not), but the whole fit
is nevertheless wrapped in a global-`.Random.seed` save/restore guard
with `set.seed(hp$seed)` inside it, so fitting is deterministic and
leaves the caller's RNG state byte-unchanged — mirroring the
sel-stablemate contract.

DEGENERATION (icp ALWAYS produces a usable, non-empty head):

- If `meta_train` is NULL, the cohort column is missing, or there are
  \\\< 2\\ VALID cohorts, the invariance screen cannot be evaluated and
  selection degrades GRACEFULLY to the PREDICTIVITY-ONLY screen (the
  pooled marginal logistic alone). `n_environments == 1L` records this
  so a caller that intends cross-cohort transfer can audit whether a
  mis-wired `meta_train` silently demoted this transfer method to pooled
  predictivity-only selection. `n_environments > 1L` means the ICP
  invariance screen ran.

- If no feature passes both screens, the selection RELAXES (documented
  in `model$selection_mode`): take the most-homogeneous (highest Q
  p-value) among the predictive features, up to `top_k` (default
  `min(5, n_predictive)`); if still empty (no predictive feature at
  all), take the single most-predictive feature. This guarantees a
  usable non-empty head while keeping the relaxation auditable.

Unlike sel-stablemate (a stability method that honestly returns the
EMPTY set when nothing is stable), icp's brief mandates an
always-non-empty head; the relaxation taken is therefore recorded in
`model$selection_mode` rather than left empty.

## Usage

``` r
fit_icp(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be uniquely named feature ids.

- y_train:

  Integer/numeric 0/1 labels aligned to `X_train`; 1 is case/disease and
  0 is control.

- meta_train:

  Optional per-sample metadata. When it contains `cohort_col`,
  non-missing values define the training environments used for the ICP
  invariance screen.

- hp:

  List of frozen hyperparameters: `cohort_col` (cohort/environment
  column in `meta_train`, default `"accession"`); `max_features`
  (candidate cap by training variance, default `50L`, integer \\\ge
  2\\); `alpha_pred` (predictivity Wald p-value threshold, default
  `0.05`, in `(0, 1)`); `alpha_inv` (invariance Cochran's-Q p-value
  threshold; a feature is invariant when its Q p-value EXCEEDS this,
  default `0.10`, in `(0, 1)`); `min_env_n` (minimum rows for a cohort
  to be a valid environment, default `8L`, integer \\\ge 4\\); `top_k`
  (relaxation cap on most-homogeneous-among-predictive features when the
  strict set is empty, default `5L`, integer \\\ge 1\\); `min_selected`
  (minimum present selected features in a specimen for a non-neutral
  score, default `1L`, integer \\\ge 1\\); and `seed` (default `1L`,
  used to restore the global RNG state). Resolved once at fit and not
  tuned inside `fit_icp()`.

## Value

A plain list of class `icp_model` containing the frozen
`selected_features` (the predictive-and-invariant rCLR feature ids), the
logistic `intercept` and named `coefficients`, the candidate
`feature_universe` (the frozen rCLR universe), `n_environments` (the
auditable meta-engagement signal: `> 1L` = ICP invariance screen ran,
`== 1L` = degenerated predictivity-only), `selection_mode` (one of
`"invariant_predictive"`, `"predictivity_only"`, `"relaxed_top_k"`,
`"relaxed_single_most_predictive"`, recording any relaxation taken),
`cohort_col`, `min_selected`, a `degenerate` flag (`TRUE` only if NO
feature could be selected even after relaxation, in which case scoring
returns the neutral `0`), and the resolved `hp`.

## References

Peters J, Buhlmann P, Meinshausen N. (2016) Causal inference by using
invariant prediction: identification and confidence intervals. *Journal
of the Royal Statistical Society: Series B* 78(5): 947-1012.

Cochran WG. (1954) The combination of estimates from different
experiments. *Biometrics* 10(1): 101-129.

Aitchison J. (1986) *The Statistical Analysis of Compositional Data*.
Chapman and Hall.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 180L; p <- 12L
env <- rep(c("e1", "e2", "e3"), length.out = n)
X <- matrix(stats::rgamma(n * p, shape = 5), n, p,
            dimnames = list(NULL, paste0("f", seq_len(p))))
lin <- log(X[, 1]) + log(X[, 2]) + log(X[, 3])
y <- stats::rbinom(n, 1, plogis(scale(lin)))
meta <- data.frame(accession = env, stringsAsFactors = FALSE)
model <- fit_icp(X, y, meta_train = meta)
score <- score_icp(model, X)
} # }
```
