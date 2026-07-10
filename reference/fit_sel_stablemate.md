# Fit StableMate stable-predictor discriminator (transfer at training, single-sample at inference)

Learns a frozen logistic discriminator over the subset of robust-CLR
(rCLR) features that StableMate identifies as both PREDICTIVE and STABLE
across the training environments (cohorts) named by
`meta_train[[cohort_col]]`. The StableMate idea is that predictors whose
label-association is consistent across cohorts are the batch-robust
ones, while predictors that are informative in a single cohort only are
spurious; the stable set is therefore the transfer-robust feature panel.
This makes `fit_sel_stablemate` a transfer-estimand method at TRAINING
(selection uses the cross-environment structure), but its deployment is
fully SINGLE-SAMPLE: the frozen model is a linear logistic score over
the selected features evaluated on one specimen's own rCLR vector. The
transfer aspect is entirely in the SELECTION.

Pipeline: (1) transform `X_train` to per-sample rCLR over a frozen
candidate universe (the `max_features` highest-variance training
features, frozen, to keep StableMate's bootstrap of `K` selections
tractable); (2) call `StableMate::stablemate` with the cohort factor as
the environment to obtain the predictivity and stability ensembles, then
extract the stable-and-predictive set with the package's own documented
selection rule (`print.stablemate`: predictors significant in the
predictivity ensemble AND in the stability ensemble at `sigthresh`); (3)
fit a frozen `glm(y ~ ., binomial)` over the rCLR of the SELECTED stable
features and store the intercept + coefficients. The entire StableMate
call is wrapped in a global-`.Random.seed` save/restore guard with
`set.seed(hp$seed)` inside it (StableMate's stochastic stepwise
selection `st2e` uses the RNG), so fitting is deterministic and leaves
the caller's RNG state untouched. StableMate is run on the in-process
SEQUENTIAL path (`ncore = NULL`), which executes under our seed in the
main process with no parallel RNG streams; the `ncore = 1L` path is a
foreach/SNOW cluster path with independent RNG streams and a namespace
export bug, so it is deliberately NOT used.

DEGENERATION: if `meta_train` is NULL, the cohort column is missing, or
there is a single usable environment, StableMate cannot assess
cross-environment stability. Fitting then degrades GRACEFULLY to a
predictivity-only selection: a single pooled environment is passed to
StableMate and the predictivity-significant set is used (stability is
undefined with one environment). `n_environments == 1L` records this so
a caller that intends cross-cohort transfer can audit whether a
mis-wired `meta_train` silently demoted this transfer method to pooled
predictivity-only selection. `n_environments` is the auditable
meta-engagement signal: `> 1L` means env-aware stable selection ran.

## Usage

``` r
fit_sel_stablemate(X_train, y_train, meta_train = NULL, hp = list())
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
  non-missing values define the training environments used for env-aware
  stable selection.

- hp:

  List of frozen hyperparameters: `cohort_col` (cohort/environment
  column in `meta_train`, default `"accession"`); `K` (StableMate
  bootstrap selections, default `100L`, integer \\\ge 2\\);
  `max_features` (candidate cap by training variance, default `50L`,
  integer \\\ge 2\\); `min_selected` (minimum present selected features
  in a specimen for a non-neutral score, default `1L`, integer \\\ge
  1\\); `sigthresh` (StableMate significance threshold for the
  predictivity/stability selection, default `0.9`, in `(0, 1)`); and
  `seed` (default `1L`, used to make the StableMate selection
  reproducible while restoring the global RNG state). Resolved once at
  fit and not tuned inside `fit_sel_stablemate()`.

## Value

A plain list of class `sel_stablemate_model` containing the frozen
`selected_features` (the stable-and-predictive rCLR feature ids), the
logistic `intercept` and named `coefficients`, the candidate
`feature_universe` (the frozen rCLR universe), `n_environments` (the
auditable meta-engagement signal: `> 1L` = env-aware stable selection,
`== 1L` = degenerated predictivity-only), `cohort_col`, `min_selected`,
a `degenerate` flag (`TRUE` when no feature was selected, in which case
scoring returns the neutral `0`), and the resolved `hp`.

## References

Deng Y, Liao Y, Wong NKP, Bhuva DD, Lin Y, Cao Y. (2024) StableMate: a
statistical method to select stable predictors in omics data. *bioRxiv*.

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
model <- fit_sel_stablemate(X, y, meta_train = meta)
score <- score_sel_stablemate(model, X)
} # }
```
