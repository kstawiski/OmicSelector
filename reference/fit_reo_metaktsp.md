# Fit REO-MetaKTSP meta-analytic pair-order discriminator

Learns a meta-analytic k-Top-Scoring-Pairs discriminator across training
cohorts named by `meta_train[[cohort_col]]`. Each unordered feature pair
receives an avgTSP score: the simple equal-cohort-weight mean of
per-cohort case-control differences in `feature_i < feature_j`. This is
a transfer-estimand method at training because pair selection uses
cross-cohort information, but it remains single-sample at inference
because
[`score_reo_metaktsp`](https://kstawiski.github.io/OmicSelector/reference/score_reo_metaktsp.md)
uses only one specimen's own pair order. When no usable cohort metadata
is available, or only one non-missing cohort is present, fitting
gracefully degenerates to ordinary single-cohort
[`fit_reo_ktsp`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_ktsp.md)
behavior over the training rows. Rows with missing cohort labels are
dropped from meta-analytic selection, cohorts with only cases or only
controls are skipped, and a zero-usable-cohort meta split falls back to
the single-cohort path. At score time, no usable retained pairs returns
the neutral vote fraction `0.5`.

## Usage

``` r
fit_reo_metaktsp(X_train, y_train, meta_train = NULL, hp = list())
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
  non-missing values define the training cohorts used for avgTSP pair
  selection.

- hp:

  List of frozen hyperparameters: `k`, the number of oriented pairs
  requested (default `11L`), and `cohort_col`, a single non-empty
  character string naming the cohort/study column in `meta_train`
  (default `"accession"`). These are resolved once during fitting and
  are not tuned inside `fit_reo_metaktsp()`.

## Value

A plain list of class `reo_metaktsp_model` containing retained oriented
`pairs`, `feature_universe`, retained pair count `k`, `cohort_col`,
number of cohorts used in selection `n_cohorts`, and the resolved `hp`.
`n_cohorts` is the auditable meta-engagement signal: `n_cohorts > 1`
means avgTSP selection ran across multiple cohorts, while
`n_cohorts == 1` means fitting degenerated to single-cohort k-TSP
(NULL/missing/single-level `cohort_col`, or a zero-usable-cohort
fallback). A caller that intends cross-cohort meta-analysis should check
`n_cohorts > 1` so a mis-wired `meta_train` cannot silently demote this
transfer method to pooled within-cohort k-TSP.

## References

Kim S, Lin C-W, Tseng GC. (2016) MetaKTSP: a meta-analytically derived
k-top-scoring-pair classifier for robust prediction of human cancer
subtypes. *Bioinformatics* 32: 1966-1973. PMID: 27153719.

Geman D, d'Avignon C, Naiman DQ, Winslow RL. (2004) Classifying gene
expression profiles from pairwise mRNA comparisons. *Statistical
Applications in Genetics and Molecular Biology* 3: Article19.

Tan AC, Naiman DQ, Xu L, Winslow RL, Geman D. (2005) Simple decision
rules for classifying human cancers from gene expression profiles.
*Bioinformatics* 21: 3896-3904.

## Examples

``` r
if (FALSE) { # \dontrun{
X <- matrix(stats::rgamma(60 * 12, shape = 2), nrow = 60,
            dimnames = list(NULL, paste0("miR-", seq_len(12))))
y <- rep(c(0, 1), times = 30)
meta <- data.frame(accession = rep(paste0("GSE", 1:3), each = 20))
X[y == 0, "miR-1"] <- X[y == 0, "miR-1"] + 5
X[y == 1, "miR-2"] <- X[y == 1, "miR-2"] + 5
model <- fit_reo_metaktsp(X, y, meta_train = meta, hp = list(k = 3L))
score <- score_reo_metaktsp(model, X)
} # }
```
