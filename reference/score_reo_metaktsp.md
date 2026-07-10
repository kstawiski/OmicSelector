# Score REO-MetaKTSP within-sample pair-order votes

Scores each row of `X` independently with the frozen oriented pairs
learned by
[`fit_reo_metaktsp`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_metaktsp.md).
The model is trained meta-analytically across cohorts via
`meta_train[[cohort_col]]` using avgTSP with equal cohort weight, but
score-time inference is single-sample: each retained pair votes 1 only
when `feature_a < feature_b` within that specimen. The graceful
single-cohort degeneration, missing-cohort-row drop, and
single-class-cohort skip are fit-time behaviors only. If no retained
pair has both features present in `X`, scoring returns the neutral vote
fraction `0.5` for every row.

## Usage

``` r
score_reo_metaktsp(model, X, meta = NULL)
```

## Arguments

- model:

  A `reo_metaktsp_model` object returned by
  [`fit_reo_metaktsp`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_metaktsp.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Numeric vector of one finite vote-fraction score per row of `X`. Scores
use only each row's own values plus the frozen model, with no
scored-batch centering, quantiles, or renormalization.

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
score_reo_metaktsp(model, X[1, , drop = FALSE])
} # }
```
