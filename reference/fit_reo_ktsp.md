# Fit REO-kTSP within-sample pair-order discriminator

Learns a frozen k-Top-Scoring-Pairs (k-TSP) discriminator from training
samples by calling
[`os_ktsp_fit`](https://kstawiski.github.io/OmicSelector/reference/os_ktsp_fit.md).
The fitted primitive stores oriented feature pairs so that larger vote
fractions indicate the case class: each retained pair votes 1 for a
specimen when `feature_a < feature_b` within that specimen. No test data
or test-time batch statistics are used.

## Usage

``` r
fit_reo_ktsp(X_train, y_train, meta_train = NULL, hp = list())
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

  List of frozen hyperparameters: `k`, the number of oriented pairs
  requested from
  [`os_ktsp_fit`](https://kstawiski.github.io/OmicSelector/reference/os_ktsp_fit.md)
  (default `11L`). `k` is resolved once during fitting and is not tuned
  inside `fit_reo_ktsp()`.

## Value

A plain list of class `reo_ktsp_model` containing the fitted
`os_ktsp_model` in `ktsp`, `feature_universe`, retained pair count `k`,
and the resolved `hp`.

## References

Geman D, d'Avignon C, Naiman DQ, Winslow RL. (2004) Classifying gene
expression profiles from pairwise mRNA comparisons. *Statistical
Applications in Genetics and Molecular Biology* 3: Article19.

Tan AC, Naiman DQ, Xu L, Winslow RL, Geman D. (2005) Simple decision
rules for classifying human cancers from gene expression profiles.
*Bioinformatics* 21: 3896-3904.

## Examples

``` r
if (FALSE) { # \dontrun{
X <- matrix(stats::rgamma(40 * 12, shape = 2), nrow = 40,
            dimnames = list(NULL, paste0("miR-", seq_len(12))))
y <- rep(c(0, 1), each = 20)
X[y == 0, "miR-1"] <- X[y == 0, "miR-1"] + 5
X[y == 1, "miR-2"] <- X[y == 1, "miR-2"] + 5
model <- fit_reo_ktsp(X, y, hp = list(k = 3L))
score <- score_reo_ktsp(model, X)
} # }
```
