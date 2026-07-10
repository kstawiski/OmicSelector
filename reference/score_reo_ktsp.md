# Score REO-kTSP within-sample pair-order votes

Scores each row of `X` independently with the frozen oriented pairs
learned by
[`fit_reo_ktsp`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_ktsp.md).
For the pairs whose two features are present in `X`, the score is the
k-TSP vote fraction: the mean of `feature_a < feature_b` votes over
usable pairs. This is the standard k-TSP decision statistic and supplies
a continuous margin score for DeLong comparisons. If no retained pair
has both features present, the function returns the neutral vote
fraction `0.5` for every row.

## Usage

``` r
score_reo_ktsp(model, X, meta = NULL)
```

## Arguments

- model:

  A `reo_ktsp_model` object returned by
  [`fit_reo_ktsp`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_ktsp.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Numeric vector of one finite vote-fraction score per row of `X`. Scores
use only each row's own values plus the frozen model, with no
scored-batch centering, quantiles, or renormalization. On a fully
present feature set the result is identical to `predict(model$ktsp, X)`.

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
score_reo_ktsp(model, X[1, , drop = FALSE])
} # }
```
