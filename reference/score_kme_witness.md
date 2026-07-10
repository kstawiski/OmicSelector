# Score kernel mean-embedding witness-at-a-point discriminator

Scores each row of `X` independently with the frozen model learned by
[`fit_kme_witness`](https://kstawiski.github.io/OmicSelector/reference/fit_kme_witness.md).
For one specimen \\x\\, the score is the witness function evaluated at
that point: \$\$mean_i k(phi(x), phi(a_i^{case})) - mean_j k(phi(x),
phi(a_j^{control})),\$\$ where \\k(a,b) = exp(-gamma \|\|a-b\|\|^2)\\.
At scoring time, the present feature set is
`intersect(model$feature_universe, colnames(X))`. Both the specimen and
the frozen raw anchors are closed to proportions over exactly that
present set and then centred-log-ratio transformed, which keeps
partial-overlap distances consistent (and the score exactly invariant to
per-sample scaling) without using any scored-batch statistic. The frozen
bandwidth `model$gamma` is applied unchanged in this reduced
present-feature space; it is not re-tuned to the overlap, which is what
keeps scoring single-sample. Re-deriving the anchor representation over
`present` uses only frozen training data, so the score depends only on
the specimen and the frozen model. If fewer than `model$hp$min_features`
features are present, the documented neutral score `0` is returned for
every row.

## Usage

``` r
score_kme_witness(model, X, meta = NULL)
```

## Arguments

- model:

  A `kme_witness_model` object returned by
  [`fit_kme_witness`](https://kstawiski.github.io/OmicSelector/reference/fit_kme_witness.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. Larger values are more
case-like because the witness is case mean embedding minus control mean
embedding. Scoring uses only each row's own values plus the frozen
anchors, `gamma`, and `pc`.

## References

Gretton A, Borgwardt KM, Rasch MJ, Scholkopf B, Smola A. (2012) A kernel
two-sample test. *Journal of Machine Learning Research* 13: 723-773.

Sutherland DJ, Tung HY, Strathmann H, De S, Ramdas A, Smola A, Gretton
A. (2017) Generative models and model criticism via optimized maximum
mean discrepancy. *International Conference on Learning
Representations*.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(60 * 30, shape = 2), nrow = 60,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 30)
X[y == 1, 1:5] <- X[y == 1, 1:5] + 20
X[y == 0, 6:10] <- X[y == 0, 6:10] + 20
model <- fit_kme_witness(X, y)
score_kme_witness(model, X[1, , drop = FALSE])
} # }
```
