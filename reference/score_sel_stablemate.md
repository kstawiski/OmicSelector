# Score StableMate stable-predictor discriminator (single-sample)

Scores each row of `X` independently with the frozen stable-feature
logistic learned by
[`fit_sel_stablemate`](https://kstawiski.github.io/OmicSelector/reference/fit_sel_stablemate.md).
For one specimen its per-sample rCLR is computed over the frozen
candidate universe (present features only), restricted to the frozen
selected stable features, and the frozen logistic linear predictor
\$\$\eta = intercept + \sum\_{j \in S} \beta_j \\ z_j\$\$ is returned
(NOT [`plogis()`](https://rdrr.io/r/stats/Logistic.html) of it); larger
is more case-like. Here \\z_j\\ is the specimen's own rCLR coordinate
for selected feature \\j\\. Because rCLR centres each specimen by its
own geometric mean over present parts, the score is EXACTLY invariant to
per-sample positive scaling and uses only that specimen's own values, so
the method is single-sample deployable and row-equivariant. Scoring uses
no random numbers and no scored-batch statistics.

Missing selected features contribute `0` (their rCLR coordinate is
absent from the specimen): the linear predictor is summed only over
selected features present in `X`. If fewer than `model$min_selected`
selected features are present in a specimen, or the model selected
nothing at fit time (`model$degenerate`), that specimen receives the
documented neutral score `0`.

## Usage

``` r
score_sel_stablemate(model, X, meta = NULL)
```

## Arguments

- model:

  A `sel_stablemate_model` object returned by
  [`fit_sel_stablemate`](https://kstawiski.github.io/OmicSelector/reference/fit_sel_stablemate.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method (the transfer/environment information is used
  only at fit time).

## Value

Plain finite numeric vector of length `nrow(X)`. Scores use only each
row's own values plus the frozen model, with no scored-batch centering,
quantiles, or renormalization.

## References

Deng Y, Liao Y, Wong NKP, Bhuva DD, Lin Y, Cao Y. (2024) StableMate: a
statistical method to select stable predictors in omics data. *bioRxiv*.

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
score_sel_stablemate(model, X[1, , drop = FALSE])
} # }
```
