# Score specimens with a frozen Nystrom diffusion-map discriminator

Scores each row of `X` independently with a frozen `man_nystrom_model`.
For each specimen, the scorer computes a per-sample rCLR vector over the
model features present in `X`, derives a diffusion-map geometry from the
frozen landmarks restricted to those same present features, embeds that
specimen by its own Nystrom affinity row to the frozen landmarks,
applies frozen coordinate standardization, and evaluates the frozen
linear head.

## Usage

``` r
score_man_nystrom(model, X, meta = NULL)
```

## Arguments

- model:

  A `man_nystrom_model` object returned by
  [`fit_man_nystrom`](https://kstawiski.github.io/OmicSelector/reference/fit_man_nystrom.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. Larger values are more
case-like. If fewer than `model$hp$min_features` model features are
present, the documented neutral score `0` is returned for every row.

## Details

At inference, no statistic is estimated from the scored rows / test
batch: no test-batch centering, renormalization, bandwidth update, graph
update, or head refit is performed. When the present feature set is a
strict subset of the training universe, the landmark diffusion geometry
IS re-derived over the present columns – but solely from the frozen
training landmarks, never from any scored row – so row-equivariance is
preserved; at full feature overlap this re-derivation reproduces the
fit-time geometry exactly. Partial-overlap scoring is therefore an
approximate graceful degradation (the scored row's per-sample rCLR is
re-centred over the present features while the frozen landmark
coordinates are sub-columned, and the re-derived eigenbasis need not
align with the fit-time head), not an exact projection. The rCLR
representation is exactly invariant to multiplying one specimen by a
positive scalar because \\\log(c v_j)-mean_k\log(c
v_k)=\log(v_j)-mean_k\log(v_k)\\ over that specimen's own positive
present features.

## References

Coifman RR, Lafon S. (2006) Diffusion maps. *Applied and Computational
Harmonic Analysis* 21: 5-30.

Bengio Y, Paiement JF, Vincent P, Delalleau O, Le Roux N, Ouimet M.
(2004) Out-of-sample extensions for LLE, Isomap, MDS, Eigenmaps, and
Spectral Clustering. *Advances in Neural Information Processing Systems*
16.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(60 * 30, shape = 3), nrow = 60,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 30)
X[y == 1, 1:6] <- X[y == 1, 1:6] + 20
X[y == 0, 7:12] <- X[y == 0, 7:12] + 20
model <- fit_man_nystrom(X, y)
score_man_nystrom(model, X[1, , drop = FALSE])
} # }
```
