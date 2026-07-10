# Score REO-singscore within-sample rank signature

Scores each row of `X` independently with the frozen signature learned
by
[`fit_reo_singscore`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_singscore.md).
For one sample, features in the frozen training universe that are
present in `X` are ranked in increasing abundance using
`rank(x, ties.method = "average") / length(x)`. The up term is the mean
fractional rank of present up features; the down term is the mean
fractional rank of present down features. If a signature direction has
no present features, it contributes the neutral expected fractional rank
`0.5`. With `use_down = TRUE`, the returned score is
`mean_rank(up) - mean_rank(down)`; otherwise it is `mean_rank(up)`.

## Usage

``` r
score_reo_singscore(model, X, meta = NULL)
```

## Arguments

- model:

  A `reo_singscore_model` object returned by
  [`fit_reo_singscore`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_singscore.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Numeric vector of one finite score per row of `X`. Scores use only each
row's own values plus the frozen model, with no scored-batch centering,
quantiles, or renormalization. If no model-universe features are shared
with `X`, scoring stops explicitly.

## References

Foroutan M, Bhuva DD, Lyu R, Horan K, Cursons J, Davis MJ. (2018) Single
sample scoring of molecular phenotypes. *BMC Bioinformatics* 19: 404.
PMID 30400809.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(40 * 30, shape = 2), nrow = 40,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 20)
X[y == 1, 1:5] <- X[y == 1, 1:5] + 5
model <- fit_reo_singscore(X, y)
score_reo_singscore(model, X[1, , drop = FALSE])
} # }
```
