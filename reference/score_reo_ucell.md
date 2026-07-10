# Score REO-UCell within-sample rank signature

Scores each row of `X` independently with the frozen signature learned
by
[`fit_reo_ucell`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_ucell.md).
For one sample, features in the frozen training universe that are
present in `X` are ranked in decreasing abundance using
`rank(-x, ties.method = "average")`. For a signature \\S\\ of size
\\n\\, ranks are capped at `maxRank_eff + 1`, where
`maxRank_eff = min(model$maxRank, number of present universe features)`.
The UCell-style score is \$\$1 - \\ \sum\_{f \in S} \min(r_f,
maxRank_eff + 1) - n(n+1)/2 \\ / (n \times maxRank_eff).\$\$ If no up or
down signature features are present, that direction contributes 0. With
`use_down = TRUE`, the returned score is `UCell(up) - UCell(down)`;
otherwise it is `UCell(up)`.

Because both the ranking universe and `maxRank_eff` depend on which
model-universe features are present in `X`, scores are only comparable
across specimens evaluated on the same feature panel; score each cohort
in a single comparison with one consistent set of features.

## Usage

``` r
score_reo_ucell(model, X, meta = NULL)
```

## Arguments

- model:

  A `reo_ucell_model` object returned by
  [`fit_reo_ucell`](https://kstawiski.github.io/OmicSelector/reference/fit_reo_ucell.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Numeric vector of one finite score per row of `X`. Scores use only each
row's own values plus the frozen model, with no scored-batch centering,
quantiles, or renormalization.

## References

Andreatta M, Carmona SJ. (2021) UCell: Robust and scalable single-cell
gene signature scoring. *Computational and Structural Biotechnology
Journal* 19: 3796-3798.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(40 * 30, shape = 2), nrow = 40,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 20)
X[y == 1, 1:5] <- X[y == 1, 1:5] + 5
model <- fit_reo_ucell(X, y)
score_reo_ucell(model, X[1, , drop = FALSE])
} # }
```
