# Score sliced random-projection Gaussian LRT discriminator

Scores each row of `X` independently with a model learned by
[`fit_ot_slicedlrt`](https://kstawiski.github.io/OmicSelector/reference/fit_ot_slicedlrt.md).
Each specimen is transformed to its own rCLR vector over the features
shared with the frozen training universe, projected onto deterministic
random unit directions, and assigned the average per-direction Gaussian
log-likelihood ratio of the case class versus the control class. Larger
values are more case-like.

## Usage

``` r
score_ot_slicedlrt(model, X, meta = NULL)
```

## Arguments

- model:

  An `ot_slicedlrt_model` object returned by
  [`fit_ot_slicedlrt`](https://kstawiski.github.io/OmicSelector/reference/fit_ot_slicedlrt.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be uniquely named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. If fewer than
`model$hp$min_features` model-universe features are present in `X`,
returns a neutral vector of zeros.

## Details

At scoring time, the present feature set is
`intersect(model$feature_universe, colnames(X))`. The sliced Gaussian
representation is then re-derived exactly once from frozen case/control
anchors restricted to that present set plus the frozen seed. A scored
row contributes only its own rCLR vector; no statistic is estimated from
any other row of `X`. This row-local computation is compatible with
singleton deployment and is invariant to row order, duplication, and
batch composition.

The default `hp$min_features` is `3L`. Note that a present overlap of
exactly one feature is degenerate: the single-feature per-sample rCLR is
identically zero (\\\log v - \log v = 0\\), so every projection is zero
and the score is a constant neutral `0` with no discrimination. Keep
`min_features` \\\ge 2\\ (the default `3L` already does) for a
discriminating configuration.

For one row with projection coordinates \\s_m\\, case moments
\\\mu\_{1m}, \sigma^2\_{1m}\\, and control moments \\\mu\_{0m},
\sigma^2\_{0m}\\, the score is \$\$\frac{1}{M} \sum_m
\left\[-\frac{(s_m-\mu\_{1m})^2}{2\sigma^2\_{1m}}
-\frac{1}{2}\log\sigma^2\_{1m}
+\frac{(s_m-\mu\_{0m})^2}{2\sigma^2\_{0m}}
+\frac{1}{2}\log\sigma^2\_{0m}\right\].\$\$

## References

Bonneel N, Rabin J, Peyre G, Pfister H. (2015) Sliced and Radon
Wasserstein barycenters of measures. *Journal of Mathematical Imaging
and Vision* 51: 22-45.

Neyman J, Pearson ES. (1933) On the problem of the most efficient tests
of statistical hypotheses. *Philosophical Transactions of the Royal
Society of London. Series A* 231: 289-337.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(80 * 30, shape = 2), nrow = 80,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 40)
X[y == 1, 1:6] <- X[y == 1, 1:6] + 25
X[y == 0, 7:12] <- X[y == 0, 7:12] + 25
model <- fit_ot_slicedlrt(X, y)
score_ot_slicedlrt(model, X[1, , drop = FALSE])
} # }
```
