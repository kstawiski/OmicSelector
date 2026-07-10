# Score the CoDaCoRe stagewise log-ratio-balance discriminator (pure base R)

Scores each row of `X` independently with the frozen model from
[`fit_coda_codacore`](https://kstawiski.github.io/OmicSelector/reference/fit_coda_codacore.md),
in PURE base R (no python). Each query is aligned to the frozen feature
universe; for each frozen discrete balance \\b\\ the specimen's discrete
ILR balance is computed over its OWN strictly-positive support in \\A_b
/ B_b\\ (a side with no positive feature contributes the neutral value
`0`); the score is the frozen logistic linear predictor
\\\mathrm{intercept} + \sum_b w_b\\\mathrm{bal}\_b(x)\\. Larger = more
case-like.

Queries with fewer than `model$hp$min_features` universe features
present (a column-overlap floor), or empty positive support over the
frozen universe on the aligned ORIGINAL abundances
(`!any(X_use[i, ] > 0)`), return the neutral score `0`. A FLAT
all-equal-positive composition gives every balance value `0` but is a
VALID specimen, scored to the intercept (NOT floored). The score of a
row depends only on that row and the frozen model, uses no random
numbers and no scored-batch statistics, and is exactly invariant to
per-specimen positive scaling.

## Usage

``` r
score_coda_codacore(model, X, meta = NULL)
```

## Arguments

- model:

  A `coda_codacore_model` object from
  [`fit_coda_codacore`](https://kstawiski.github.io/OmicSelector/reference/fit_coda_codacore.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with named feature columns.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`; larger values are more
case-like.

## References

Gordon-Rodriguez E, Quinn TP, Cunningham JP. (2021) Learning sparse
log-ratios for high-throughput sequencing data. *Bioinformatics*
38(1):157-163.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_coda_codacore(X, y)
score_coda_codacore(model, X[1, , drop = FALSE])
} # }
```
