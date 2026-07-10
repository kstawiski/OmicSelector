# Score MFDFA spectrum (frac-mfdfa) within-sample discriminator

Scores each row of `X` independently with the frozen model learned by
[`fit_frac_mfdfa`](https://kstawiski.github.io/OmicSelector/reference/fit_frac_mfdfa.md).
Model-universe features present in `X` are placed in the frozen
canonical order; for one specimen its within-sample rCLR landscape is
resampled to the frozen length, transformed into a fixed MFDFA spectrum
descriptor using the frozen scales and q grid, standardized by the
frozen training center/scale, and passed through the frozen ridge-LDA
head. Larger scores are more case-like.

Scoring uses only each row's own values plus the frozen model. If fewer
than `model$hp$min_features` model-universe features are present in `X`,
or a specimen's own positive support is below that floor, the documented
neutral score `0` is returned.

## Usage

``` r
score_frac_mfdfa(model, X, meta = NULL)
```

## Arguments

- model:

  A `frac_mfdfa_model` returned by
  [`fit_frac_mfdfa`](https://kstawiski.github.io/OmicSelector/reference/fit_frac_mfdfa.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. Larger values are more
case-like.

## References

Kantelhardt JW, Zschiegner SA, Koscielny-Bunde E, et al. (2002)
Multifractal detrended fluctuation analysis of nonstationary time
series. *Physica A* 316(1-4): 87-114.

Peng CK, Buldyrev SV, Havlin S, et al. (1994) Mosaic organization of DNA
nucleotides. *Physical Review E* 49(2): 1685-1689.

## Examples

``` r
if (FALSE) { # \dontrun{
score_frac_mfdfa(model, X[1, , drop = FALSE])
} # }
```
