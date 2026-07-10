# Flag samples whose dominance score exceeds a calibrated threshold

Convenience wrapper: returns a logical vector of length `nrow(x)`
indicating which samples have a dominance score exceeding `threshold`.
Use
[`fit_dominance_threshold`](https://kstawiski.github.io/OmicSelector/reference/fit_dominance_threshold.md)
to calibrate the threshold from a training reference cohort.

## Usage

``` r
ws_dominance_flag(x, top_k = 2L, threshold = 0.6, pseudocount = 0.5)
```

## Arguments

- x:

  Samples x features matrix.

- top_k:

  Integer; number of top features to include in the dominance numerator.
  Default 2.

- threshold:

  Numeric. Default 0.60.

- pseudocount:

  Additive pseudocount. Default 0.5.

## Value

Logical vector of length `nrow(x)`.
