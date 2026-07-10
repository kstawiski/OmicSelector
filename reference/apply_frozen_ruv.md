# Apply frozen RUV correction to new samples

Removes the unwanted-variation component estimated by
[`fit_frozen_ruv`](https://kstawiski.github.io/OmicSelector/reference/fit_frozen_ruv.md)
from one or more new samples. The latent factors are approximated from
the frozen pseudo-inverse learned during training.

## Usage

``` r
apply_frozen_ruv(x, fit)
```

## Arguments

- x:

  Numeric vector (single sample) or matrix (samples \\\times\\
  features).

- fit:

  A `frozen_ruv_fit` object from
  [`fit_frozen_ruv`](https://kstawiski.github.io/OmicSelector/reference/fit_frozen_ruv.md).

## Value

Same shape as `x`, with the unwanted-variation component removed and
training feature means re-added. The attribute `fit_summary` summarises
the fit parameters.
