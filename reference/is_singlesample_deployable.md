# Check whether a deployment scorer is single-sample deployable

User-facing wrapper over
[`singlesample_is_row_equivariant`](https://kstawiski.github.io/OmicSelector/reference/singlesample_is_row_equivariant.md).
Supply representative probe rows in `X_probe`; the check asks whether
each specimen's singleton score equals the same row's score inside a
batch, with the model state unchanged during scoring.

## Usage

``` r
is_singlesample_deployable(deployable, X_probe = NULL)
```

## Arguments

- deployable:

  A `singlesample_deployable` object from
  [`deploy_singlesample`](https://kstawiski.github.io/OmicSelector/reference/deploy_singlesample.md).

- X_probe:

  Numeric matrix/data frame of probe rows. If `NULL`, the check cannot
  be run; callers should supply probe rows from the intended panel or
  validation surface.

## Value

`TRUE` when the probe check passes, otherwise `FALSE`.
