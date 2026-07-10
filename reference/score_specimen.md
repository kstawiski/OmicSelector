# Score new specimens with a frozen single-sample deployment object

Scores one or more incoming specimens using only `x_new` plus the frozen
parameters stored in `deployable`. No scored batch, co-resident
reference cohort, or test-time batch statistic is required.

## Usage

``` r
score_specimen(deployable, x_new, meta = NULL)
```

## Arguments

- deployable:

  A `singlesample_deployable` object from
  [`deploy_singlesample`](https://kstawiski.github.io/OmicSelector/reference/deploy_singlesample.md).

- x_new:

  A single named numeric vector, or a numeric matrix/data frame with one
  or more specimen rows and feature columns.

- meta:

  Optional per-specimen metadata, one row per row of `x_new`.

## Value

Numeric score vector with one value per incoming specimen.
