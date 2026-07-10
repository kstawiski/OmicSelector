# Create an MDE/Margin Record

Create an MDE/Margin Record

## Usage

``` r
os_mde_margin_record(
  n,
  n_clusters,
  method_ids,
  null_draws,
  mde80,
  w90,
  floor = 0.02
)
```

## Arguments

- n:

  Effective sample count.

- n_clusters:

  Effective cluster count.

- method_ids:

  Character vector defining the fixed family.

- null_draws:

  Number of null draws.

- mde80:

  Minimum detectable effect at 80 percent power.

- w90:

  Width of the relevant 90 percent confidence interval.

- floor:

  Minimum floor margin.

## Value

List with class `"os_mde_margin_record"`.
