# Eligible methods for the single-sample selector

Returns the frozen tier-R1 method bank that can be considered by the
experimental single-sample selector. The bank contains only
within-cohort baselines and discriminators accepted by
\[deploy_singlesample()\]; transfer estimands, conditional routes,
negative controls, and reticulate routes are excluded.

## Usage

``` r
singlesample_selector_candidates(roster = singlesample_method_roster())
```

## Arguments

- roster:

  Method roster from \[singlesample_method_roster()\].

## Value

Character vector in frozen roster order.
