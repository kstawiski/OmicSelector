# Construct fixed eligibility-only complete-support panels

Builds complete method-by-unit panels without reading outcome or
performance columns. At each coverage threshold, methods are retained
when eligible in at least the threshold proportion of primary units;
units are then retained only when every retained method is eligible.
Rank summaries are prohibited when fewer than \`min_methods\` methods or
\`min_units\` units remain.

## Usage

``` r
singlesample_complete_support_panels(
  eligibility,
  roster = singlesample_method_roster(),
  thresholds = c(0.5, 0.8, 0.95),
  primary_units = NULL,
  method_col = "method_id",
  unit_col = "unit_id",
  eligible_col = "eligible",
  min_methods = 5L,
  min_units = 5L
)
```

## Arguments

- eligibility:

  Data frame containing method, unit, and eligibility columns.

- roster:

  Frozen method roster; defaults to \[singlesample_method_roster()\].
  Roster row order is preserved.

- thresholds:

  Coverage proportions. Defaults to 0.50, 0.80, and 0.95.

- primary_units:

  Optional ordered vector of primary unit identifiers. If \`NULL\`,
  first-appearance order in \`eligibility\` is used.

- method_col, unit_col, eligible_col:

  Column names in \`eligibility\`.

- min_methods, min_units:

  Minimum panel dimensions required to permit rank and top-k summaries.

## Value

A list with \`summary\`, \`method_membership\`, and \`unit_membership\`
data frames. Panel construction uses no columns beyond the identifiers
and Boolean eligibility flag.

## Details

The eligibility table must contain exactly one Boolean cell for every
within-method by primary-unit combination. Additional known transfer
methods or non-primary units are ignored.
