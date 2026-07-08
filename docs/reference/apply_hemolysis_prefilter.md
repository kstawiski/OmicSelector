# Apply a Hemolysis Pre-Filter

Scores new samples with a fit created by \[fit_hemolysis_prefilter()\]
and returns a keep/reject decision for each sample.

## Usage

``` r
apply_hemolysis_prefilter(fit, X_new, feature_names = NULL)
```

## Arguments

- fit:

  A fitted \`"omicselector_hemolysis_prefilter"\` object.

- X_new:

  Numeric vector or matrix of new samples.

- feature_names:

  Optional feature names for \`X_new\`.

## Value

A \`data.frame\` with \`hemolysis_score\`, \`keep\`, \`rejected\`, and
\`threshold\`.
