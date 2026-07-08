# Fit a Hemolysis Pre-Filter

Tunes a Blondal-style sample rejection threshold using training controls
only. The threshold is selected by contrasting clean controls against
synthetic proxy-hemolyzed controls and maximizing balanced accuracy.

## Usage

``` r
fit_hemolysis_prefilter(
  X_train_controls,
  numerator = "MIMAT0000092",
  denominator = "MIMAT0006764",
  feature_names = NULL,
  sensitivity = .paper1_hemolysis_coefficients(),
  synthetic_strengths = c(0.5, 1, 2),
  synthetic_repeats = 3L,
  threshold_grid = seq(0.5, 0.995, length.out = 100),
  seed = 42L,
  floor_log2 = -8
)
```

## Arguments

- X_train_controls:

  Numeric matrix of healthy-control samples only.

- numerator:

  Name of the numerator feature. Defaults to \`"MIMAT0000092"\`.

- denominator:

  Name of the denominator feature. Defaults to \`"MIMAT0006764"\`.

- feature_names:

  Optional feature names. If omitted, \`colnames()\` are used.

- sensitivity:

  Named positive numeric vector of hemolysis-sensitivity coefficients.
  Defaults to the locked eLife-14 proxy coefficients.

- synthetic_strengths:

  Proxy perturbation strengths used during tuning.

- synthetic_repeats:

  Number of synthetic repeats per strength.

- threshold_grid:

  Quantiles used to generate threshold candidates.

- seed:

  Integer seed for synthetic perturbations.

- floor_log2:

  Pseudo-linear floor used by the synthetic perturbation.

## Value

An object of class \`"omicselector_hemolysis_prefilter"\`.
