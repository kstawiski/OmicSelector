# Score new specimens with a frozen single-sample deployment object

Scores one or more incoming specimens using only `x_new` plus the frozen
parameters stored in `deployable`. No scored batch, co-resident
reference cohort, or test-time batch statistic is required.

## Usage

``` r
score_specimen(deployable, x_new, meta = NULL, detected = NULL, strict = FALSE)
```

## Arguments

- deployable:

  A `singlesample_deployable` object from
  [`deploy_singlesample`](https://kstawiski.github.io/OmicSelector/reference/deploy_singlesample.md).

- x_new:

  A single named numeric vector, or a numeric matrix/data frame with one
  or more specimen rows and feature columns. With an attached
  preprocessor, supply raw counts or Ct values; Ct character markers are
  accepted and converted to non-detects.

- meta:

  Optional per-specimen metadata, one row per row of `x_new`. Legacy
  metadata routes require their training metadata field (kit, biofluid,
  technology or rin). Missing values raise \`MISSING_METADATA\`; unseen
  categorical levels raise \`UNSEEN_METADATA_LEVEL\`. Hemolysis
  residuals accept a finite \`hemo_score\` or require both frozen
  hemolysis anchors.

- detected:

  Optional logical matrix aligned to `x_new` (or a named logical vector
  for one specimen). TRUE means detected before imputation. Supply this
  mask when preprocessing replaces nondetects by positive values that
  cannot be distinguished from measured values. No detection status is
  inferred from scored-batch variation or from a constant specimen
  profile. With a preprocessor, align this mask to the raw input,
  including extra columns. It is combined with preprocessing detection
  by logical AND.

- strict:

  Logical scalar; defaults to `FALSE`. If `TRUE`, k-TSP and pair-ratio
  require every fitted-pair feature, trimmed rCLR requires the frozen
  training feature universe, and balance ILR requires every frozen
  balance member. Other routes and custom adapters reject with
  `STRICT_UNSUPPORTED`. Existing input guards apply in both modes.
  Registered custom scorers declaring \`strict_support = TRUE\` receive
  the flag and apply their documented strict policy.

## Value

Numeric score vector with one value per incoming specimen. For k-TSP,
pair-ratio, trimmed-rCLR and balance-ILR panels, the scorer's `support`
attribute is retained; see
[`score_reo_ktsp`](https://kstawiski.github.io/OmicSelector/reference/score_reo_ktsp.md),
[`score_reo_pairratio`](https://kstawiski.github.io/OmicSelector/reference/score_reo_pairratio.md),
[`score_ws_rclr_panel`](https://kstawiski.github.io/OmicSelector/reference/score_ws_rclr_panel.md)
and
[`score_ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/score_ws_balance_ilr.md).
The default preserves each scorer's missing-feature policy. Use
`strict = TRUE` here to require complete support. Invalid strict values
raise `INVALID_STRICT`; missing features raise
`MISSING_REQUIRED_FEATURES`, and absent frozen universe metadata raises
`MISSING_FITTED_UNIVERSE`.

## Details

Invalid input raises a `singlesample_input_error` with machine-readable
`code` and affected `rows`. Codes include `EMPTY_INPUT`,
`INVALID_FEATURE_NAMES`, `ALL_REQUIRED_ABSENT`, `ALL_NONDETECT`,
`NONFINITE_REQUIRED_VALUE`, `EMPTY_FITTED_SUPPORT`,
`CONSTANT_FITTED_SUPPORT`, and `NO_USABLE_SUPPORT`. Required NA, NaN or
infinite abundance values are rejected before scoring in both default
and strict modes. Raw Ct missing values are converted to non-detects
when a Ct preprocessor is attached. A batch containing an invalid row is
rejected as a whole. Checks use only each specimen and frozen training
information. Positive imputation cannot be recognized without `detected`
or an appropriate frozen `nondetect_limit`. Refit legacy models to
record training support status.
