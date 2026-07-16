# Fit a training-frozen trimmed-rCLR signed-panel score

Fits the canonical trimmed-rCLR reference used by the OmicSelector
single-sample benchmark. Each training specimen is transformed
independently with \[ws_rclr_trimmed()\]. Features are ranked by the
absolute standardized case-control mean difference in the training data,
the first \`panel_size\` features are retained, and each feature
receives the sign of its training case-minus-control median difference.
No held-out value or outcome is used.

## Usage

``` r
fit_ws_rclr_panel(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric training matrix with specimens in rows and named features in
  columns.

- y_train:

  Binary training labels coded 0/1.

- meta_train:

  Unused; accepted for the canonical roster fit contract.

- hp:

  Optional named list. Supported fields are \`panel_size\` (default 20),
  \`pseudocount\`, \`trim_upper\` (0.10), \`trim_lower\` (0.05),
  \`exclude_features\`, \`zero_policy\` (\`"pseudocount"\`), and
  \`min_centering_size\` (8).

## Value

A \`ws_rclr_panel_model\` containing the frozen feature panel, signs,
training effects, and transform settings.

## See also

\[score_ws_rclr_panel()\], \[ws_rclr_trimmed()\]
