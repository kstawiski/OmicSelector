# Score specimens with a frozen trimmed-rCLR signed panel

Applies \[ws_rclr_trimmed()\] independently to each incoming specimen
and sums the transformed values over the training-frozen feature panel
after applying the training-frozen signs. Extra features are ignored.
Missing panel features fail closed rather than being imputed from a
scored batch.

## Usage

``` r
score_ws_rclr_panel(model, X, meta = NULL)
```

## Arguments

- model:

  A \`ws_rclr_panel_model\` returned by \[fit_ws_rclr_panel()\].

- X:

  Numeric specimens-by-features matrix.

- meta:

  Unused; accepted for the canonical roster score contract.

## Value

One finite numeric score per specimen.
