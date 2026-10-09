# Score specimens with a frozen trimmed-rCLR signed panel

Applies
[`ws_rclr_trimmed`](https://kstawiski.github.io/OmicSelector/reference/ws_rclr_trimmed.md)
to each row of the entire supplied matrix before selecting the frozen
panel and summing its signed values. Background columns can change the
centering set and the default specimen-specific pseudocount, so adding
or removing them can change scores even when every panel feature is
present. Use a consistent feature universe for comparable scores.
Missing panel features always fail closed.

## Usage

``` r
score_ws_rclr_panel(model, X, meta = NULL, strict = FALSE)
```

## Arguments

- model:

  A `ws_rclr_panel_model` from
  [`fit_ws_rclr_panel`](https://kstawiski.github.io/OmicSelector/reference/fit_ws_rclr_panel.md).

- X:

  Numeric specimens-by-features matrix.

- meta:

  Unused; accepted for the canonical roster score contract.

- strict:

  Logical scalar. The default `FALSE` preserves the supplied- matrix
  transform and permits missing non-panel training columns for backward
  compatibility. Set `TRUE` to reject any missing column from the frozen
  training feature universe, including background columns. Extra columns
  are still transformed in both modes; strict mode checks absence only.

## Value

One finite numeric score per specimen, with a `support` data-frame
attribute containing one row per score. `n_available` counts training-
universe columns present, `n_fitted` counts all training-universe
columns, and `full_support` reports whether these counts agree. These
are feature counts, not pair counts or the number of features retained
for centering. They describe column availability, not detection or assay
quality. A TRUE flag does not rule out effects of extra background
columns. Save `attr(score, "support")` before subsetting or numeric
coercion, which can discard attributes.

## Examples

``` r
set.seed(7)
X <- matrix(rexp(20 * 8), 20, dimnames = list(NULL, paste0("f", 1:8)))
model <- fit_ws_rclr_panel(X, rep(0:1, each = 10),
                          hp = list(panel_size = 1, min_centering_size = 1,
                                    trim_upper = 0, trim_lower = 0))
full <- score_ws_rclr_panel(model, X, strict = TRUE)
panel_only <- X[, model$panel, drop = FALSE]
partial <- score_ws_rclr_panel(model, panel_only)
attr(partial, "support")  # 1 available / 8 fitted features; full_support = FALSE
#>    n_available n_fitted full_support
#> 1            1        8        FALSE
#> 2            1        8        FALSE
#> 3            1        8        FALSE
#> 4            1        8        FALSE
#> 5            1        8        FALSE
#> 6            1        8        FALSE
#> 7            1        8        FALSE
#> 8            1        8        FALSE
#> 9            1        8        FALSE
#> 10           1        8        FALSE
#> 11           1        8        FALSE
#> 12           1        8        FALSE
#> 13           1        8        FALSE
#> 14           1        8        FALSE
#> 15           1        8        FALSE
#> 16           1        8        FALSE
#> 17           1        8        FALSE
#> 18           1        8        FALSE
#> 19           1        8        FALSE
#> 20           1        8        FALSE
stopifnot(any(abs(as.numeric(full) - as.numeric(partial)) > 1e-8))
tryCatch(score_ws_rclr_panel(model, panel_only, strict = TRUE),
         singlesample_input_error = function(e) e$code)
#> [1] "MISSING_REQUIRED_FEATURES"
```
