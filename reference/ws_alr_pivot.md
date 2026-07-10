# Additive log-ratio with a frozen pivot pool

Within-sample additive log-ratio (ALR) transformation using a frozen
pivot pool of housekeeping-like miRNAs as the geometric-mean
denominator. \$\$x\_{\mathrm{alr},i} = \log(x_i + c) -
\log(\mathrm{gmean}(x\_{\mathrm{pivot}} + c))\$\$ Anchors the
transformation against a stable internal reference rather than choosing
a single feature as denominator (the classical Aitchison ALR).

Fail-closed by default: refuses to compute when fewer than
`min_pivot_present` pivots are present in the input panel, unless the
operator explicitly opts in to a global rCLR fallback. Silent
substitution of all-feature centering for a curated denominator can mask
haemolysis or platelet contamination - flagged by codex Round 1 review.

## Usage

``` r
ws_alr_pivot(
  x,
  pivot_features = NULL,
  pseudocount = 0.5,
  allow_global_fallback = FALSE,
  min_pivot_present = 3L,
  verbose = FALSE
)
```

## Arguments

- x:

  Samples x features matrix (rows = samples). Must have column names to
  identify pivot features.

- pivot_features:

  Character vector of feature names to use as pivot pool. If `NULL`,
  defaults to
  [`ws_default_pivot_pool()`](https://kstawiski.github.io/OmicSelector/reference/ws_default_pivot_pool.md).

- pseudocount:

  Additive pseudocount. Default 0.5.

- allow_global_fallback:

  Logical. If `TRUE`, falls back to a global rCLR (all features) when
  fewer than `min_pivot_present` pivots are present. Default `FALSE`
  (fail-closed).

- min_pivot_present:

  Minimum pivot features required. Default 3.

- verbose:

  Logical. If `TRUE` and fallback is taken, emits a message.

## Value

Same shape as `x`; the ALR-transformed values, with attributes
`pivot_features`, `pseudocount`, and `fallback_used`.

## Details

**Backwards-compatibility / single-sample-pipeline note (v2.4.0).** The
package default is `allow_global_fallback = FALSE` (fail-closed): the
function errors out rather than silently substituting a global rCLR
centering when fewer than `min_pivot_present` pivots are available. This
is the conservative choice for end users running interactive analyses,
because it prevents haemolysis-contaminated cells from being silently
rescued by an all-feature centering that includes the very contaminants
they were meant to flag.

The single-sample pipeline (within-sample compositional methods)
describes the fallback as having been "activated" with the affected
cells explicitly flagged in the per-cell results. The pipeline opts in
via `allow_global_fallback = TRUE` on every call. This alignment is
intentional: the package keeps the safer default for ad-hoc use, while
the single-sample pipeline has the explicit logging discipline to track
every cell that took the fallback path.
