# Isometric log-ratio balances on a frozen partition tree

Computes orthonormal isometric log-ratio (ILR) balances on a frozen
sequential binary partition. For each partition (G+, G-) the balance is
\$\$b = \sqrt{\frac{\|G+\|\\\|G-\|}{\|G+\| + \|G-\|}}
\log\\\left(\frac{\mathrm{gmean}(x\_{G+})}{\mathrm{gmean}(x\_{G-})}\right).\$\$
The default partition
([`ws_default_sbp`](https://kstawiski.github.io/OmicSelector/reference/ws_default_sbp.md))
encodes circulating miRNA biology: erythrocyte-vs-rest,
platelet-vs-rest, let-7 family within-cluster, miR-17~92, miR-200,
miR-371~373, and a dominant-vs-tail abundance balance computed at
runtime.

Tolerant to feature absence: features missing from the panel are simply
removed from the numerator/denominator before computation; a balance
becomes NA only if either side becomes empty.

## Usage

``` r
ws_balance_ilr(
  x,
  balances = ws_default_sbp(),
  pseudocount = NULL,
  aggregate = c("gmean", "trimmed_gmean"),
  min_balance_coverage = 0.8
)
```

## Arguments

- x:

  Numeric vector with miRNA names, OR a matrix (samples x features).
  Must have feature names - names are used to match the partition
  entries.

- balances:

  Named list of balances, each with components `numerator` and
  `denominator`. Defaults to
  [`ws_default_sbp()`](https://kstawiski.github.io/OmicSelector/reference/ws_default_sbp.md).

- pseudocount:

  Additive pseudocount before logging. Default 1e-6 of sample sum.

- aggregate:

  "gmean" (default) or "trimmed_gmean" (10% trim each end).

- min_balance_coverage:

  Numeric in \[0, 1\]. Minimum fraction of the eight balances that must
  be computable from the input panel for the sample / cohort to be
  deemed eligible. When fewer than this fraction of balances have both
  numerator and denominator features present, the returned vector /
  matrix is filled with NA and carries the attribute `coverage_failed`
  set to TRUE. Default 0.8 (manuscript Methods, 'Within-sample
  compositional methods').

## Value

Named numeric vector of balance values (single sample), or matrix
samples x balances (matrix input). Carries attributes `coverage`
(numeric in \[0, 1\]) and `coverage_failed` (logical) indicating whether
the 80% gate was satisfied.

## Examples

``` r
if (FALSE) { # \dontrun{
x <- rlnorm(60, meanlog = 5, sdlog = 1.2)
names(x) <- ws_default_pivot_pool()[1:6]
ws_balance_ilr(x)  # NA where biology dictionary names don't match panel
} # }
```
