# Block-aware two-stage BH-FDR for specimen-shared cohort clusters

Naïve BH-FDR within a modality family treats specimen-shared cohorts
(e.g., Toray clusters with \\\sim 94.5\\\\ specimen reuse) as
independent tests, which is anti-conservative due to strongly correlated
test statistics. This function implements a two-stage correction:

1.  Within each provenance block, combine cohort-level p-values via
    Fisher's chi-square statistic. Two reference distributions are
    reported side by side:

    - *conservative reference*: \\S/c\\ against \\\chi^2\_{2k}\\;
      intentionally conservative; used for the manuscript's headline
      counts.

    - *textbook reference* (Brown 1975 / Kost-McDermott 2002): \\S/c\\
      against \\\chi^2\_{2k/c}\\ with \\df\\ rescaled by the same
      inflation factor.

2.  Apply BH-FDR across block-level p-values, separately for each
    reference.

The inflation factor is \\c(k) = 1 + (k - 1)\rho\\ with \\k\\ the number
of cohorts in the block and \\\rho\\ the assumed within-block
correlation (default \\\rho = 0.25\\ reproduces the manuscript's
headline).

## Usage

``` r
singlesample_bh_fdr_correct_blocked(results, block_id, rho = 0.25)
```

## Arguments

- results:

  List of
  [`singlesample_matched_null_benchmark`](https://kstawiski.github.io/OmicSelector/reference/singlesample_matched_null_benchmark.md)
  return objects.

- block_id:

  Character vector — provenance block identifier per cohort (e.g.,
  `"Toray-cluster"`, `"Affy-singleton"`).

- rho:

  Numeric in \[0, 1\]. Assumed within-block correlation; controls the
  inflation factor \\c(k) = 1 + (k - 1)\rho\\. Default 0.25 reproduces
  the manuscript headline. Set to 0 to recover vanilla Fisher.

## Value

A `data.frame` with columns: `cohort_idx`, `block_id`, `p_emp`,
`p_block` (conservative), `p_block_textbook`, `q_block_BH`
(conservative), `q_block_BH_textbook`, `q_cohort_within_block`.
