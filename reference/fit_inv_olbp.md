# Fit ordinal Local Binary Pattern (inv-olbp) within-sample discriminator

Learns a frozen, single-sample discriminator from the *ordinal texture*
of a specimen's abundance profile. The construction has three frozen
stages, all estimated from training data only:

1.  **Canonical feature order.** Features are arranged once, in
    *decreasing* training column mean of `X_train` (ties broken by
    feature index, so the order is deterministic and
    locale-independent). This fixed permutation, `order_perm`, defines a
    circular sequence of feature positions shared by every specimen.

2.  **Ordinal Local Binary Pattern (OLBP) histogram.** For one specimen
    its values are read off in the canonical order and converted to
    within-sample ranks (`rank(.,"average")`). At each position the
    centre rank is compared to its `P` circular neighbours (offsets
    `c(-(P/2):-1, 1:(P/2))`); each comparison contributes one bit (`1`
    when the centre rank is \\\ge\\ the neighbour rank), and the `P`
    bits are packed – in the listed offset order, `offsets[k]` carrying
    weight \\2^{k-1}\\ – into an integer code in \\\[0, 2^P-1\]\\. The
    length-\\2^P\\ normalised histogram of these codes over all
    positions is the specimen's texture descriptor. Because it is built
    from within-sample ranks it is invariant to any per-sample
    monotone-increasing transform (library size, input amount, log/qPCR
    vs microarray vs NGS scaling).

3.  **Frozen linear head.** Per-sample histograms of the training matrix
    are standardised by a frozen `center`/`scale` (each bin's training
    mean and standard deviation, the latter floored at `eps`;
    zero-variance bins are deactivated and contribute `0`). A logistic
    regression (`glm.fit`, binomial) of the labels on the standardised
    active bins gives the head coefficients; if that fit errors, fails
    to converge, returns non-finite coefficients, or exhibits complete
    or quasi-complete separation (fitted probabilities pinned at 0 or 1,
    which would otherwise freeze a numerically unstable divergent fit),
    a diagonal-LDA mean-difference direction (case-minus-control mean on
    the standardised active bins, thresholded at the class midpoint) is
    used instead. Both orient *larger = more case-like*.

Disease reshapes *which* miRNAs dominate relative to their canonical
neighbours; the OLBP histogram is a compact, batch-robust summary of
that within-sample shape. No test data and no test-time batch statistic
are used: the feature order, neighbour offsets, histogram binning,
standardisation, and head are all frozen here from training only.

## Usage

``` r
fit_inv_olbp(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be uniquely named feature ids.

- y_train:

  Integer/numeric 0/1 labels aligned to `X_train`; 1 is case/disease and
  0 is control.

- meta_train:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

- hp:

  List of frozen hyperparameters: `neighbors` the number of circular
  OLBP neighbours \\P\\ (default `4L`; must be an even integer in
  `[2, 8]`, giving \\2^P\\ histogram bins), `min_features` the
  feature-overlap floor at scoring (default `max(neighbors + 1L, 4L)`,
  positive integer), and `eps` the positive standard-deviation floor
  applied to active histogram bins in the head (only genuinely
  zero-variance bins are deactivated and contribute `0`; default
  `1e-8`). These are resolved once during fitting and are not tuned
  inside `fit_inv_olbp()`. No randomness is used, so fitting is
  deterministic.

## Value

A plain list of class `inv_olbp_model` containing `feature_universe`,
the frozen `order_perm` (feature names in canonical order),
`neighbor_offsets`, `n_codes` (\\= 2^P\\), the frozen `head` (`type` –
`"glm"` or `"diag_lda"` – `center`, `scale`, `active`, `coef`,
`intercept`), and the resolved `hp`.

## Details

The per-specimen score (see
[`score_inv_olbp`](https://kstawiski.github.io/OmicSelector/reference/score_inv_olbp.md))
is the frozen linear predictor on the standardised OLBP histogram,
\$\$S(x) = \beta_0 + \sum\_{b} \beta_b \\ \frac{h_b(x) - c_b}{s_b},\$\$
where \\h(x)\\ is the specimen's normalised code histogram, \\c, s\\ are
the frozen bin centres/scales, and \\\beta\\ the frozen head. Larger
\\S\\ means more case-like. Every quantity except \\h(x)\\ is frozen at
fit time, and \\h(x)\\ depends only on \\x\\'s own ranks, so the score
is computable from a single specimen and is unchanged by per-sample
scaling.

## References

Ojala T, Pietikainen M, Maenpaa T. (2002) Multiresolution gray-scale and
rotation invariant texture classification with local binary patterns.
*IEEE Transactions on Pattern Analysis and Machine Intelligence* 24(7):
971-987.

Geman D, d'Avignon C, Naiman DQ, Winslow RL. (2004) Classifying gene
expression profiles from pairwise mRNA comparisons. *Statistical
Applications in Genetics and Molecular Biology* 3: Article19.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
p <- 40
X <- matrix(stats::rgamma(120 * p, shape = 20, rate = 0.2), nrow = 120,
            dimnames = list(NULL, paste0("miR-", sprintf("%03d", seq_len(p)))))
y <- rep(c(0, 1), each = 60)
# cases elevate a mid block, reordering the within-sample profile
X[y == 1, 11:16] <- X[y == 1, 11:16] + 150
model <- fit_inv_olbp(X, y)
score <- score_inv_olbp(model, X)
} # }
```
