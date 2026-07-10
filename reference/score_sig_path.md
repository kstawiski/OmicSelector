# Score the truncated path-signature discriminator (pure base R, row-by-row)

Scores each row of `X` independently with the frozen model from
[`fit_sig_path`](https://kstawiski.github.io/OmicSelector/reference/fit_sig_path.md).
Each query is mapped to the per-sample robust CLR over the frozen
feature universe (absent universe features carry the neutral rCLR value
`0`), lifted to the time-augmented 2D path, its closed-form truncated
path signature computed in PURE BASE R ONE ROW AT A TIME, standardized
with the FROZEN training statistics, and passed through the FROZEN
ridge-logistic linear head. The score is \\\mathrm{intercept} + \sum_d
w_d \cdot \mathrm{std\\sig}\_d(x)\\; larger = more case-like. The
signature is a deterministic pure-R function of the row alone, so a
row's score is bit-identical whether scored alone or in any batch.

Queries with fewer than `model$hp$min_features` universe features
present (a column-overlap floor), or empty positive support over the
frozen universe on the (pre-rCLR) aligned ORIGINAL abundances
(`!any(X_use[i, ] > 0)`), return the neutral score `0`. A FLAT
all-equal-positive composition maps to the rCLR origin (all-zero \\v\\)
but is a VALID specimen scored normally (NOT floored): the time
coordinate still moves, so its path signature is non-degenerate. The
score of a row depends only on that row and the frozen model, and is
exactly invariant to per-specimen positive scaling.

## Usage

``` r
score_sig_path(model, X, meta = NULL)
```

## Arguments

- model:

  A `sig_path_model` object from
  [`fit_sig_path`](https://kstawiski.github.io/OmicSelector/reference/fit_sig_path.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with named feature columns.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`; larger values are more
case-like.

## References

Chevyrev I, Kormilitzin A. (2016) A primer on the signature method in
machine learning. *arXiv:1603.03788*.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_sig_path(X, y)
score_sig_path(model, X[1, , drop = FALSE])
} # }
```
