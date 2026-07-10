# Score the GASF-image CNN single-sample discriminator (forced row-by-row)

Scores each row of `X` independently with the frozen model from
[`fit_img_gasfcnn`](https://kstawiski.github.io/OmicSelector/reference/fit_img_gasfcnn.md).
The python CNN is REBUILT from the exported R-side `state_dict`
(`load_state_dict`, [`eval()`](https://rdrr.io/r/base/eval.html) -\>
FROZEN BatchNorm, float64), then each query is mapped to the per-sample
robust CLR over the frozen feature universe (absent universe features
carry the neutral rCLR value `0`), rendered to its GASF image under the
FROZEN bounds, and forwarded through the CNN ONE ROW AT A TIME. The
score is the eval-mode logit; larger = more case-like. Forcing the
\\n=1\\ forward path uniformly (float64) makes a row's logit
bit-identical whether scored alone or in any batch.

Queries with fewer than `model$hp$min_features` universe features
present (a column-overlap floor), or empty positive support over the
frozen universe on the (pre-rCLR) aligned abundances
(`!any(X_use[i, ] > 0)`), return the neutral score `0`. A FLAT
all-equal-positive composition maps to the rCLR origin but is a VALID
specimen and is scored normally. The score of a row depends only on that
row and the frozen model, and is exactly invariant to per-specimen
positive scaling.

## Usage

``` r
score_img_gasfcnn(model, X, meta = NULL)
```

## Arguments

- model:

  An `img_gasfcnn_model` object from
  [`fit_img_gasfcnn`](https://kstawiski.github.io/OmicSelector/reference/fit_img_gasfcnn.md).

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

Wang Z, Oates T. (2015) Imaging Time-Series to Improve Classification
and Imputation. *IJCAI*; arXiv:1506.00327.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_img_gasfcnn(X, y)
score_img_gasfcnn(model, X[1, , drop = FALSE])
} # }
```
