# Dataset-Identity AUC ("Bias Floor")

Fits a logistic regression using ONLY cohort/dataset dummy variables as
predictors and reports its AUC against the case/control label. This is
the empirical lower bound on how much of an apparent classifier AUC can
be explained by cohort identity alone, before any biomarker is used.

If the bias floor is high (e.g., \> 0.70), any pooled classifier AUC on
the same data inherits that share of discrimination from cohort
provenance rather than biology. A responsible biomarker claim must show
AUC *materially* above the bias floor, on data where the floor itself
has been mitigated.

## Usage

``` r
os_bias_floor_auc(y, cohort, cv_folds = 5L, seed = 42L)
```

## Arguments

- y:

  Binary outcome vector (0/1 or factor with two levels). Length n.

- cohort:

  Vector or factor of cohort/dataset identity. Length n.

- cv_folds:

  Integer. If non-NULL, also computes K-fold cross-validated AUC.
  Default 5. Set to NULL to skip CV.

- seed:

  Random seed for CV folds. Default 42L.

## Value

A list with class `"os_bias_floor"`:

- apparent_auc:

  Logistic regression AUC on the training data (upper bound — optimistic
  by design)

- cv_auc:

  Cross-validated AUC, if requested

- per_cohort:

  Data frame: cohort, n, case rate, mean predicted probability

- interpretation:

  Character string summarizing the floor

## Examples

``` r
if (FALSE) { # \dontrun{
y <- c(rep(1, 317), rep(0, 200), rep(1, 40), rep(0, 100))
cohort <- c(rep("A", 317), rep("A", 200), rep("B", 40), rep("B", 100))
bf <- os_bias_floor_auc(y, cohort)
print(bf)
} # }
```
