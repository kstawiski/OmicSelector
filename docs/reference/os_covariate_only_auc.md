# Covariate-Only AUC

Fits a logistic regression using a single covariate (typically age or
sex) as the only predictor of case status. Quantifies how much of the
apparent classifier signal could be explained by demographic confounding
alone.

## Usage

``` r
os_covariate_only_auc(y, covariate, covariate_name = "covariate")
```

## Arguments

- y:

  Binary outcome vector.

- covariate:

  Numeric or factor covariate (e.g., age in years).

- covariate_name:

  Human-readable label for reporting.

## Value

A list:

- auc:

  AUC of covariate-only classifier

- direction:

  "higher values = cases" or "higher values = controls"

- effect_size:

  Cohen's d (numeric) or Cramer's V (factor)
