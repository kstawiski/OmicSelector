# Single-sample deployment with OmicSelector

``` r

library(OmicSelector)
#> OmicSelector 2.6.6 - Leakage-Resistant Biomarker Discovery
#> Use OmicPipeline$new() to start. Run ?OmicSelector for help.
```

## Motivation

> **Deployment problem.** Classical batch-aware corrections such as
> ComBat, quantile normalization, or RUV are designed around a batch or
> reference set. That is a mismatch for a clinical-style handoff where
> one new specimen arrives alone and must be scored without co-resident
> test samples.
>
> **Honest caveat.** The value of the API below is *deployability*: it
> scores one incoming specimen from its own panel values plus frozen
> fit-time parameters, with no co-resident test batch. It is not a claim
> that `ws-balance-ilr`, or any single-sample method, out-discriminates
> a batch-corrected pipeline — choose and validate the scoring method on
> your own data.

## Make a tiny synthetic compositional panel

This vignette uses a self-contained synthetic count-like matrix. The
feature names come from the default within-sample balance dictionary so
the example exercises the same panel shape as the deployment tests, but
no data are downloaded.

``` r

set.seed(101)

sbp <- ws_default_sbp()
sentinels <- c("ALL_NON_RBC", "ALL_NON_PLT",
               "TOP_K_BY_ABUNDANCE", "TAIL_BY_ABUNDANCE")
features <- unique(unlist(lapply(sbp, function(balance) {
  c(balance$numerator, balance$denominator)
}), use.names = FALSE))
features <- features[!features %in% sentinels]
features <- features[seq_len(30L)]

n <- 40L
y <- rep(c(0L, 1L), each = n / 2L)
X <- matrix(
  stats::rgamma(n * length(features), shape = 20, rate = 2),
  nrow = n,
  dimnames = list(paste0("S", seq_len(n)), features)
)
X[y == 1L, features[1:5]] <- X[y == 1L, features[1:5]] + 80
X[y == 0L, features[6:10]] <- X[y == 0L, features[6:10]] + 80

holdout <- nrow(X)
X_train <- X[-holdout, , drop = FALSE]
y_train <- y[-holdout]
x_new <- X[holdout, , drop = FALSE]

dim(X_train)
#> [1] 39 30
dim(x_new)
#> [1]  1 30
```

## Fit once, then freeze the scorer

[`deploy_singlesample()`](https://kstawiski.github.io/OmicSelector/reference/deploy_singlesample.md)
fits one rostered method and stores the frozen model state needed for
future singleton scoring. The default `ws-balance-ilr` method is a
compositional deployment default, not a certified benchmark winner.

``` r

dep <- deploy_singlesample(
  X_train,
  y_train,
  method = "ws-balance-ilr"
)

print(dep)
#> OmicSelector single-sample deployable
#>   method_id: ws-balance-ilr
#>   n_train: 39
#>   gate: verified row-equivariant scoring on 39 probe rows (single specimen == same row in batch; frozen model)
#>   caveat: deployability only; not a claim of superior discrimination over a batch-corrected pipeline.
```

## Score one incoming specimen

[`score_specimen()`](https://kstawiski.github.io/OmicSelector/reference/score_specimen.md)
needs only the held-out specimen row and the frozen deployment object.
No test batch, co-resident reference cohort, or test-time batch
correction is supplied.

``` r

score <- score_specimen(dep, x_new)
score
#> [1] 8.42726
#> attr(,"support")
#>   n_available n_fitted full_support
#> 1          30       30         TRUE
```

## Check the row-equivariance guarantee

The deployability gate checks the property needed for singleton use: a
specimen’s score is the same when computed alone as when computed as a
row in a batch, with the frozen model unchanged.

``` r

is_singlesample_deployable(dep, X_probe = X_train)
#> [1] TRUE
```

## Freeze conversion from raw counts or Ct values

A frozen preprocessor records the feature names and conversion
constants. Counts are divided by each specimen’s own total over
available frozen features and multiplied by one million. Ct values
become `2^(40 - Ct)` by default. Missing values, non-numeric Ct markers,
and Ct values at or above 40 become abundance 1 with detection FALSE.
Scoring rejects an all-non-detect specimen.

``` r

counts_train <- round(X_train)
counts_new <- round(x_new)
count_prep <- fit_singlesample_preprocessor(counts_train, type = "counts")
count_dep <- deploy_singlesample(counts_train, y_train, method = "reo-ktsp",
                                preprocessor = count_prep)
score_specimen(count_dep, counts_new)
#> [1] 1
#> attr(,"support")
#>   n_available n_fitted full_support
#> 1          11       11         TRUE
is_singlesample_deployable(count_dep, counts_train)
#> [1] TRUE

ct_train <- 40 - log2(X_train)
ct_new <- 40 - log2(x_new)
ct_prep <- fit_singlesample_preprocessor(ct_train, type = "ct")
ct_dep <- deploy_singlesample(ct_train, y_train, method = "reo-ktsp",
                             preprocessor = ct_prep)
score_specimen(ct_dep, ct_new)
#> [1] 1
#> attr(,"support")
#>   n_available n_fitted full_support
#> 1          11       11         TRUE
is_singlesample_deployable(ct_dep, ct_train)
#> [1] TRUE
```

Extra columns are dropped. Missing frozen columns remain absent; for
CPM, missing columns change the total used for conversion. Use
`strict = TRUE` to require the scorer’s complete fitted support. Any
supplied `detected` mask must match the raw input; it is combined with
the preprocessor’s mask by logical AND. The numerical checks above do
not establish clinical performance.

## Register a custom scorer

Registration checks that a specimen has the same finite score alone and
in different batches, and that scoring leaves the fitted model
unchanged. The example freezes the training mean of the first feature.
It is a software example, not a validated clinical score.

``` r

fit_custom <- function(X_train, y_train, ...) {
  list(panel = colnames(X_train)[1L], center = mean(X_train[, 1L]))
}
score_custom <- function(model, X, meta = NULL, ...) {
  as.numeric(X[, model$panel] - model$center)
}
register_singlesample_method(
  "example_centered_feature", fit_custom, score_custom,
  estimand = "within", probe = list(X = X_train, y = y_train)
)
registered_singlesample_methods()[, c("method_id", "row_source", "audited")]
#>                  method_id row_source audited
#> 1 example_centered_feature       user    TRUE
custom_dep <- deploy_singlesample(
  X_train, y_train, method = "example_centered_feature"
)
stopifnot(is_singlesample_deployable(custom_dep, X_train))
stopifnot(identical(score_specimen(custom_dep, x_new),
                    score_specimen(custom_dep, rbind(x_new, X_train))[1L]))
unregister_singlesample_method("example_centered_feature")
score_specimen(custom_dep, x_new) # Existing deployments retain their scorer.
#> [1] 39.19399
```

Routes are registered for the current R session. Use
`estimand = "within"` for within-cohort benchmark helpers. Strict
requests require a callback that accepts `strict` and registration with
`strict_support = TRUE`. Probe auditing does not establish performance
on other panels or clinical readiness.

## Inspect the method roster, not a ranking

The deployment wrapper is intentionally fail-closed. Roughly 30 rostered
single-sample methods are deployable through this fit-freeze-score
pattern, and methods needing cohort context, unavailable backends, or
non-scalar outputs reject cleanly rather than pretending to be
singleton-safe. Use
[`singlesample_method_roster()`](https://kstawiski.github.io/OmicSelector/reference/singlesample_method_roster.md)
to inspect the current roster and do not treat the roster order as a
performance ranking.

``` r

roster <- singlesample_method_roster()
head(roster[, c("method_id", "tier", "dep_route")], 8)
#>             method_id tier dep_route
#> 1     ws-rclr-trimmed   R1    base-R
#> 2        ws-alr-pivot   R1    base-R
#> 3     ws-mad-logratio   R1    base-R
#> 4        ws-dominance   R1    base-R
#> 5          frozen-ruv   R1    base-R
#> 6    robust-pca-resid   R1    base-R
#> 7     hemolysis-resid   R1    base-R
#> 8 kit-stratified-rclr   R1    base-R
```

The deployability contract is what the roster guarantees, not a ranking:
a deployable method can score one incoming specimen without a new batch
or a test-time batch-correction step. It does not certify that any
method provides superior discrimination, so benchmark the method choice
on your own cohorts.
