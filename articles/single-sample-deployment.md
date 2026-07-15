# Single-sample deployment with OmicSelector

``` r

library(OmicSelector)
#> OmicSelector 2.6.1 - Leakage-Resistant Biomarker Discovery
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
```

## Check the row-equivariance guarantee

The deployability gate checks the property needed for singleton use: a
specimen’s score is the same when computed alone as when computed as a
row in a batch, with the frozen model unchanged.

``` r

is_singlesample_deployable(dep, X_probe = X_train)
#> [1] TRUE
```

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
