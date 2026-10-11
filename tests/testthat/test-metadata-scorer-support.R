.metadata_support_data <- function() {
  set.seed(101)
  features <- c(.ta_default_pivot_pool(), paste0('F', 1:12))
  X <- matrix(rexp(80 * length(features)) + 1, 80,
              dimnames = list(NULL, features))
  y <- rep(0:1, each = 40)
  X[y == 1, 7:9] <- X[y == 1, 7:9] + 3
  list(X = X, y = y, meta = data.frame(
    technology = rep(c('NGS', 'Toray', 'qPCR', 'NanoString'), 20),
    biofluid = rep(c('serum-clot', 'plasma-EDTA'), each = 2, length.out = 80)))
}

test_that('unsupported metadata scorers raise input codes, not NA score vectors', {
  d <- .metadata_support_data()
  cases <- list(
    list(id = 'cross-tech-harmonized-rclr', code = 'INSUFFICIENT_CROSS_TECH_ANCHORS'),
    list(id = 'biofluid-anchor-rclr', code = 'INSUFFICIENT_BIOFLUID_ANCHORS'),
    list(id = 'tech-residualized-alr', code = 'INSUFFICIENT_ALR_PIVOTS'))
  for (case in cases) {
    X <- d$X; meta <- d$meta
    if (case$id == 'cross-tech-harmonized-rclr') meta$technology <- 'NGS'
    if (case$id == 'biofluid-anchor-rclr') meta$biofluid <- 'serum-clot'
    if (case$id == 'tech-residualized-alr') X <- X[, -c(3:6), drop = FALSE]
    err <- tryCatch(deploy_singlesample(X, d$y, case$id, meta_train = meta,
                                        verify = FALSE), error = identity)
    expect_s3_class(err, 'singlesample_input_error')
    expect_identical(err$code, case$code, info = case$id)
  }
})

test_that('valid fitted support stays finite and specimen-local', {
  d <- .metadata_support_data()
  for (id in c('cross-tech-harmonized-rclr', 'biofluid-anchor-rclr',
               'tech-residualized-alr')) {
    dep <- deploy_singlesample(d$X, d$y, id, meta_train = d$meta, verify = FALSE)
    before <- singlesample_state_fingerprint(dep)
    batch <- score_specimen(dep, d$X[1:5, ], d$meta[1:5, ])
    one <- vapply(1:5, function(i) as.numeric(score_specimen(
      dep, d$X[i, , drop = FALSE], d$meta[i, , drop = FALSE])), numeric(1))
    expect_true(all(is.finite(batch)), info = id)
    expect_equal(as.numeric(batch), one, tolerance = 1e-10, info = id)
    expect_identical(before, singlesample_state_fingerprint(dep), info = id)
    bad <- d$X[1, , drop = FALSE]; bad[1, 1] <- Inf
    err <- tryCatch(score_specimen(dep, bad, d$meta[1, , drop = FALSE]), error = identity)
    expect_identical(err$code, 'NONFINITE_REQUIRED_VALUE', info = id)
  }
})

test_that('direct predictors reject lost frozen anchors and pivots', {
  d <- .metadata_support_data(); panel <- colnames(d$X)[7:9]
  m <- d$meta; m$disease <- d$y; m$cohort <- 'training'
  cross <- fit_cross_tech_harmonized_rclr(d$X, m, panel)
  bio <- fit_biofluid_anchor_rclr(d$X, m, panel)
  tech <- fit_tech_residualized_alr(d$X, m, panel)
  cross$anchor_features <- cross$anchor_features[1:2]
  bio$anchor_features <- bio$anchor_features[1:2]
  tech$pivot_features <- tech$pivot_features[1:2]
  code <- function(expr) tryCatch(force(expr), error = function(e) e$code)
  expect_identical(code(predict_cross_tech_harmonized_rclr(cross, d$X[1:2, ], m[1:2, ])),
                   'INSUFFICIENT_CROSS_TECH_ANCHORS')
  expect_identical(code(predict_biofluid_anchor_rclr(bio, d$X[1:2, ], m[1:2, ])),
                   'INSUFFICIENT_BIOFLUID_ANCHORS')
  expect_identical(code(predict_tech_residualized_alr(tech, d$X[1:2, ], m[1:2, ])),
                   'INSUFFICIENT_ALR_PIVOTS')
})
