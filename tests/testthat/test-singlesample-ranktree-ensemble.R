ranktree_example <- function() {
  set.seed(812)
  X <- matrix(stats::rgamma(72 * 10, shape = 2), nrow = 72,
              dimnames = list(paste0("s", seq_len(72)),
                              c("miR-1", "miR.2", "miR_3",
                                paste0("miR-", 4:10))))
  y <- rep(c(0, 1), each = 36)
  X[y == 1, 1:3] <- X[y == 1, 1:3] + 5
  list(X = X, y = y)
}

test_that("external rank-tree registry cannot expand the frozen primary roster", {
  primary <- singlesample_method_roster()
  external <- singlesample_external_competitor_roster()
  expect_equal(nrow(primary), 74L)
  expect_equal(nrow(external), 2L)
  expect_setequal(external$method_id, c("reo-rankforest", "reo-rankboost"))
  expect_true(all(!external$primary_family))
  expect_length(intersect(primary$method_id, external$method_id), 0L)
})

test_that("rank-tree hyperparameters freeze deployable datrank=FALSE", {
  rf <- OmicSelector:::.reo_rankforest_resolve_hp(list(ntree = 25L))
  boost <- OmicSelector:::.reo_rankboost_resolve_hp(list(ntree = 30L))
  expect_false(rf$datrank)
  expect_false(boost$datrank)
  expect_equal(rf$rf_cores, 1L)
  expect_equal(rf$seed, 1L)
  expect_equal(rf$max_screen_features, 128L)
  expect_equal(boost$nodedepth, 3L)
  expect_equal(boost$max_screen_features, 128L)
  expect_equal(
    OmicSelector:::.reo_ranktree_dimreduce_quantile(2548L, 128L, TRUE),
    1 - 128 / 2548
  )
  expect_equal(
    OmicSelector:::.reo_ranktree_dimreduce_quantile(100L, 128L, TRUE),
    0.75
  )
  expect_false(
    OmicSelector:::.reo_ranktree_dimreduce_quantile(2548L, 128L, FALSE)
  )
  expect_error(
    OmicSelector:::.reo_rankforest_resolve_hp(list(datarank = TRUE)),
    "unknown hp field"
  )
  expect_error(
    OmicSelector:::.reo_rankboost_resolve_hp(list(weights = 1)),
    "unknown hp field"
  )
  expect_error(
    OmicSelector:::.reo_rankboost_resolve_hp(list(max_screen_features = 0)),
    "positive integer"
  )
})

test_that("rank-tree probability extraction preserves singleton array shape", {
  value <- array(
    c(0.8, 0.2), dim = c(1L, 2L, 1L),
    dimnames = list("sample", c("control", "case"), "20")
  )
  expect_equal(
    OmicSelector:::.reo_ranktree_extract_case_probability(
      value, 1L, "case", "test"
    ),
    0.2
  )
})

test_that("rank-tree backend pair contract recovers the training feature order", {
  backend <- structure(list(
    var.names = c(
      "feature_0000002_less_than_feature_0000001",
      "feature_0000002_less_than_feature_0000003",
      "feature_0000001_less_than_feature_0000003"
    )
  ), class = "gbm")
  expect_identical(
    OmicSelector:::.reo_ranktree_backend_feature_order(
      backend, sprintf("feature_%07d", 1:3), "test"
    ),
    c("feature_0000002", "feature_0000001", "feature_0000003")
  )
  broken <- backend
  broken$var.names <- broken$var.names[-1]
  expect_error(
    OmicSelector:::.reo_ranktree_backend_feature_order(
      broken, sprintf("feature_%07d", 1:3), "test"
    ),
    "incomplete or inconsistent"
  )
})

test_that("rank-tree wrappers fail closed without their exact audited backend", {
  if (requireNamespace("ranktreeEnsemble", quietly = TRUE) &&
      identical(as.character(utils::packageVersion("ranktreeEnsemble")), "0.24")) {
    skip("exact audited backend is installed")
  }
  dat <- ranktree_example()
  expect_error(fit_reo_rankforest(dat$X, dat$y),
               "ranktreeEnsemble 0.24|requires ranktreeEnsemble 0.24")
})

test_that("random rank forest is deterministic and singleton-equivariant", {
  skip_if_not_installed("ranktreeEnsemble", minimum_version = "0.24")
  skip_if_not(identical(as.character(utils::packageVersion("ranktreeEnsemble")),
                        "0.24"))
  dat <- ranktree_example()
  set.seed(99)
  before <- get(".Random.seed", envir = globalenv())
  old_rf_cores <- getOption("rf.cores")
  options(rf.cores = 3L)
  on.exit(options(rf.cores = old_rf_cores), add = TRUE)
  model <- fit_reo_rankforest(
    dat$X, dat$y, hp = list(ntree = 20L, seed = 7L)
  )
  expect_true(length(model$backend_safe_feature_names) < ncol(dat$X))
  expected_pairs <- unlist(lapply(
    seq_len(length(model$backend_safe_feature_names) - 1L),
    function(i) {
      paste0(
        model$backend_safe_feature_names[[i]], "_less_than_",
        model$backend_safe_feature_names[
          (i + 1L):length(model$backend_safe_feature_names)
        ]
      )
    }
  ), use.names = FALSE)
  expect_setequal(
    model$backend$xvar.names,
    expected_pairs
  )
  model_bytes <- serialize(model, NULL, version = 3L)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  expect_equal(getOption("rf.cores"), 3L)
  batch <- score_reo_rankforest(model, dat$X[1:12, , drop = FALSE])
  single <- vapply(seq_len(12), function(i) {
    score_reo_rankforest(model, dat$X[i, , drop = FALSE])
  }, numeric(1))
  expect_equal(batch, single, tolerance = 1e-12)
  perm <- c(8:12, 1:7)
  expect_equal(score_reo_rankforest(model, dat$X[perm, , drop = FALSE]),
               batch[perm], tolerance = 1e-12)
  duplicated <- score_reo_rankforest(
    model, rbind(dat$X[1:4, , drop = FALSE], dat$X[2, , drop = FALSE])
  )
  expect_equal(duplicated[5], batch[2], tolerance = 1e-12)
  selected_original <- model$feature_names[
    match(model$backend_safe_feature_names, model$safe_feature_names)
  ]
  expect_equal(
    score_reo_rankforest(model, dat$X[1:12, selected_original, drop = FALSE]),
    batch, tolerance = 1e-12
  )
  model_again <- fit_reo_rankforest(
    dat$X, dat$y, hp = list(ntree = 20L, seed = 7L)
  )
  expect_equal(score_reo_rankforest(model_again, dat$X[1:12, , drop = FALSE]),
               batch, tolerance = 1e-12)
  expect_error(
    score_reo_rankforest(
      model, dat$X[1:3, setdiff(colnames(dat$X), selected_original[[1L]]), drop = FALSE]
    ),
    "lacks frozen backend feature"
  )
  expect_identical(serialize(model, NULL, version = 3L), model_bytes)
})

test_that("boosted rank trees are deterministic and singleton-equivariant", {
  skip_if_not_installed("ranktreeEnsemble", minimum_version = "0.24")
  skip_if_not(identical(as.character(utils::packageVersion("ranktreeEnsemble")),
                        "0.24"))
  dat <- ranktree_example()
  set.seed(101)
  before <- get(".Random.seed", envir = globalenv())
  model <- fit_reo_rankboost(
    dat$X, dat$y,
    hp = list(ntree = 20L, seed = 11L, nodedepth = 2L, nodesize = 3L)
  )
  expect_true(length(model$backend_safe_feature_names) < ncol(dat$X))
  model_bytes <- serialize(model, NULL, version = 3L)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  batch <- score_reo_rankboost(model, dat$X[1:12, , drop = FALSE])
  single <- vapply(seq_len(12), function(i) {
    score_reo_rankboost(model, dat$X[i, , drop = FALSE])
  }, numeric(1))
  expect_equal(batch, single, tolerance = 1e-12)
  perm <- c(8:12, 1:7)
  expect_equal(score_reo_rankboost(model, dat$X[perm, , drop = FALSE]),
               batch[perm], tolerance = 1e-12)
  duplicated <- score_reo_rankboost(
    model, rbind(dat$X[1:4, , drop = FALSE], dat$X[2, , drop = FALSE])
  )
  expect_equal(duplicated[5], batch[2], tolerance = 1e-12)
  selected_original <- model$feature_names[
    match(model$backend_safe_feature_names, model$safe_feature_names)
  ]
  expect_equal(
    score_reo_rankboost(model, dat$X[1:12, selected_original, drop = FALSE]),
    batch, tolerance = 1e-12
  )
  expect_error(
    score_reo_rankboost(
      model, dat$X[1:3, setdiff(colnames(dat$X), selected_original[[1L]]), drop = FALSE]
    ),
    "lacks frozen backend feature"
  )
  model_again <- fit_reo_rankboost(
    dat$X, dat$y,
    hp = list(ntree = 20L, seed = 11L, nodedepth = 2L, nodesize = 3L)
  )
  expect_equal(score_reo_rankboost(model_again, dat$X[1:12, , drop = FALSE]),
               batch, tolerance = 1e-12)
  expect_identical(serialize(model, NULL, version = 3L), model_bytes)
})
