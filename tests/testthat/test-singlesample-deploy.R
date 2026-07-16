make_singlesample_deploy_data <- function(n = 40L, seed = 101L) {
  set.seed(seed)
  sbp <- ws_default_sbp()
  sentinels <- c("ALL_NON_RBC", "ALL_NON_PLT", "TOP_K_BY_ABUNDANCE",
                 "TAIL_BY_ABUNDANCE")
  features <- unique(unlist(lapply(sbp, function(b) {
    c(b$numerator, b$denominator)
  }), use.names = FALSE))
  features <- features[!features %in% sentinels]
  features <- features[seq_len(30L)]

  y <- rep(c(0, 1), each = n / 2L)
  X <- matrix(stats::rgamma(n * length(features), shape = 20, rate = 2),
              nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  X[y == 1L, features[1:5]] <- X[y == 1L, features[1:5]] + 80
  X[y == 0L, features[6:10]] <- X[y == 0L, features[6:10]] + 80
  list(X = X, y = y)
}

expect_clean_deploy_outcome <- function(X, y, method) {
  res <- tryCatch(
    list(ok = TRUE,
         value = suppressWarnings(
           deploy_singlesample(X, y, method = method, verify = FALSE)
         )),
    error = function(e) list(ok = FALSE, message = conditionMessage(e))
  )

  if (isTRUE(res$ok)) {
    dep <- res$value
    expect_s3_class(dep, "singlesample_deployable")
    scores <- score_specimen(dep, X[seq_len(2L), , drop = FALSE])
    expect_type(scores, "double")
    expect_length(scores, 2L)
    expect_true(all(is.finite(scores)))
    expect_true(is_singlesample_deployable(dep, X_probe = X))
  } else {
    expect_match(
      res$message,
      "deploy_singlesample:.*(not directly deployable|requires|does not return one score|reticulate|tier R2)"
    )
    expect_false(grepl("unused argument", res$message, fixed = TRUE))
  }

  invisible(res)
}

classify_singlesample_deploy_smoke <- function(X, y) {
  roster <- singlesample_method_roster()
  methods <- roster$method_id
  failures <- character()
  status <- setNames(rep(NA_character_, length(methods)), methods)
  reason <- setNames(rep(NA_character_, length(methods)), methods)

  for (method in methods) {
    res <- tryCatch(
      list(ok = TRUE,
           value = suppressWarnings(
             deploy_singlesample(X, y, method = method, verify = FALSE)
           )),
      error = function(e) list(ok = FALSE, message = conditionMessage(e))
    )

    if (isTRUE(res$ok)) {
      dep <- res$value
      if (!inherits(dep, "singlesample_deployable")) {
        failures <- c(failures, paste(method, "returned a non-deployable object"))
        next
      }
      score_res <- tryCatch(score_specimen(dep, X[1L, ]),
                            error = function(e) e)
      if (inherits(score_res, "error")) {
        failures <- c(failures, paste(method, "score_specimen errored:",
                                      conditionMessage(score_res)))
        next
      }
      if (!is.numeric(score_res) || length(score_res) != 1L ||
          !is.finite(score_res)) {
        failures <- c(failures, paste(method, "score_specimen returned",
                                      paste(capture.output(str(score_res)),
                                            collapse = " ")))
        next
      }
      status[[method]] <- "working"
      reason[[method]] <- "working"
      next
    }

    msg <- res$message
    clear <- grepl(
      "not directly deployable|requires|tier R2|reticulate|does not return one score",
      msg
    )
    cryptic <- grepl(
      "unused argument|could not find function|array with dim",
      msg,
      fixed = FALSE
    )
    if (!clear || cryptic) {
      failures <- c(failures, paste(method, "cryptic error:", msg))
      next
    }
    status[[method]] <- "rejected"
    if (grepl("tier R2|reticulate", msg)) {
      reason[[method]] <- "tier R2/reticulate"
    } else if (grepl("does not return one score", msg)) {
      reason[[method]] <- "non-scalar scorer"
    } else if (grepl("not directly deployable", msg)) {
      reason[[method]] <- "not directly deployable"
    } else {
      reason[[method]] <- "requires backend/function/signature"
    }
  }

  list(status = status, reason = reason, failures = failures)
}

test_that("deploy_singlesample fits default ws-balance-ilr deployment object", {
  dat <- make_singlesample_deploy_data()

  dep <- deploy_singlesample(dat$X, dat$y)

  expect_s3_class(dep, "singlesample_deployable")
  expect_s3_class(dep$model, "ws_balance_ilr_model")
  expect_equal(dep$method_id, "ws-balance-ilr")
  expect_true(is.function(dep$score_fn))
  expect_equal(dep$n_train, nrow(dat$X))
  expect_true(is.numeric(dep$fit_time))
  expect_true(isTRUE(dep$verified))
  expect_equal(dep$verify_n, nrow(dat$X))
})

test_that("score_specimen gives the same score for singleton and batch rows", {
  dat <- make_singlesample_deploy_data()
  dep <- deploy_singlesample(dat$X, dat$y)

  single <- score_specimen(dep, dat$X[1L, ])
  batch <- score_ws_balance_ilr(dep$model, dat$X)[1L]
  delta <- max(abs(single - batch))

  expect_lt(delta, 1e-9)
  expect_equal(single, batch, tolerance = 1e-12)
})

test_that("method selection errors are explicit", {
  dat <- make_singlesample_deploy_data()

  expect_error(
    deploy_singlesample(dat$X, dat$y, method = "not-a-method",
                        verify = FALSE),
    "Available methods:.*ws-balance-ilr"
  )

  roster <- singlesample_method_roster()
  torch_method <- roster$method_id[grepl("torch", roster$dep_route,
                                         ignore.case = TRUE)][1L]
  expect_false(is.na(torch_method))
  expect_error(
    deploy_singlesample(dat$X, dat$y, method = torch_method,
                        verify = FALSE),
    "tier R2.*Python backend / venv.*out of scope"
  )
})

test_that("tier R2 reticulate methods are rejected before Python fitting", {
  dat <- make_singlesample_deploy_data()

  for (method in c("ai-tabpfn", "tab-tabicl")) {
    expect_error(
      deploy_singlesample(dat$X, dat$y, method = method, verify = FALSE),
      "tier R2.*Reticulate methods.*Python backend / venv"
    )
  }
})

test_that("non-canonical fit signatures do not leak unused-argument errors", {
  dat <- make_singlesample_deploy_data()

  invisible(expect_clean_deploy_outcome(dat$X, dat$y, "frozen-ruv"))
  invisible(expect_clean_deploy_outcome(dat$X, dat$y, "ws-dominance"))
})

test_that("trimmed-rCLR reference deploys through its frozen signed panel", {
  dat <- make_singlesample_deploy_data()

  deployed <- deploy_singlesample(
    dat$X, dat$y, method = "ws-rclr-trimmed", verify = FALSE)
  expect_s3_class(deployed, "singlesample_deployable")
  expect_equal(length(score_specimen(deployed, dat$X[1:3, , drop = FALSE])), 3L)
})

test_that("verify gate and user-facing deployability check pass", {
  dat <- make_singlesample_deploy_data()
  dep <- deploy_singlesample(dat$X, dat$y, verify = TRUE)

  expect_true(isTRUE(dep$verified))
  expect_true(is_singlesample_deployable(dep, X_probe = dat$X))
})

test_that("deployment input validation fails explicitly", {
  dat <- make_singlesample_deploy_data()
  dep <- deploy_singlesample(dat$X, dat$y)

  expect_error(
    is_singlesample_deployable(dep, X_probe = NULL),
    "X_probe is required"
  )
  expect_error(
    score_specimen(dep, as.numeric(dat$X[1L, ])),
    "feature names"
  )
  expect_error(
    score_specimen(dep, dat$X[1L:2L, , drop = FALSE],
                   meta = data.frame(batch = "A")),
    "meta must have one row per row of x_new"
  )
  expect_error(
    deploy_singlesample(dat$X, dat$y, meta_train = data.frame(batch = "A"),
                        verify = FALSE),
    "meta_train must have one row per row of X_train"
  )
})

test_that("score_specimen does not mutate the deployment object", {
  dat <- make_singlesample_deploy_data()
  dep <- deploy_singlesample(dat$X, dat$y)
  before <- serialize(dep, connection = NULL)

  invisible(score_specimen(dep, dat$X[2L, ]))

  expect_identical(serialize(dep, connection = NULL), before)
})

test_that("score_specimen accepts vector, matrix, and data-frame inputs", {
  dat <- make_singlesample_deploy_data()
  dep <- deploy_singlesample(dat$X, dat$y)

  s_vec <- score_specimen(dep, dat$X[3L, ])
  s_mat <- score_specimen(dep, dat$X[3L, , drop = FALSE])
  s_df <- score_specimen(dep, as.data.frame(dat$X[3L:4L, , drop = FALSE]))

  expect_equal(s_vec, s_mat, tolerance = 0)
  expect_equal(s_df, score_ws_balance_ilr(dep$model, dat$X[3L:4L, ]),
               tolerance = 1e-12)
})

test_that("score_specimen keeps NA handling with the underlying scorer", {
  dat <- make_singlesample_deploy_data()
  dep <- deploy_singlesample(dat$X, dat$y)
  x_bad <- dat$X[1L, , drop = FALSE]
  x_bad[1L, 1L] <- NA_real_

  expect_error(score_specimen(dep, x_bad), "finite numeric values")
  expect_error(score_ws_balance_ilr(dep$model, x_bad),
               "finite numeric values")
})

test_that("print method reports deployment framing without superiority claims", {
  dat <- make_singlesample_deploy_data()
  dep <- deploy_singlesample(dat$X, dat$y)

  expect_output(print(dep), "method_id: ws-balance-ilr")
  expect_output(print(dep), "deployability only")
  expect_output(print(dep), "batch-corrected pipeline")
})

test_that("all roster methods deploy or fail with deploy-specific errors", {
  dat <- make_singlesample_deploy_data(n = 20L)

  smoke <- classify_singlesample_deploy_smoke(dat$X, dat$y)
  counts <- table(smoke$status, useNA = "no")
  reasons <- table(smoke$reason[smoke$status == "rejected"], useNA = "no")
  message(
    "singlesample deploy smoke counts: ",
    paste(names(counts), as.integer(counts), sep = "=", collapse = ", "),
    "; rejected reasons: ",
    paste(names(reasons), as.integer(reasons), sep = "=", collapse = ", ")
  )

  expect_equal(smoke$failures, character())
  expect_equal(sum(counts), nrow(singlesample_method_roster()))
})

test_that("deploy guard invariant: tier=='R2' is exactly the reticulate dep_route set", {
  # deploy_singlesample() rejects roster tier 'R2' as the reticulate/Python-backend
  # set; this invariant locks the docstring + guard against future manifest drift
  # (an R2 non-reticulate row, or an R1 reticulate row, would silently break both).
  m <- singlesample_method_roster()
  expect_true(all((m$tier == "R2") == grepl("reticulate", m$dep_route)),
              info = "tier R2 must equal the reticulate dep_route set for the deploy guard to be correct")
})
