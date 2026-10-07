test_that("StableMate fallback uses the first candidate after active libraries without changing libPaths", {
  fn <- .ss_load_installed_namespace
  env <- new.env(parent = environment(fn))
  environment(fn) <- env
  loaded <- character()
  calls <- 0L
  env$isNamespaceLoaded <- function(package) FALSE
  env$.ss_login_home <- function() "/fixture"
  env$list.dirs <- function(...) c("/fixture/R/x86_64-pc-linux-gnu-library/4.5",
                                  "/fixture/R/x86_64-pc-linux-gnu-library/4.6")
  env$dir.exists <- function(path) rep(TRUE, length(path))
  env$loadNamespace <- function(package, lib.loc) {
    expect_identical(package, "StableMate")
    loaded <<- lib.loc
    calls <<- calls + 1L
  }
  before <- .libPaths()
  expect_true(fn("StableMate"))
  expect_identical(loaded, unique(c(before, "/fixture/R/omicselector")))
  expect_identical(calls, 1L)
  expect_identical(.libPaths(), before)
  env$dir.exists <- function(path) grepl("/4.5/StableMate$", path)
  expect_true(fn("StableMate"))
  expect_identical(loaded, unique(c(before, "/fixture/R/x86_64-pc-linux-gnu-library/4.5")))
  env$dir.exists <- function(path) rep(TRUE, length(path))
  env$loadNamespace <- function(...) { calls <<- calls + 1L; stop("incompatible namespace") }
  calls <- 0L
  expect_false(fn("StableMate"))
  expect_identical(calls, 1L)
  expect_identical(.libPaths(), before)
  env$dir.exists <- function(path) rep(FALSE, length(path))
  expect_false(fn("StableMate"))
  expect_identical(calls, 1L)
  env$isNamespaceLoaded <- function(package) TRUE
  env$.ss_login_home <- function() stop("already loaded: lookup must not run")
  expect_true(fn("StableMate"))
})

test_that("StableMate fitter tries the fallback only after lookup fails and stops if still unavailable", {
  fn <- fit_sel_stablemate
  env <- new.env(parent = environment(fn))
  environment(fn) <- env
  available <- FALSE
  calls <- 0L
  succeeds <- FALSE
  env$requireNamespace <- function(...) available
  env$.ss_load_installed_namespace <- function(package) {
    expect_identical(package, "StableMate")
    calls <<- calls + 1L
    available <<- succeeds
  }
  env$.reo_check_matrix <- function(...) stop("dependency passed; matrix validation reached")
  expect_error(fn(NULL, NULL), "package 'StableMate' is required", fixed = TRUE)
  expect_identical(calls, 1L)
  succeeds <- TRUE
  expect_error(fn(NULL, NULL), "matrix validation reached", fixed = TRUE)
  expect_identical(calls, 2L)
  expect_error(fn(NULL, NULL), "matrix validation reached", fixed = TRUE)
  expect_identical(calls, 2L)
})
