# Keep Python compilation and torch caches local to this test run. Configure
# these before any test imports torch, and restore the caller's environment.
local({
  cache_dir <- tempfile("omicselector-test-torch-")
  dir.create(cache_dir)
  # Torch precompiled headers use Python's default temp directory directly.
  cache_vars <- c("TORCHINDUCTOR_CACHE_DIR", "TRITON_CACHE_DIR", "TMPDIR")
  previous <- Sys.getenv(cache_vars, unset = NA_character_)
  names(previous) <- cache_vars
  do.call(Sys.setenv, as.list(stats::setNames(
    c(file.path(cache_dir, c("inductor", "triton")), cache_dir), cache_vars
  )))
  withr::defer({
    unlink(cache_dir, recursive = TRUE, force = TRUE)
    present <- !is.na(previous)
    if (any(present)) do.call(Sys.setenv, as.list(previous[present]))
    if (any(!present)) Sys.unsetenv(names(previous)[!present])
  }, envir = testthat::teardown_env())
})
