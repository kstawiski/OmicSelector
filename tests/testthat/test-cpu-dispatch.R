test_that("portable dispatch defaults preserve explicit caller settings", {
  keys <- c("ATEN_CPU_CAPABILITY", "MKL_ENABLE_INSTRUCTIONS",
    "ONEDNN_MAX_CPU_ISA", "DNNL_MAX_CPU_ISA", "OPENBLAS_CORETYPE")
  withr::local_envvar(setNames(rep(NA_character_, length(keys)), keys))
  skip_if_not(tolower(Sys.info()[["machine"]]) %in% c("x86_64", "amd64"))
  OmicSelector:::.ss_configure_cpu_dispatch(avx512 = FALSE)
  expect_identical(unname(Sys.getenv(keys, unset = NA_character_)),
                   c("avx2", "AVX2", "AVX2", "AVX2", NA_character_))
  OmicSelector:::.ss_configure_cpu_dispatch(avx512 = TRUE)
  expect_identical(Sys.getenv("OPENBLAS_CORETYPE"), "Zen")
  Sys.setenv(ATEN_CPU_CAPABILITY = "avx512", OPENBLAS_CORETYPE = "")
  OmicSelector:::.ss_configure_cpu_dispatch(avx512 = TRUE)
  expect_identical(Sys.getenv("ATEN_CPU_CAPABILITY"), "avx512")
  expect_identical(Sys.getenv("OPENBLAS_CORETYPE"), "")
})

test_that("capability and settings mismatches warn or fail before scoring", {
  fitted <- OmicSelector:::.ss_cpu_dispatch_provenance("AVX2")
  model <- list(provenance = list(cpu_dispatch = fitted))
  check <- OmicSelector:::.ss_check_cpu_dispatch
  expect_true(check(model, active = fitted))
  active <- fitted
  active$capability <- "AVX512"
  expect_warning(expect_false(check(model, active = active)), "capability=AVX512")
  expect_error(check(model, active = active, strict = TRUE), "CPU dispatch differs")
  active <- fitted
  active$settings[["MKL_ENABLE_INSTRUCTIONS"]] <- "AVX512"
  expect_warning(check(model, active = active), "MKL_ENABLE_INSTRUCTIONS=AVX512")
  withr::local_options(OmicSelector.cpu_dispatch_strict = TRUE)
  expect_error(check(model, active = active), "CPU dispatch differs")
  expect_error(check(list(), active = fitted), "unrecorded")
  expect_error(check(model, active = fitted, strict = NA), "strict must be")
  expect_identical(model$provenance$cpu_dispatch, fitted)
})
