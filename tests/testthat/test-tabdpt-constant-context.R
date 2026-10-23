# A constant training coordinate may become detected in a held-out specimen.
# Run on CPU with a small context; never download a checkpoint in CI.
test_that("TabDPT scores a query detected in a training-constant coordinate", {
  skip_if_not_installed("reticulate")
  old_offline <- Sys.getenv("HF_HUB_OFFLINE", unset = NA_character_)
  Sys.setenv(HF_HUB_OFFLINE = "1")
  on.exit(if (is.na(old_offline)) Sys.unsetenv("HF_HUB_OFFLINE") else
    Sys.setenv(HF_HUB_OFFLINE = old_offline), add = TRUE)
  skip_if_not(reticulate::py_module_available("tabdpt"), "Python TabDPT unavailable")
  hub <- reticulate::import("huggingface_hub", convert = TRUE)
  checkpoint <- hub$try_to_load_from_cache("Layer6/TabDPT", "tabdpt1_2.safetensors")
  skip_if_not(is.character(checkpoint) && length(checkpoint) == 1L &&
    file.exists(checkpoint), "Offline TabDPT checkpoint unavailable")

  set.seed(809L)
  X <- exp(matrix(rnorm(12L * 140L), 12L, 140L,
    dimnames = list(NULL, paste0("f", seq_len(140L)))))
  X[, 135:140] <- 0
  y <- rep(0:1, 6L)
  query <- X[1:2, , drop = FALSE]
  query[, 135:140] <- exp(seq_len(12L) / 4)
  model <- fit_tabdpt(X, y, hp = list(device = "cpu"))
  expect_true(all(model$context_rclr[, 135:140] == 0))
  frozen <- singlesample_state_fingerprint(model)
  scores <- score_tabdpt(model, query)
  expect_true(all(is.finite(scores)))
  expect_true(all(scores >= 0 & scores <= 1))
  expect_identical(scores, vapply(1:2, function(i)
    score_tabdpt(model, query[i, , drop = FALSE]), numeric(1)))
  expect_identical(score_tabdpt(model, query[2:1, , drop = FALSE]), rev(scores))
  expect_identical(singlesample_state_fingerprint(model), frozen)
})
