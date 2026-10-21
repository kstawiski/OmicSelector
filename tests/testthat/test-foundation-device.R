.foundation_torch <- function() {
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("torch"), "Python torch unavailable")
  reticulate::import("torch", convert = FALSE)
}

.foundation_backend <- function(method, gpu = FALSE) {
  torch <- .foundation_torch()
  if (gpu) skip_if_not(reticulate::py_to_r(torch$cuda$is_available()), "CUDA unavailable")
  skip_if_not(reticulate::py_module_available(method), paste(method, "unavailable"))
  torch
}

.foundation_data <- function() {
  set.seed(902)
  X <- exp(matrix(rnorm(24 * 140), 24, 140,
                  dimnames = list(NULL, paste0("f", seq_len(140)))))
  list(X = X, y = rep(0:1, 12))
}

test_that("TabDPT CPU projection is deterministic and leaves RNG untouched", {
  torch <- .foundation_torch()
  mod <- .ss_foundation_module()
  X <- matrix(seq_len(24 * 140) %% 137, 24, 140)
  before <- torch$random$get_rng_state()$clone()
  a <- mod$cpu_projection(X, 128L, 17L)
  torch$manual_seed(819L)
  other <- reticulate::py_to_r(torch$random$get_rng_state()$numpy())
  b <- mod$cpu_projection(X, 128L, 17L)
  expect_identical(reticulate::py_to_r(a$numpy()), reticulate::py_to_r(b$numpy()))
  expect_identical(other, reticulate::py_to_r(torch$random$get_rng_state()$numpy()))
  torch$random$set_rng_state(before)
  expect_equal(dim(reticulate::py_to_r(a$numpy())), c(140L, 24L))
})

test_that("full precision restores the caller's kernel policy", {
  torch <- .foundation_torch()
  old <- reticulate::py_to_r(torch$backends$cuda$matmul$allow_tf32)
  ctx <- .ss_foundation_module()$full_precision()
  ctx$`__enter__`()
  expect_false(reticulate::py_to_r(torch$backends$cuda$matmul$allow_tf32))
  expect_false(reticulate::py_to_r(torch$backends$mha$get_fastpath_enabled()))
  ctx$`__exit__`(NULL, NULL, NULL)
  expect_identical(reticulate::py_to_r(torch$backends$cuda$matmul$allow_tf32), old)
})

for (.method in c("tabpfn", "tabdpt")) local({
  method <- .method
  test_that(paste(method, "records float32 and preserves frozen CPU state"), {
    .foundation_backend(method)
    d <- .foundation_data()
    model <- get(paste0("fit_", method))(d$X, d$y, hp = list(device = "cpu"))
    expect_identical(model$inference_settings$precision, "float32")
    expect_identical(model$inference_settings$attention, "math")
    if (method == "tabdpt") {
      expect_false(model$inference_settings$use_flash)
      expect_identical(model$inference_settings$projection_device, "cpu")
      expect_equal(nrow(model$projection), 140L)
    } else {
      expect_false(model$inference_settings$memory_saving_mode)
      expect_identical(model$inference_settings$memory_saving_policy,
        "training_rows_or_features_gt_512")
      expect_identical(model$inference_settings$encoder_device, "cpu")
      expect_identical(model$inference_settings$column_embedding_device, "cpu")
      expect_false(model$inference_settings$gpu_preprocessing)
      expect_identical(model$inference_settings$query_preprocessing, "singleton")
      expect_identical(model$inference_settings$fit_mode, "fit_preprocessors")
    }
    if (method == "tabpfn") {
      np <- reticulate::import("numpy", convert = FALSE)
      z <- np$array(.tabpfn_rclr_matrix(d$X[1:2, model$feature_universe, drop = FALSE]), dtype = "float64", copy = TRUE)
      member <- reticulate::py_get_item(model$clf$executor_$ensemble_members, 0L)
      batched <- member$transform_X_test(z)
      one <- member$transform_X_test(np$array(.tabpfn_rclr_matrix(d$X[1, model$feature_universe, drop = FALSE]), dtype = "float64", copy = TRUE))
      expect_identical(reticulate::py_to_r(batched)[1, ], as.numeric(reticulate::py_to_r(one)))
    }
    before <- singlesample_state_fingerprint(model)
    model$clf$to("cpu")
    expect_identical(singlesample_state_fingerprint(model), before)
    score <- get(paste0("score_", method))
    a <- score(model, d$X[1:2, , drop = FALSE])
    expect_identical(singlesample_state_fingerprint(model), before)
    expect_identical(a, vapply(1:2, function(i) score(model, d$X[i, , drop = FALSE]), numeric(1)))
  })
  test_that(paste(method, "CPU and GPU fits and round trip agree"), {
    .foundation_backend(method, gpu = TRUE)
    d <- .foundation_data()
    fit <- get(paste0("fit_", method)); score <- get(paste0("score_", method))
    cpu <- fit(d$X, d$y, hp = list(device = "cpu"))
    gpu <- fit(d$X, d$y, hp = list(device = "cuda"))
    X <- d$X[1:2, , drop = FALSE]
    a <- score(cpu, X); b <- score(gpu, X)
    expect_lt(max(abs(a - b)), 1e-5)
    expect_identical(a >= .5, b >= .5)
    if (method == "tabdpt") expect_identical(cpu$projection, gpu$projection)
    before <- singlesample_state_fingerprint(cpu)
    cpu$clf$to("cuda")
    moved <- score(cpu, X)
    expect_lt(max(abs(a - moved)), 1e-5)
    cpu$clf$to("cpu")
    expect_identical(singlesample_state_fingerprint(cpu), before)
  })
})

test_that("foundation copies isolate buffers while retaining parameter values", {
  torch <- .foundation_torch()
  probe <- reticulate::py_run_string("
import torch

def check_foundation_copy(module):
    net = torch.nn.Sequential(torch.nn.Linear(3, 3), torch.nn.BatchNorm1d(3)).eval()
    cloned = module.private_copy(net)
    assert torch.equal(net[0].weight, cloned[0].weight)
    assert net[0].weight.data_ptr() == cloned[0].weight.data_ptr()
    assert net[1].running_mean.data_ptr() != cloned[1].running_mean.data_ptr()
    cloned[1].running_mean.add_(7)
    assert torch.equal(net[1].running_mean, torch.zeros(3))
    assert net[0].weight.requires_grad
    assert not cloned[0].weight.requires_grad
    return True
", convert = FALSE)
  expect_true(reticulate::py_to_r(probe$check_foundation_copy(.ss_foundation_module())))
})
