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
  expect_equal(dim(reticulate::py_to_r(a$numpy())), c(140L, 23L))
})

test_that("centered projection fixes rank and signs across CPU thread counts", {
  .foundation_torch()
  probe <- reticulate::py_run_string("
import numpy as np
import torch
from unittest.mock import patch

def check_projection_rank(module):
    previous_threads = torch.get_num_threads()
    rng = torch.random.get_rng_state().clone()
    x = np.random.default_rng(912).normal(size=(24, 140)).astype(np.float32)
    try:
        projections = []
        for threads in (1, 2, 4, 8):
            torch.set_num_threads(threads)
            with patch.object(torch.linalg, 'svd', wraps=torch.linalg.svd) as svd:
                v, margin = module.cpu_projection(x, 128, 17, return_sign_margin=True)
                assert svd.call_args.kwargs == dict(full_matrices=False)
                assert svd.call_args.args[0].dtype == torch.float64
                torch.testing.assert_close(svd.call_args.args[0].mean(dim=0), torch.zeros(140, dtype=torch.float64), rtol=0, atol=1e-15)
            assert tuple(v.shape) == (140, 23)
            anchors = v.abs().argmax(dim=0)
            assert (v[anchors, torch.arange(v.shape[1])] >= 0).all()
            top = v.abs().topk(2, dim=0).values
            assert margin == float((top[0] - top[1]).min())
            assert margin > 0
            assert torch.equal(rng, torch.random.get_rng_state())
            projections.append(v)
        for v in projections[1:]:
            torch.testing.assert_close(v, projections[0], rtol=0, atol=1e-5)
        for n, p, q in ((129, 140, 128), (150, 140, 128), (24, 10, 10)):
            z = np.random.default_rng(81).normal(size=(n, p)).astype(np.float32)
            v = module.cpu_projection(z, 128, 17)
            centered = torch.as_tensor(z, dtype=torch.float64)
            centered = centered - centered.mean(dim=0, keepdim=True)
            old = torch.linalg.svd(centered, full_matrices=False)[2][:q].T.float()
            assert tuple(v.shape) == (p, q)
            anchors = old.abs().argmax(dim=0)
            signs = torch.where(old[anchors, torch.arange(old.shape[1])] < 0, -1.0, 1.0)
            assert torch.equal(v, old * signs)
            assert torch.equal(v, module.cpu_projection(z, 128, 999))
        # A duplicate training row removes another determined direction.
        duplicate = x.copy()
        duplicate[-1] = duplicate[0]
        reference = None
        for threads in (1, 2, 4, 8):
            torch.set_num_threads(threads)
            v, info = module.cpu_projection(duplicate, 128, 17, return_metadata=True)
            assert info['rank'] == info['q'] == 22
            assert v.shape == (140, 22) and info['min_retained_ratio'] > 1
            assert info['sign_min_margin'] > 0
            if reference is not None:
                torch.testing.assert_close(v, reference, rtol=0, atol=1e-5)
            reference = v
        # Constant data have rank zero; no arbitrary direction is retained.
        v, info = module.cpu_projection(np.ones((4, 140)), 128, 17, return_metadata=True)
        assert v.shape == (140, 0) and info['rank'] == info['q'] == 0
        assert info['min_retained_ratio'] is None
        # An exactly rank-two matrix with dependent rows and columns.
        z = np.arange(12, dtype=float).reshape(4, 3)
        z = np.column_stack((z, np.array([0., 1., 0., 1.])))
        v, info = module.cpu_projection(z, 128, 17, return_metadata=True)
        assert info['rank'] == info['q'] == np.linalg.matrix_rank(z - z.mean(axis=0)) == 2
        original_svd = torch.linalg.svd
        def reversed_svd(*args, **kwargs):
            u, singular, vh = original_svd(*args, **kwargs)
            return -u, singular, -vh
        baseline = module.cpu_projection(x, 128, 17)
        with patch.object(torch.linalg, 'svd', side_effect=reversed_svd):
            assert torch.equal(baseline, module.cpu_projection(x, 128, 17))
        # A loading tie uses the first feature; the reported margin is zero.
        def tied_svd(*args, **kwargs):
            return None, torch.tensor([1.0], dtype=torch.float64), torch.tensor([[-0.5, 0.5]], dtype=torch.float64)
        with patch.object(torch.linalg, 'svd', side_effect=tied_svd):
            v, margin = module.cpu_projection(np.ones((2, 2)), 1, 17, return_sign_margin=True)
            assert torch.equal(v[:, 0], torch.tensor([0.5, -0.5])) and margin == 0
        for n in (0, 1):
            try:
                module.cpu_projection(np.zeros((n, 140)), 128, 17)
            except ValueError as error:
                assert 'at least two training rows' in str(error)
            else:
                raise AssertionError('insufficient training rows accepted')
        return True
    finally:
        torch.set_num_threads(previous_threads)
        torch.random.set_rng_state(rng)
", convert = FALSE)
  expect_true(reticulate::py_to_r(probe$check_projection_rank(.ss_foundation_module())))
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
      expect_identical(model$inference_settings$fit_mode, "fit_with_cache")
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

test_that("copying and fingerprinting release estimator storage without cyclic GC", {
  .foundation_torch()
  probe <- reticulate::py_run_string("
import gc
import weakref
import torch

def check_immediate_release(foundation, state):
    enabled = gc.isenabled()
    gc.collect()
    gc.disable()
    try:
        net = torch.nn.Linear(3, 3)
        net.register_buffer('context', torch.ones(1024))
        before = state.state_fingerprint(net)
        cloned = foundation.private_copy(net)
        cloned_ref = weakref.ref(cloned)
        buffer_ref = weakref.ref(cloned.context)
        del cloned
        assert cloned_ref() is None, 'private estimator retained after scoring'
        assert buffer_ref() is None, 'private context storage retained after scoring'
        assert state.state_fingerprint(net) == before
        net_ref = weakref.ref(net)
        context_ref = weakref.ref(net.context)
        del net
        assert net_ref() is None, 'fingerprint retained the fitted estimator'
        assert context_ref() is None, 'fingerprint retained context storage'

        # Fail-closed hashing must release objects on the error path too.
        tensor = torch.ones(1024)
        tensor_ref = weakref.ref(tensor)
        invalid = {'context': tensor, 'unsupported': iter([])}
        try:
            state.state_fingerprint(invalid)
        except TypeError:
            pass
        else:
            raise AssertionError('unsupported state was accepted')
        del tensor, invalid
        assert tensor_ref() is None, 'failed fingerprint retained context storage'
        return True
    finally:
        if enabled:
            gc.enable()
", convert = FALSE)
  expect_true(reticulate::py_to_r(probe$check_immediate_release(
    .ss_foundation_module(), .ss_python_state_module())))
})


test_that("TabPFN cached and uncached frozen contexts give equal scores", {
  .foundation_backend("tabpfn")
  d <- .foundation_data()
  cached <- fit_tabpfn(d$X, d$y, hp = list(device = "cpu"))
  expect_identical(cached$inference_settings$context_cache_device, "cpu")
  expect_identical(cached$inference_settings$context_cache_precision, "float32")
  frozen <- reticulate::import_from_path("tabpfn_frozen",
    path = system.file("python", package = "OmicSelector"), convert = FALSE)
  probe <- reticulate::py_run_string("
from unittest.mock import patch
def without_context_cache(module):
    return patch.object(module, 'cache_context', lambda clf: clf)
def change_cached_tensor(clf):
    for module in clf.executor_.models[0].modules():
        for name, value in module.named_buffers(recurse=False):
            if name in ('_k_cache', '_v_cache', '_kv_cache') and value is not None:
                import torch
                with torch.inference_mode():
                    value.reshape(-1)[0] += 1
                return True
    return False
", convert = FALSE)
  ctx <- probe$without_context_cache(frozen)
  ctx$`__enter__`()
  uncached <- tryCatch(fit_tabpfn(d$X, d$y, hp = list(device = "cpu")),
    finally = ctx$`__exit__`(NULL, NULL, NULL))
  q <- d$X[1:3, , drop = FALSE]
  a <- score_tabpfn(uncached, q)
  b <- score_tabpfn(cached, q)
  expect_lt(max(abs(a - b)), 1e-5)
  expect_identical(a >= .5, b >= .5)
  before <- singlesample_state_fingerprint(cached)
  cached$clf$to("cpu")
  expect_identical(singlesample_state_fingerprint(cached), before)
  cached$hp$score_batch <- TRUE
  expect_equal(score_tabpfn(cached, q), b, tolerance = 1e-5)
  cached$hp$score_batch <- FALSE
  expect_identical(singlesample_state_fingerprint(cached), before)
  expect_true(reticulate::py_to_r(probe$change_cached_tensor(cached$clf)))
  expect_false(identical(singlesample_state_fingerprint(cached), before))
})
