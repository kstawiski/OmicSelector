test_that("deployed foundation calls release private Python contexts after each row", {
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("torch"), "Python torch unavailable")
  probe <- reticulate::py_run_string("
import numpy as np
import torch
import weakref

class DeployedMemoryProbe:
    buffers = []
    def __init__(self):
        self.context = torch.ones(4096)
    def __deepcopy__(self, memo):
        result = DeployedMemoryProbe()
        self.buffers.append(weakref.ref(result.context))
        return result
    def predict_proba(self, x, **kwargs):
        return np.tile([0.4, 0.6], (len(x), 1))
    @classmethod
    def retained_bytes(cls):
        return sum(r().numel() * r().element_size()
                   for r in cls.buffers if r() is not None)
", convert = FALSE)
  X <- matrix(c(1, 3, 2), 1, dimnames = list(NULL, c("a", "b", "c")))
  for (method in c("tabpfn", "tabdpt")) {
    model <- structure(list(clf = probe$DeployedMemoryProbe(),
      feature_universe = colnames(X),
      hp = list(seed = 1L, min_features = 1L, score_batch = TRUE, temperature = 1),
      case_col = 2L), class = paste0(method, "_model"))
    dep <- list(model = model,
      method_id = if (method == "tabpfn") "ai-tabpfn" else "tab-tabdpt")
    retained <- numeric(50)
    values <- numeric(50)
    for (i in seq_len(50)) {
      values[i] <- .ss_python_score(dep, X, NULL, strict = FALSE)
      retained[i] <- reticulate::py_to_r(probe$DeployedMemoryProbe$retained_bytes())
    }
    expect_equal(values, rep(0.6, 50))
    # Bound live Python storage, independent of the allocator's RSS policy.
    expect_identical(retained, numeric(50), info = method)
    expect_equal(.ss_python_score(dep, X[rep(1L, 50), , drop = FALSE],
      NULL, strict = FALSE), rep(0.6, 50))
    expect_equal(reticulate::py_to_r(probe$DeployedMemoryProbe$retained_bytes()), 0)
  }
})

test_that("direct calls and non-foundation deployments release copies on return and error", {
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("torch"), "Python torch unavailable")
  probe <- reticulate::py_run_string("
import numpy as np
import torch
import weakref

class ScoringMemoryProbe:
    copies = []
    fail = False
    promote = None
    def __init__(self):
        self.context = torch.ones(4096)
    def __deepcopy__(self, memo):
        result = ScoringMemoryProbe()
        self.copies.append(weakref.ref(result))
        return result
    def predict_proba(self, x, **kwargs):
        if self.promote is not None:
            type(self).promote()
        if self.fail:
            raise RuntimeError('deliberate prediction failure')
        return np.tile([0.4, 0.6], (len(x), 1))
    @classmethod
    def retained_bytes(cls):
        return sum(r().context.numel() * r().context.element_size()
                   for r in cls.copies if r() is not None)
", convert = FALSE)
  X <- matrix(c(1, 3, 2), 1, dimnames = list(NULL, c("a", "b", "c")))
  ids <- c(tabpfn = "ai-tabpfn", tabdpt = "tab-tabdpt", tabicl = "tab-tabicl")
  for (method in names(ids)) {
    for (route in c("direct", "deployed")) {
      gc()
      model <- structure(list(clf = probe$ScoringMemoryProbe(),
        feature_universe = colnames(X),
        hp = list(seed = 1L, min_features = 1L, score_batch = TRUE, temperature = 1),
        case_col = 2L), class = paste0(method, "_model"))
      dep <- list(model = model, method_id = ids[[method]])
      score <- function() if (route == "direct") {
        singlesample_score_call(ids[[method]], model, X)
      } else .ss_python_score(dep, X, NULL, strict = FALSE)
      values <- retained <- numeric(50)
      for (i in seq_len(50)) {
        values[i] <- score()
        retained[i] <- reticulate::py_to_r(probe$ScoringMemoryProbe$retained_bytes())
      }
      expect_equal(values, rep(0.6, 50))
      expect_identical(retained, numeric(50), info = paste(method, route))
      # A collection while prediction is active promotes its live R wrappers.
      # A young-generation collection alone then cannot release the copy.
      probe$ScoringMemoryProbe$promote <- function() invisible(gc())
      expect_equal(score(), 0.6)
      expect_equal(reticulate::py_to_r(probe$ScoringMemoryProbe$retained_bytes()), 0)
      probe$ScoringMemoryProbe$promote <- NULL
      probe$ScoringMemoryProbe$fail <- TRUE
      expect_error(score(), "deliberate prediction failure")
      expect_equal(reticulate::py_to_r(probe$ScoringMemoryProbe$retained_bytes()), 0,
                   info = paste(method, route, "error"))
      probe$ScoringMemoryProbe$fail <- FALSE
      expect_equal(score(), 0.6)
      expect_equal(reticulate::py_to_r(probe$ScoringMemoryProbe$retained_bytes()), 0)
    }
  }
})

test_that("generic Python deployments release cycles and partially created clones", {
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("torch"), "Python torch unavailable")
  probe <- reticulate::py_run_string("
import weakref
import numpy as np
class GenericMemoryProbe:
    copies = []
    def __init__(self):
        self.buffer = np.ones(4096)
        self.cycle = self
    def __deepcopy__(self, memo):
        result = GenericMemoryProbe()
        self.copies.append(weakref.ref(result))
        return result
    def score(self):
        return 0.6
    @classmethod
    def retained_bytes(cls):
        return sum(r().buffer.nbytes for r in cls.copies if r() is not None)
class FailingMemoryClone:
    def __deepcopy__(self, memo):
        result = GenericMemoryProbe()
        GenericMemoryProbe.copies.append(weakref.ref(result))
        raise RuntimeError('deliberate clone failure')
", convert = FALSE)
  testthat::local_mocked_bindings(singlesample_score_call = function(method_id, model, ...) {
    reticulate::py_to_r(model$first$score())
  }, .package = "OmicSelector")
  X <- matrix(1, 1, 1, dimnames = list(NULL, "a"))
  dep <- list(method_id = "memory-probe", model = list(first = probe$GenericMemoryProbe()))
  for (i in seq_len(20)) {
    expect_equal(.ss_python_score(dep, X, NULL, FALSE), 0.6)
    expect_equal(reticulate::py_to_r(probe$GenericMemoryProbe$retained_bytes()), 0)
  }
  dep$model$second <- probe$FailingMemoryClone()
  expect_error(.ss_python_score(dep, X, NULL, FALSE), "deliberate clone failure")
  expect_equal(reticulate::py_to_r(probe$GenericMemoryProbe$retained_bytes()), 0)
  dep$model <- list(second = probe$FailingMemoryClone())
  expect_error(.ss_python_score(dep, X, NULL, FALSE), "deliberate clone failure")
  expect_equal(reticulate::py_to_r(probe$GenericMemoryProbe$retained_bytes()), 0)
})
