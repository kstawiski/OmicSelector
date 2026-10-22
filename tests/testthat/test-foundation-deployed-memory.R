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
