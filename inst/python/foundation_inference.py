"""Device-neutral preprocessing and private, full-precision inference helpers."""
from contextlib import contextmanager
import copy
import types
import numpy as np


def cpu_projection(x, max_features, seed, return_sign_margin=False,
                   return_metadata=False):
    """Fit centered PCA by exact float64 CPU SVD with canonical column signs.

    The seed argument is retained for callers; exact SVD consumes no random
    draws. Numerical rank uses the default NumPy matrix_rank tolerance.
    Return float32 directions and optionally their sign margin or diagnostics.
    """
    import torch
    x = torch.as_tensor(np.array(x, copy=True), device="cpu", dtype=torch.float64)
    if x.shape[0] < 2:
        raise ValueError("Centered PCA requires at least two training rows")
    _, singular, vh = torch.linalg.svd(
        x - x.mean(dim=0, keepdim=True), full_matrices=False)
    tolerance = float(singular[0]) * max(x.shape) * torch.finfo(torch.float64).eps
    rank = int((singular > tolerance).sum())
    q = min(rank, x.shape[1], int(max_features))
    projection = vh[:q].T.to(dtype=torch.float32)
    # Orient each direction by its largest absolute loading. argmax chooses
    # the first feature in an exact tie, independent of the initial SVD sign.
    magnitudes = projection.abs()
    anchors = magnitudes.argmax(dim=0)
    columns = torch.arange(projection.shape[1])
    signs = torch.where(projection[anchors, columns] < 0, -1.0, 1.0)
    projection = (projection * signs).contiguous()
    largest = magnitudes.topk(min(2, projection.shape[0]), dim=0).values
    margins = largest[0] - largest[1] if largest.shape[0] == 2 else largest[0]
    # No retained direction means there is no sign decision to resolve.
    margin = float(margins.min()) if q else float("inf")
    metadata = dict(rank=rank, q=q, sign_min_margin=margin,
                    rank_tolerance=tolerance,
                    min_retained_ratio=float(singular[q - 1]) / tolerance if q else None)
    if return_metadata:
        return projection, metadata
    return (projection, margin) if return_sign_margin else projection


@contextmanager
def full_precision():
    """Select float32 math kernels, restoring the caller's backend policy."""
    import torch
    from torch.nn.attention import SDPBackend, sdpa_kernel
    matmul = torch.backends.cuda.matmul.allow_tf32
    cudnn = torch.backends.cudnn.allow_tf32
    fastpath = torch.backends.mha.get_fastpath_enabled()
    try:
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        torch.backends.mha.set_fastpath_enabled(False)
        with sdpa_kernel(SDPBackend.MATH), torch.no_grad():
            yield
    finally:
        torch.backends.cuda.matmul.allow_tf32 = matmul
        torch.backends.cudnn.allow_tf32 = cudnn
        torch.backends.mha.set_fastpath_enabled(fastpath)


def private_copy(estimator):
    """Copy mutable inference state while sharing read-only neural parameters.

    Buffers, preprocessing, and query caches remain private. The supported
    foundation-model forward passes do not write parameters. In particular we
    never cache a second estimator that could become stale after a user edit.
    """
    import torch
    memo = {}
    seen = set()
    def visit(x):
        if id(x) in seen:
            return
        seen.add(id(x))
        if isinstance(x, torch.nn.Parameter):
            # A distinct Tensor wrapper permits freezing without changing the
            # fitted parameter's requires_grad flag; storage is read-only here.
            memo[id(x)] = torch.nn.Parameter(x.detach(), requires_grad=False)
        elif isinstance(x, dict):
            for v in x.values():
                visit(v)
        elif isinstance(x, (list, tuple)):
            for v in x:
                visit(v)
        elif hasattr(x, "__dict__") and not isinstance(x, (type, types.ModuleType, types.FunctionType)):
            for v in vars(x).values():
                visit(v)
    try:
        visit(estimator)
        return copy.deepcopy(estimator, memo)
    finally:
        # The recursive closure forms a cycle. deepcopy fills memo with the
        # entire estimator graph, including large context tensors. Release it
        # now rather than retaining original and private models until cyclic GC.
        memo.clear()
        seen.clear()


def predict(estimator, x, **kwargs):
    with full_precision():
        return estimator.predict_proba(x, **kwargs)
