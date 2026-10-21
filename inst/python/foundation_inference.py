"""Device-neutral preprocessing and private, full-precision inference helpers."""
from contextlib import contextmanager
import copy
import types
import numpy as np


def cpu_projection(x, max_features, seed):
    """Fit randomized PCA with an isolated, explicitly seeded CPU generator.

    The QR/SVD construction is the same range-finder used by pca_lowrank.
    Neither the global generator nor a CUDA generator supplies random draws.
    """
    import torch
    x = torch.as_tensor(np.array(x, copy=True), device="cpu", dtype=torch.float32)
    q = min(x.shape[0], x.shape[1], int(max_features))
    centered = x - x.mean(dim=0, keepdim=True)
    generator = torch.Generator(device="cpu").manual_seed(int(seed))
    probe = torch.randn(x.shape[1], q, generator=generator, device="cpu", dtype=x.dtype)
    basis = torch.linalg.qr(centered @ probe, mode="reduced").Q
    for _ in range(2):
        basis = torch.linalg.qr(centered.T @ basis, mode="reduced").Q
        basis = torch.linalg.qr(centered @ basis, mode="reduced").Q
    _, _, vh = torch.linalg.svd(basis.T @ centered, full_matrices=False)
    return vh.T.contiguous()


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
    visit(estimator)
    return copy.deepcopy(estimator, memo)


def predict(estimator, x, **kwargs):
    with full_precision():
        return estimator.predict_proba(x, **kwargs)
