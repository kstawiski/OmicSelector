"""Sinkhorn OT helpers for OmicSelector paper benchmarks.

The functions in this file are intentionally small and dependency-light.  POT
(`ot`) is used when it is available, but the internal NumPy Sinkhorn routine is
kept as the reproducible fallback used by the R wrapper tests.
"""

from __future__ import annotations

import math
from typing import Iterable, Mapping, Sequence

import numpy as np

try:  # pragma: no cover - depends on local optional environment
    import ot as _pot

    _POT_AVAILABLE = True
except Exception:  # pragma: no cover - exercised when POT is absent
    _pot = None
    _POT_AVAILABLE = False


def pot_available() -> bool:
    """Return whether the optional POT package was importable."""

    return bool(_POT_AVAILABLE)


def _as_cloud_list(cohort_clouds: Sequence[np.ndarray] | Mapping[str, np.ndarray]) -> list[np.ndarray]:
    if isinstance(cohort_clouds, Mapping):
        vals = list(cohort_clouds.values())
    else:
        vals = list(cohort_clouds)
    out: list[np.ndarray] = []
    for x in vals:
        arr = np.asarray(x, dtype=float)
        if arr.ndim != 2 or arr.shape[0] == 0 or arr.shape[1] == 0:
            continue
        keep = np.all(np.isfinite(arr), axis=1)
        if np.any(keep):
            out.append(arr[keep])
    if not out:
        raise ValueError("No non-empty finite cohort clouds supplied.")
    return out


def _normalise_weights(w: np.ndarray, n: int) -> np.ndarray:
    if w is None:
        out = np.ones(n, dtype=float) / float(n)
    else:
        out = np.asarray(w, dtype=float).reshape(-1)
        if out.size != n:
            raise ValueError(f"Weight length {out.size} does not match expected length {n}.")
        out[~np.isfinite(out) | (out < 0)] = 0.0
        s = float(out.sum())
        out = np.ones(n, dtype=float) / float(n) if s <= 0 else out / s
    return out


def _sqeuclidean_cost(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x2 = np.sum(x * x, axis=1, keepdims=True)
    y2 = np.sum(y * y, axis=1, keepdims=True).T
    c = x2 + y2 - 2.0 * np.dot(x, y.T)
    np.maximum(c, 0.0, out=c)
    finite = c[np.isfinite(c)]
    if finite.size == 0:
        raise ValueError("Non-finite Sinkhorn cost matrix.")
    med = float(np.median(finite[finite > 0])) if np.any(finite > 0) else 1.0
    if not math.isfinite(med) or med <= 0:
        med = 1.0
    return c / med


def sinkhorn_plan(
    source_weights: np.ndarray,
    target_weights: np.ndarray,
    cost: np.ndarray,
    epsilon: float = 0.1,
    max_iter: int = 250,
    tol: float = 1e-7,
    use_pot: bool = True,
) -> np.ndarray:
    """Compute an entropic transport plan.

    The cost is assumed to be finite and non-negative.  When POT is importable,
    this delegates to ``ot.sinkhorn``; otherwise it uses a standard scaling
    iteration with small floors to avoid division by zero.
    """

    a = _normalise_weights(source_weights, len(source_weights))
    b = _normalise_weights(target_weights, len(target_weights))
    c = np.asarray(cost, dtype=float)
    eps = max(float(epsilon), 1e-6)

    if use_pot and _POT_AVAILABLE:  # pragma: no cover - optional local path
        try:
            plan = _pot.sinkhorn(a, b, c, reg=eps, numItermax=int(max_iter), stopThr=float(tol))
            plan = np.asarray(plan, dtype=float)
            if np.all(np.isfinite(plan)) and plan.sum() > 0:
                return plan
        except Exception:
            pass

    kernel = np.exp(-c / eps)
    kernel[~np.isfinite(kernel)] = 0.0
    kernel = np.maximum(kernel, 1e-300)
    u = np.ones_like(a)
    v = np.ones_like(b)
    last_u = u.copy()
    for _ in range(int(max_iter)):
        kv = kernel @ v
        kv[kv <= 0] = 1e-300
        u = a / kv
        ktu = kernel.T @ u
        ktu[ktu <= 0] = 1e-300
        v = b / ktu
        if np.max(np.abs(u - last_u)) < tol:
            break
        last_u = u.copy()
    plan = (u[:, None] * kernel) * v[None, :]
    plan[~np.isfinite(plan)] = 0.0
    total = plan.sum()
    if total <= 0:
        return np.outer(a, b)
    return plan / total


def _initial_support(
    clouds: list[np.ndarray],
    n_atoms: int,
    rng: np.random.Generator,
) -> np.ndarray:
    pooled = np.vstack(clouds)
    n_atoms = int(min(max(1, n_atoms), pooled.shape[0]))
    replace = pooled.shape[0] < n_atoms
    idx = rng.choice(pooled.shape[0], size=n_atoms, replace=replace)
    support = pooled[idx].copy()
    if support.shape[0] < n_atoms:
        extra = rng.choice(pooled.shape[0], size=n_atoms - support.shape[0], replace=True)
        support = np.vstack([support, pooled[extra]])
    return support


def free_support_sinkhorn_barycenter(
    cohort_clouds: Sequence[np.ndarray] | Mapping[str, np.ndarray],
    cohort_weights: np.ndarray | None = None,
    n_atoms: int = 1000,
    epsilon: float = 0.1,
    reg_lambda: float = 0.5,
    barycenter_iter: int = 4,
    max_sinkhorn_iter: int = 250,
    tol: float = 1e-7,
    seed: int = 42,
    use_pot: bool = True,
) -> dict:
    """Approximate a free-support entropic W2 barycenter.

    Each iteration transports the current barycenter support to every training
    cohort cloud, row-normalizes each coupling, averages the mapped atom
    locations by cohort weights, and applies a damped update.  This is an
    intentionally conservative CPU implementation for reviewer benchmarking,
    not a replacement for specialized large-scale OT solvers.
    """

    clouds = _as_cloud_list(cohort_clouds)
    dim = clouds[0].shape[1]
    if any(x.shape[1] != dim for x in clouds):
        raise ValueError("All cohort clouds must have the same feature dimension.")
    rng = np.random.default_rng(int(seed))
    support = _initial_support(clouds, int(n_atoms), rng)
    atom_weights = np.ones(support.shape[0], dtype=float) / float(support.shape[0])
    cw = _normalise_weights(cohort_weights, len(clouds))
    damping = min(max(float(reg_lambda), 0.0), 1.0)
    if damping <= 0:
        damping = 0.5

    for _ in range(int(max(1, barycenter_iter))):
        mapped_sum = np.zeros_like(support)
        for weight, cloud in zip(cw, clouds):
            target_weights = np.ones(cloud.shape[0], dtype=float) / float(cloud.shape[0])
            cost = _sqeuclidean_cost(support, cloud)
            plan = sinkhorn_plan(
                atom_weights,
                target_weights,
                cost,
                epsilon=epsilon,
                max_iter=max_sinkhorn_iter,
                tol=tol,
                use_pot=use_pot,
            )
            row_mass = plan.sum(axis=1, keepdims=True)
            row_mass[row_mass <= 0] = 1.0
            mapped = (plan @ cloud) / row_mass
            mapped_sum += float(weight) * mapped
        support = (1.0 - damping) * support + damping * mapped_sum

    return {
        "atoms": support,
        "weights": atom_weights,
        "backend": "python_pot_sinkhorn" if (use_pot and _POT_AVAILABLE) else "python_internal_sinkhorn",
        "pot_available": bool(_POT_AVAILABLE),
        "n_atoms": int(support.shape[0]),
        "n_cohorts": int(len(clouds)),
    }


def project_to_barycenter(
    source_cloud: np.ndarray,
    barycenter_atoms: np.ndarray,
    barycenter_weights: np.ndarray | None = None,
    epsilon: float = 0.1,
    max_sinkhorn_iter: int = 250,
    tol: float = 1e-7,
    use_pot: bool = True,
    single_sample_kernel_fallback: bool = True,
) -> dict:
    """Project source samples onto a fitted barycenter by barycentric mapping."""

    x = np.asarray(source_cloud, dtype=float)
    atoms = np.asarray(barycenter_atoms, dtype=float)
    if x.ndim != 2 or atoms.ndim != 2:
        raise ValueError("source_cloud and barycenter_atoms must be 2D arrays.")
    if x.shape[1] != atoms.shape[1]:
        raise ValueError("Source and barycenter feature dimensions differ.")
    b = _normalise_weights(barycenter_weights, atoms.shape[0])
    cost = _sqeuclidean_cost(x, atoms)

    if x.shape[0] == 1 and single_sample_kernel_fallback:
        eps = max(float(epsilon), 1e-6)
        weights = b * np.exp(-cost[0] / eps)
        weights[~np.isfinite(weights)] = 0.0
        if weights.sum() <= 0:
            weights = b.copy()
        weights = weights / weights.sum()
        projection = weights.reshape(1, -1) @ atoms
        return {
            "projection": projection,
            "backend": "single_sample_entropic_kernel",
            "plan": weights.reshape(1, -1),
        }

    a = np.ones(x.shape[0], dtype=float) / float(x.shape[0])
    plan = sinkhorn_plan(
        a,
        b,
        cost,
        epsilon=epsilon,
        max_iter=max_sinkhorn_iter,
        tol=tol,
        use_pot=use_pot,
    )
    row_mass = plan.sum(axis=1, keepdims=True)
    row_mass[row_mass <= 0] = 1.0
    projection = (plan @ atoms) / row_mass
    return {
        "projection": projection,
        "backend": "python_pot_sinkhorn" if (use_pot and _POT_AVAILABLE) else "python_internal_sinkhorn",
        "plan": plan,
    }
