"""Deterministic state hashing and read-only inference for fitted Python objects."""
import copy
import hashlib
import inspect
import random
import types
from pathlib import PurePath
from contextlib import contextmanager
import numpy as np


def state_fingerprint(obj):
    # Traverse estimator state, including nested torch parameters/buffers,
    # preprocessing and training context. Never hash pointer addresses or pickle
    # storage IDs. Repeated references and cycles retain deterministic indices.
    import torch
    h = hashlib.sha256()
    seen = {}
    keepalive = []  # Prevent temporary state_dict object IDs from being reused.
    unsupported = set()

    def put(value):
        data = value if isinstance(value, bytes) else str(value).encode('utf-8')
        h.update(len(data).to_bytes(8, 'big'))
        h.update(data)

    def walk(x):
        put(type(x).__module__ + '.' + type(x).__qualname__)
        if x is None or isinstance(x, (str, bytes, bool, int, float, complex, np.generic)):
            put(repr(x)); return
        if isinstance(x, (types.ModuleType, types.FunctionType, types.BuiltinFunctionType, type)):
            put(getattr(x, '__module__', '')); put(getattr(x, '__qualname__', getattr(x, '__name__', ''))); return
        if isinstance(x, np.dtype):
            put(x.str); walk(x.descr); return
        ident = id(x)
        if ident in seen:
            put('reference'); put(seen[ident]); return
        seen[ident] = len(seen)
        keepalive.append(x)
        if isinstance(x, torch.Tensor):
            put(str(x.dtype)); put(tuple(x.shape)); put(x.requires_grad)
            a = x.detach().cpu().contiguous()
            if a.is_sparse:
                a = a.to_dense()
            put(a.reshape(-1).view(torch.uint8).numpy().tobytes()); return
        if isinstance(x, np.ndarray):
            put(str(x.dtype)); put(x.shape)
            if x.dtype.hasobject:
                walk(x.tolist())
            else:
                put(np.ascontiguousarray(x).tobytes())
            return
        if isinstance(x, dict):
            for key in sorted(x, key=lambda k: (type(k).__name__, str(k))):
                walk(key); walk(x[key])
            return
        if isinstance(x, (list, tuple)):
            for item in x: walk(item)
            return
        if isinstance(x, (set, frozenset)):
            for item in sorted(x, key=str): walk(item)
            return
        if isinstance(x, torch.nn.Module):
            walk(x.state_dict())
            walk([(name, m.training) for name, m in x.named_modules()])
            # Architecture/configuration and non-persistent buffers also affect
            # inference. Include ordinary attributes, avoiding duplicate graph.
            walk({k: v for k, v in vars(x).items() if k not in ('_parameters', '_modules', '_buffers')})
            walk(x._buffers)
            for name, child in x.named_children(): put(name); walk(child)
            return
        if isinstance(x, np.random.RandomState): walk(x.get_state()); return
        if isinstance(x, np.random.Generator): walk(x.bit_generator.state); return
        if isinstance(x, torch.Generator): walk(x.get_state()); return
        if hasattr(x, '__dict__'):
            walk(vars(x)); return
        # Immutable configuration types have stable text, but unknown opaque
        # objects must fail closed rather than pretending their state was read.
        if isinstance(x, (slice, range, PurePath, torch.device, torch.dtype)):
            put(str(x)); return
        unsupported.add(type(x).__module__ + '.' + type(x).__qualname__)

    try:
        walk(obj)
        if unsupported:
            raise TypeError('Unsupported state types: ' + ', '.join(sorted(unsupported)))
        return h.hexdigest()
    finally:
        # walk's recursive closure survives until cyclic GC. Its keepalive list
        # must not retain a fitted estimator or tensor storage after hashing.
        keepalive.clear()
        seen.clear()


def freeze(obj):
    import torch
    seen = set()
    def visit(x):
        if id(x) in seen: return
        seen.add(id(x))
        if isinstance(x, torch.nn.Module):
            x.eval()
            for p in x.parameters(): p.requires_grad_(False)
        elif isinstance(x, dict):
            for v in x.values(): visit(v)
        elif isinstance(x, (list, tuple)):
            for v in x: visit(v)
        elif hasattr(x, '__dict__') and not isinstance(x, (types.ModuleType, type)):
            for v in vars(x).values(): visit(v)
    visit(obj)


def clone_for_inference(obj, freeze_state=True):
    import torch
    # Estimator prediction can populate caches or write query-specific fields.
    # A private copy retains fitted parameters/context and discards these writes.
    # deepcopy can call backend __setstate__: TabICL reconstructs a randomly
    # initialized network before loading checkpoint weights. Preserve all torch
    # generators even if copying/freezing fails, without reseeding the copy.
    devices = list(range(torch.cuda.device_count())) if torch.cuda.is_available() else []
    with torch.random.fork_rng(devices=devices):
        cloned = copy.deepcopy(obj)
        # Low-level estimators keep their existing evaluation/gradient flags so
        # copying cannot select a different prediction path. Deployment explicitly
        # freezes its copies and supplies the no-grad inference context.
        if freeze_state:
            freeze(cloned)
        return cloned


@contextmanager
def inference(seed):
    import torch
    r, n = random.getstate(), np.random.get_state()
    cpu = torch.random.get_rng_state()
    cuda = torch.cuda.get_rng_state_all() if torch.cuda.is_available() else None
    try:
        random.seed(int(seed)); np.random.seed(int(seed)); torch.manual_seed(int(seed))
        with torch.no_grad():
            yield
    finally:
        random.setstate(r); np.random.set_state(n); torch.random.set_rng_state(cpu)
        if cuda is not None: torch.cuda.set_rng_state_all(cuda)
