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
    # Hash content, not the alias graph reconstructed by pickle. Memoized
    # subtree hashes avoid repeatedly reading shared tensor storage; a repeated
    # reference contributes the same hash as an equal independent value.
    import torch
    memo, active, keepalive, unsupported = {}, {}, [], set()

    def walk(x):
        ident = id(x)
        if ident in active:
            return hashlib.sha256(('cycle:' + str(len(active) - active[ident])).encode()).digest()
        if ident in memo:
            return memo[ident]
        h = hashlib.sha256()
        def put(value):
            data = value if isinstance(value, bytes) else str(value).encode('utf-8')
            h.update(len(data).to_bytes(8, 'big')); h.update(data)
        def child(value):
            put(walk(value))
        put(type(x).__module__ + '.' + type(x).__qualname__)
        active[ident] = len(active)
        keepalive.append(x)
        try:
            if x is None or isinstance(x, (str, bytes, bool, int, float, complex, np.generic)):
                put(repr(x))
            elif isinstance(x, (types.ModuleType, types.FunctionType, types.BuiltinFunctionType, type)):
                put(getattr(x, '__module__', ''))
                put(getattr(x, '__qualname__', getattr(x, '__name__', '')))
            elif isinstance(x, np.dtype):
                put(x.str); child(x.descr)
            elif isinstance(x, torch.Tensor):
                put(str(x.dtype)); put(tuple(x.shape))
                # Autograd ownership is reconstructed on load and does not
                # change frozen inference. Evaluation/training flags below do.
                a = x.detach().cpu()
                if a.is_sparse:
                    a = a.to_dense()
                put(a.contiguous().reshape(-1).view(torch.uint8).numpy().tobytes())
            elif isinstance(x, np.ndarray):
                put(str(x.dtype)); put(x.shape)
                if x.dtype.hasobject:
                    child(x.tolist())
                else:
                    put(np.ascontiguousarray(x).tobytes())
            elif isinstance(x, dict):
                for key in sorted(x, key=lambda k: (type(k).__name__, str(k))):
                    child(key); child(x[key])
            elif isinstance(x, (list, tuple)):
                for item in x: child(item)
            elif isinstance(x, (set, frozenset)):
                for value in sorted(walk(item) for item in x): put(value)
            elif isinstance(x, torch.nn.Module):
                child(x.state_dict())
                child([(name, m.training) for name, m in x.named_modules()])
                attributes = {k: v for k, v in vars(x).items()
                              if k not in ('_parameters', '_modules', '_buffers')}
                # TabPFN's PerFeatureEncoderLayer.__setstate__ inserts this
                # backward-compatibility default even when it was absent at fit.
                if (type(x).__module__ == 'tabpfn.architectures.base.layer'
                        and type(x).__name__ == 'PerFeatureEncoderLayer'):
                    attributes.setdefault('save_peak_mem_factor', None)
                child(attributes)
                child(x._buffers)
                for name, module in x.named_children(): put(name); child(module)
            elif isinstance(x, np.random.RandomState): child(x.get_state())
            elif isinstance(x, np.random.Generator): child(x.bit_generator.state)
            elif isinstance(x, torch.Generator): child(x.get_state())
            elif (type(x).__module__ == 'tabicl._model.inference'
                  and type(x).__name__ == 'InferenceManager'):
                # Reconfigured before every TabICL forward pass from the
                # estimator's inference_config_ (which is hashed separately).
                # Buffer pools, device streams and memory probes are runtime
                # workspaces, not the frozen context/KV cache or model weights.
                child({key: getattr(x, key) for key in ('enc_name', 'out_dim', 'out_no_seq')})
            elif (type(x).__module__.startswith('tabicl._sklearn.')
                  and type(x).__name__ in ('TabICLClassifier', 'TabICLRegressor')):
                # Only checkpoint provenance: __setstate__ either rebuilds this
                # host-specific path or loads the embedded, fully hashed weights.
                child({k: v for k, v in vars(x).items() if k != 'model_path_'})
            elif hasattr(x, '__dict__'): child(vars(x))
            elif isinstance(x, (slice, range, PurePath, torch.device, torch.dtype)):
                put(str(x))
            else:
                unsupported.add(type(x).__module__ + '.' + type(x).__qualname__)
            result = h.digest()
            memo[ident] = result
            return result
        finally:
            active.pop(ident, None)

    try:
        result = walk(obj)
        if unsupported:
            raise TypeError('Unsupported state types: ' + ', '.join(sorted(unsupported)))
        return result.hex()
    finally:
        # Recursive closures must not retain estimator/tensor storage until GC.
        keepalive.clear(); memo.clear(); active.clear()


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


class CopyTracker:
    """Weak ownership audit for one R scoring frame, including nested models."""
    def __init__(self):
        self.references = []
        self.untrackable = False

    def track(self, obj):
        import weakref
        try:
            self.references.append(weakref.ref(obj))
        except TypeError:
            # General Python-backed deployments may include non-weakrefable
            # containers. They require the conservative full R collection.
            self.untrackable = True

    def pending(self):
        return self.untrackable or any(ref() is not None for ref in self.references)

    def collect_cycles(self):
        import gc
        gc.collect()


def clear_exception_frames(error):
    """Keep exception diagnostics without retaining failed inference locals."""
    import traceback
    seen = set()
    pending = [error]
    while pending:
        error = pending.pop()
        if error is None or id(error) in seen:
            continue
        seen.add(id(error))
        traceback.clear_frames(error.__traceback__)
        pending.extend((error.__cause__, error.__context__))
