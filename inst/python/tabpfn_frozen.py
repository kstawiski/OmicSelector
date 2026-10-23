"""Canonical CPU encoders and column embeddings for the pinned TabPFN-v2 model."""
import torch
import inspect
import numpy as np
from tabpfn.preprocessing.ensemble import TabPFNEnsembleMember
from tabpfn.inference import InferenceEngineCacheKV
from foundation_inference import full_precision
from tabpfn_memory import query_limit
from tabpfn.architectures.base.attention.full_attention import MultiHeadAttention
from tabpfn.architectures.base.transformer import PerFeatureTransformer


class FrozenEnsembleMember(TabPFNEnsembleMember):
    """Canonicalize query preprocessing before its floating-point row hash.

    PCA/polynomial transforms can differ below 1e-12 with the query batch shape.
    TabPFN hashes those values into a fingerprint, amplifying tiny differences.
    Use exactly the singleton preprocessing path, then batch neural inference.
    """
    def transform_X_test(self, X):
        if X.shape[0] <= 1:
            return super().transform_X_test(X)
        rows = [super(FrozenEnsembleMember, self).transform_X_test(X[i:i + 1])
                for i in range(X.shape[0])]
        return torch.cat(rows, dim=0) if isinstance(rows[0], torch.Tensor) else np.concatenate(rows, axis=0)


# TabPFN 8.0.8's memory decorator checks .requires_grad on optional None
# arguments before checking no_grad. Keep its chunk arithmetic, but validate
# only tensor inputs. This adapter affects only this fitted model's attention.
_attention_compute = inspect.getclosurevars(MultiHeadAttention._compute).nonlocals["method"]


class FrozenMultiHeadAttention(MultiHeadAttention):
    def _compute(self, x, *args, add_input=False, allow_inplace=False,
                 save_peak_mem_factor=None, **kwargs):
        if save_peak_mem_factor is None:
            result = _attention_compute(self, x, *args, **kwargs)
            return x + result if add_input else result
        if not allow_inplace or not isinstance(save_peak_mem_factor, int) or save_peak_mem_factor <= 1:
            raise ValueError("Memory-saving attention requires inference chunks and an inplace buffer")
        tensors = [x, *self.parameters(), *(a for a in args if isinstance(a, torch.Tensor))]
        if torch.is_grad_enabled() and any(t.requires_grad for t in tensors):
            raise RuntimeError("Memory-saving attention is inference-only")
        size = (x.size(0) + save_peak_mem_factor - 1) // save_peak_mem_factor
        chunks = [torch.split(a, size) if isinstance(a, torch.Tensor)
                  else [a] * save_peak_mem_factor for a in (x, *args)]
        for block, *block_args in zip(*chunks):
            result = _attention_compute(self, block, *block_args, **kwargs)
            if add_input:
                block[:] += result
            else:
                block[:] = result
        return x


def _map_tensors(value, function):
    if isinstance(value, torch.Tensor):
        return function(value)
    if isinstance(value, dict):
        return {key: _map_tensors(item, function) for key, item in value.items()}
    if isinstance(value, tuple):
        return tuple(_map_tensors(item, function) for item in value)
    if isinstance(value, list):
        return [_map_tensors(item, function) for item in value]
    return value


class CPUEncoder(torch.nn.Module):
    """Keep data-dependent normalization on CPU even when the transformer moves."""
    def __init__(self, encoder):
        super().__init__()
        self.encoder = encoder.cpu().float()

    def _apply(self, fn, recurse=True):
        # This fitted preprocessing state intentionally stays float32 on CPU.
        return self

    def forward(self, *args, **kwargs):
        devices = []
        def to_cpu(tensor):
            devices.append(tensor.device)
            return tensor.cpu()
        cpu_args = _map_tensors(args, to_cpu)
        cpu_kwargs = _map_tensors(kwargs, to_cpu)
        result = self.encoder(*cpu_args, **cpu_kwargs)
        device = devices[0] if devices else torch.device("cpu")
        return _map_tensors(result, lambda tensor: tensor.to(device))


class FrozenPerFeatureTransformer(PerFeatureTransformer):
    def add_embeddings(self, x, y, *, data_dags, num_features, seq_len,
                       cache_embeddings=False, use_cached_embeddings=False):
        if data_dags is not None:
            raise ValueError("OmicSelector's frozen TabPFN does not accept data DAGs")
        embedding = self.cpu_column_embeddings[x.shape[2]].to(x.device, x.dtype)
        x = x + embedding[None, None]
        return x, y


def freeze(clf):
    """Freeze every ensemble's CPU-derived encoder and positional feature state."""
    for member in clf.executor_.ensemble_members:
        member.__class__ = FrozenEnsembleMember
    for model in clf.models_:
        if not isinstance(model, PerFeatureTransformer):
            raise TypeError("Expected the pinned TabPFN-v2 PerFeatureTransformer")
        for module in model.modules():
            if isinstance(module, MultiHeadAttention):
                module.__class__ = FrozenMultiHeadAttention
        model.encoder = CPUEncoder(model.encoder)
        model.y_encoder = CPUEncoder(model.y_encoder)
        if hasattr(model, "feature_positional_embedding_embeddings"):
            model.feature_positional_embedding_embeddings = CPUEncoder(
                model.feature_positional_embedding_embeddings)
        counts = {int((member.X_train.shape[1] + model.features_per_group - 1)
                      // model.features_per_group)
                  for member in clf.executor_.ensemble_members}
        model.cpu_column_embeddings = {}
        with torch.no_grad():
            for count in counts:
                x = torch.zeros((1, 1, count, model.emsize), dtype=torch.float32)
                y = torch.zeros((1, 1, model.emsize), dtype=torch.float32)
                result, _ = model.add_embeddings(x, y, data_dags=None,
                    num_features=count * model.features_per_group, seq_len=1)
                model.cpu_column_embeddings[count] = result[0, 0].clone()
        model.__class__ = FrozenPerFeatureTransformer
    return clf


class FrozenCacheKV(InferenceEngineCacheKV):
    """Canonical CPU-built context cache with an explicit inference device.

    Cached tensors are registered attention buffers; fitted encoder statistics
    stay on CPU. Device changes select where private prediction copies execute.
    """
    def iter_outputs(self, X, *, autocast, task_type, only_return_standard_out=True):
        # OmicSelector freezes CPU preprocessing before building this cache.
        # The upstream cache engine requests a GPU-preprocessing cache even
        # when GPU preprocessing is disabled; use the fitted CPU schema here.
        from dataclasses import replace
        from tabpfn.preprocessing.datamodel import FeatureModality
        from tabpfn.inference import DEFAULT_SAVE_PEAK_MEMORY_FACTOR
        if autocast:
            raise ValueError("Frozen TabPFN context requires float32 inference")
        for member, model in zip(self.ensemble_members, self.models, strict=True):
            model.to(self.device)
            x = torch.as_tensor(member.transform_X_test(X), dtype=torch.float32,
                                device=self.device).unsqueeze(1)
            cats = [member.feature_schema.indices_for(FeatureModality.CATEGORICAL)]
            options = replace(model.get_default_performance_options(),
                save_peak_memory_factor=DEFAULT_SAVE_PEAK_MEMORY_FACTOR
                if self.save_peak_mem else None)
            with torch.inference_mode():
                output = model(x, y=None, only_return_standard_out=only_return_standard_out,
                               categorical_inds=cats, performance_options=options)
            model.cpu()
            yield output if isinstance(output, dict) else output.squeeze(1), member.config

    def _move_models_to_devices(self, devices):
        self.device = devices[0]
        # Keep the retained context canonical. iter_outputs moves only the
        # scorer's private copies, and returns their neural buffers to CPU.


def cache_context(clf):
    """Build one immutable transformer context per already-frozen ensemble."""
    from copy import deepcopy
    from tabpfn.inference import SingleDeviceInferenceEngine, DEFAULT_SAVE_PEAK_MEMORY_FACTOR
    from dataclasses import replace
    from tabpfn.preprocessing.datamodel import FeatureModality
    old = clf.executor_
    if isinstance(old, FrozenCacheKV):
        return clf
    engine = FrozenCacheKV.__new__(FrozenCacheKV)
    models = []
    with full_precision(), torch.inference_mode():
        for member in old.ensemble_members:
            model = deepcopy(clf.models_[member.config._model_index]).cpu().float()
            model.cache_trainset_representation = True
            x = torch.as_tensor(member.X_train, dtype=torch.float32, device="cpu").unsqueeze(1)
            y = torch.as_tensor(member.y_train, dtype=torch.float32, device="cpu")
            cats = [member.feature_schema.indices_for(FeatureModality.CATEGORICAL)]
            options = replace(model.get_default_performance_options(),
                save_peak_memory_factor=DEFAULT_SAVE_PEAK_MEMORY_FACTOR
                if old.save_peak_mem else None)
            model(x, y, only_return_standard_out=True, categorical_inds=cats,
                  performance_options=options)
            models.append(model)
    SingleDeviceInferenceEngine.__init__(engine, models=models,
        save_peak_mem=old.save_peak_mem, dtype_byte_size=old.dtype_byte_size,
        force_inference_dtype=old.force_inference_dtype)
    engine.ensemble_members = old.ensemble_members
    engine.device = clf.devices_[0]
    # Cached inference disables the upstream attention memory-saving factor.
    # Bound its dense query-by-context workspace using the actual frozen KV
    # layouts. Keep the original R chunk partition whenever it already fits.
    layouts = []
    for model in models:
        for module in model.modules():
            if isinstance(module, MultiHeadAttention):
                cache = next((getattr(module, name, None)
                              for name in ("_kv_cache", "_k_cache")
                              if isinstance(getattr(module, name, None), torch.Tensor)), None)
                if cache is not None:
                    layouts.append((cache.shape[0], cache.shape[1],
                                    module._nhead, cache.element_size()))
    engine.query_chunk_size = query_limit(layouts)
    clf.executor_ = engine
    clf.fit_mode = "fit_with_cache"
    return clf
