"""Canonical CPU encoders and column embeddings for the pinned TabPFN-v2 model."""
import torch
from tabpfn.architectures.base.transformer import PerFeatureTransformer


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
    for model in clf.models_:
        if not isinstance(model, PerFeatureTransformer):
            raise TypeError("Expected the pinned TabPFN-v2 PerFeatureTransformer")
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
