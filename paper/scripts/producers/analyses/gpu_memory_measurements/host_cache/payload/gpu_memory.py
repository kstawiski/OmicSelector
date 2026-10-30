"""Keep the exact float32 KV cache in host RAM; transfer one ensemble slice per CUDA call."""
import functools,json
from collections import OrderedDict

def configure():
    from tabicl import TabICLClassifier
    from tabicl._model.kv_cache import TabICLCache
    if getattr(TabICLClassifier.__init__, '_pb1_host_cache', False):
        return
    original_init = TabICLClassifier.__init__
    @functools.wraps(original_init)
    def bounded(self, *args, **kwargs):
        kwargs['batch_size'] = 1
        return original_init(self, *args, **kwargs)
    bounded._pb1_host_cache = True
    TabICLClassifier.__init__ = bounded
    concat = TabICLCache.concat
    def cpu_concat(caches, dim=0):
        return concat([cache.to('cpu') for cache in caches], dim=dim)
    TabICLCache.concat = staticmethod(cpu_concat)
    slice_batch = TabICLCache.slice_batch
    def cuda_slice(self, start, end):
        return slice_batch(self, start, end).to('cuda')
    TabICLCache.slice_batch = cuda_slice
    # Called by backend deserialization, including its private inference copies.
    def restore_host_cache(self):
        if getattr(self, 'model_kv_cache_', None) is not None:
            self.model_kv_cache_ = OrderedDict((k,c.to('cpu')) for k,c in self.model_kv_cache_.items())
    TabICLClassifier._move_cache_to_device = restore_host_cache

def audit_json(clf):
    import torch
    tensors=[]
    seen=set()
    def walk(x):
        if id(x) in seen:return
        seen.add(id(x))
        if isinstance(x,torch.Tensor):tensors.append(x)
        elif isinstance(x,dict):
            for y in x.values():walk(y)
        elif isinstance(x,(list,tuple)):
            for y in x:walk(y)
        elif hasattr(x,'__dict__'):
            for y in vars(x).values():walk(y)
    walk(clf.model_kv_cache_)
    assert tensors and all(t.device.type=='cpu' and t.dtype==torch.float32 for t in tensors)
    assert clf.batch_size==1 and clf.n_estimators==8 and clf.use_amp is False and clf.kv_cache is True and clf.offload_mode is False
    assert next(clf.model_.parameters()).device.type=='cuda'
    return json.dumps(dict(cache_device='cpu',cache_tensor_count=len(tensors),cache_bytes=sum(t.numel()*t.element_size() for t in tensors),cache_dtype='torch.float32',inference_device=str(next(clf.model_.parameters()).device),ensemble_batch_size=clf.batch_size,n_estimators=clf.n_estimators,use_amp=clf.use_amp,kv_cache=clf.kv_cache,column_embedding_offload=clf.offload_mode),sort_keys=True)
