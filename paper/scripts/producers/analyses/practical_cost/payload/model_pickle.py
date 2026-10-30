"""Persist exact Python-held fitted objects without reticulate's pickle wrapper."""
import pickle
import cloudpickle

def save(obj, filename):
    with open(filename, 'wb') as stream:
        cloudpickle.dump(obj, stream, protocol=pickle.HIGHEST_PROTOCOL)

def load(filename):
    with open(filename, 'rb') as stream:
        return pickle.load(stream)
