#mpld3 hack to avoid errors where numpy arrays can't be JSON serialized.
import json

from mpld3 import _display
import numpy as np

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
_display.NumpyEncoder = NumpyEncoder
