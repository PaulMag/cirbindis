import numpy as np

def to_list(x, dtype=None):
    try:
        iter(x)
        if isinstance(x, basestring):
            raise TypeError
        x = np.array(x, dtype=dtype)
    except TypeError:
        x = np.array([x], dtype=dtype)
    return x
