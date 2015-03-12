import numpy as np

def to_list(x, dtype=None, separator=" "):
    if isinstance(x, basestring):
        x = np.array(x.split(separator), dtype=dtype)
    else:
        try:
            iter(x)
            x = np.array(x, dtype=dtype)
        except TypeError:
            x = np.array([x], dtype=dtype)
    return x
