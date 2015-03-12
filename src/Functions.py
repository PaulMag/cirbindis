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

def cart2pol(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return [r, theta]

def pol2cart(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return [x, y]
