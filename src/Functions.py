import numpy as np


def to_list(x, dtype=None, separator=" "):
    """Converts any sequence or non-sequence into an array.

    The use of a numpy array is primarely because then it is easy to also convert the type. Except for that it couls just as well be a list.

    If x is a string, split the string with the separator.
    If x is a sequence, convert it into an array.
    If x is a non-sequence, but it into a size-1 array.

    X: (anything)  Something to be converted into an array.
    dtype: (type) What type of objects the array contains. F.ex. float. If
        None, numpy will interpret the type itself.
    separator: (string) If x is a list to be splitted, this is the
        separator. Normally a space or comma.
    """

    if isinstance(x, basestring):  # If x is a string.
        x = np.array(x.split(separator), dtype=dtype)
    else:
        try:
            iter(x)  # If x is a sequence.
            x = np.array(x, dtype=dtype)
        except TypeError:  # If x is a non-sequence.
            x = np.array([x], dtype=dtype)

    return x


def cart2pol(x, y):
    """Convert cartesian coordinates into polar coordinates."""
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return [r, theta]


def pol2cart(r, theta):
    """Convert polar coordinates into cartesian coordinates."""
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return [x, y]
