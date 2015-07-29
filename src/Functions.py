import os
import numpy as np


positive = "yes", "true"
negative = "no", "false"
exclude  = "exclude", "never"
def extract_booleans(d):
    """Convert yes/true/no/false strings in a (nested) dict to actual booleans.

    d: (dict) The input dict that will have its strings converted to
        booleans.  Strings that does not have a boolean meaning will keep
        their original value.
    return: (dict) Returns the input dict with the same shape.
    """
    for key, value in d.iteritems():
        if isinstance(value, dict):
            extract_booleans(value)
        else:
            if value in positive:
                d[key] = True
            elif value in negative:
                d[key] = False
            elif value in exclude:
                d[key] = None
    return d


def to_list(x, dtype=None, separator=" "):
    """Converts any sequence or non-sequence into an array.

    The use of a numpy array is primarely used because then it is easy to
    also convert the type. Except for that it could just as well be a list.

    If x is a string, split the string with the separator.
    If x is a sequence, convert it into an array.
    If x is a non-sequence, put it into a size-1 array.

    x: (anything)  Something to be converted into an array.
    dtype: (type) What type of objects the array contains. F.ex. float. If
        None, numpy will interpret the type itself.
    separator: (string) If x is a string to be split, this is the separator.
        Usually a space or comma.
    """

    if isinstance(x, basestring):  # If x is a string.
        x = np.array(x.split(separator), dtype=dtype)
    elif isinstance(x, dict):  # If x is a dictionairy.
        x = np.array([x], dtype=dtype)
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


def make_folder(directory, warning=False):
    """Make the directory, but do nothing if it already exists. Optionally
    print a warning if it already exists, if there should be a reason for
    that.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        if warning:
            print "%s already exists." % (directory)
