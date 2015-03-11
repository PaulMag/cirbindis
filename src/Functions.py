def to_list(x):
    try:
        iter(x)
    except TypeError:
        x = [x]
    return x
