import numpy as np

class Star:

    def __init__(self, position, radius, intensity):
        self.position = np.array(position)
        self.radius = radius
        self.intensity = intensity
