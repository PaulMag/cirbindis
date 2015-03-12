import numpy as np

import Functions as func



class Star:


    def __init__(self, d=None, position=None, radius=None, intensity=None):
        """Make a star instance.

        d: (dictionairy) Must contain position, radius and intensity and can
            be provided instead of giving these other arguments individually.
        position: (float, array-like) Coordinates of the star.
        radius: (float) Radius of the star.
        intensity: (float) Intensity of the star.
        """

        if d is not None:
            try:
                position = np.array([
                    float(d["position"]["x"]),
                    float(d["position"]["y"]),
                    0,
                ])
            except KeyError, ValueError:
                position = np.array(func.pol2cart(
                    float(d["position"]["r"]),
                    float(d["position"]["theta"]),
                ) + [0])
            radius = float(d["radius"])
            intensity = float(d["intensity"])

        self.position = np.array(position)
        self.radius = radius
        self.intensity = intensity
