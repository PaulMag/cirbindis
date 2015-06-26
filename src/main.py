import numpy as np
import xmltodict

from DensityMap import DensityMap
from Star import Star
import Functions as func


if __name__ == "__main__":

    infile = open("input.xml", "r")
    input_ = xmltodict.parse(infile)["input"]
    infile.close()


    for radius_in in func.to_list(input_["radius_in"], float):

        if radius_in is None or np.isnan(radius_in):
            radius_largest = 0.
            for stardict in func.to_list(input_["star"]):
                star = Star(stardict)
                radius = np.linalg.norm(star.position) + star.radius
                if radius > radius_largest:
                    radius_largest = radius
            radius_in = radius_largest

        for radius_out in func.to_list(input_["radius_out"], float):

            dataset = DensityMap(
                filename=input_["datafile"],
                dataname=input_["dataname"],
                coordsystem=input_["system"],
                unit=input_["unit"],
                inclinations=func.to_list(input_["inclination"], float),
                radius_in=radius_in,
                radius_out=radius_out,
                diskmass=float(input_["diskmass"]),
                diskradius=float(input_["diskradius"]),
                H=float(input_["H0"]),
                kappa=float(input_["kappa"]),
            )

            for filename in func.to_list(input_["resave_as"]):
                if filename is None:
                    continue
                elif len(filename) == 0:
                    continue
                dataset.writeto(filename)

            for star in func.to_list(input_["star"]):
                dataset.add_star(star)

            dataset.make_lightcurve(
                n_angle=int(input_["azimuthsteps"]),
                n_radius=int(input_["radiussteps"]),
                unit=input_["unit"]["angle"],
                save=True,
                show=True,
                normalizations=func.to_list(input_["normalization"]),
                outfolder=input_["outfolder"],
            )
