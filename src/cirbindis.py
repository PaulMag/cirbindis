__version__ = "0.1"

import sys
import numpy as np
import xmltodict

from DensityMap import DensityMap
from Star import Star
import Functions as func


if __name__ == "__main__":

    # Attempt to read input file and tell the user if the input is wrong:
    try:
        infile = open(sys.argv[1], "r")
    except IndexError:
        print (
            "usage: python %s input_file.xml\n"
            "  Specify the (path)name of an input xml file. "
            % sys.argv[0]
        )
        sys.exit(1)
    except IOError:
        print "No such file: '%s'" % sys.argv[1]
        print (
            "usage: python %s input_file.xml\n"
            "  Specify the (path)name of an input xml file. "
            % sys.argv[0]
        )
        if sys.argv[1] in (
            "help", "Help", "-h", "-H", "--help", "man", "manual"
        ):
            print (
                "  See the usermanual for details on how to use cirbindis:\n"
                "  /doc/cirbindis_usermanual.pdf\n"
                "  or\n"
                "  https://bitbucket.org/paulmag/circumbinarydisk"
            )
        sys.exit(1)
    try:
        input_ = xmltodict.parse(infile)["input"]
    except Exception:
        print (
            "usage: python %s input_file.xml\n"
            "  Specify the (path)name of an input xml file. "
            % sys.argv[0]
        )
        raise
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

            if radius_out is None or np.isnan(radius_out):
                radius_out = np.inf

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

            dataset.set_physical_units()

            for star in func.to_list(input_["star"]):
                dataset.add_star(star)

            dataset.make_lightcurve(
                n_angle=int(input_["azimuthsteps"]),
                n_radius=int(input_["radiussteps"]),
                unit=input_["unit"]["angle"],
                save=True,
                savefig=True,
                show=True,
                normalizations=func.to_list(input_["normalization"]),
                outfolder=input_["outfolder"],
            )
