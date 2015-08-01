__version__ = "0.4.1"

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
                "  https://github.com/PaulMag/cirbindis"
            )
        if sys.argv[1] in ("version", "Version", "-v", "-V", "--version"):
            print "CirBinDis version %s" % __version__
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
    input_ = func.extract_booleans(input_)


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

            have_resaved = False

            for H0 in func.to_list(input_["H0"], float):
                for diskmass in func.to_list(input_["diskmass"], float):
                    if len(func.to_list(input_["diskmass"], float)) > 1:
                        dataname = input_["dataname"]+"__diskmass=%g" % diskmass
                    else:
                        dataname = input_["dataname"]

                    dataset = DensityMap(
                        filename=input_["datafile"],
                        dataname=dataname,
                        coordsystem=input_["system"],
                        unit=input_["unit"],
                        inclinations=func.to_list(input_["inclination"], float),
                        radius_in=radius_in,
                        radius_out=radius_out,
                        diskmass=diskmass,
                        diskradius=float(input_["diskradius"]),
                        H0=H0,
                        R0=float(input_["R0"]),
                        H_power=float(input_["H_power"]),
                        kappa=float(input_["kappa"]),
                    )

                    if not have_resaved:
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
                        n_gridz=int(input_["sylindergridz"]),
                        unit=input_["unit"]["angle"],
                        lcurve_show=input_["lightcurves"]["show_plot"],
                        lcurve_savefig=input_["lightcurves"]["save_plot"],
                        lcurve_savecsv=input_["lightcurves"]["save_csvtable"],
                        dprofile_show=input_["densityprofiles"]["show_plot"],
                        dprofile_savefig=input_["densityprofiles"]["save_plot"],
                        normalizations=func.to_list(input_["normalization"]),
                        outfolder=input_["outfolder"],
                    )
