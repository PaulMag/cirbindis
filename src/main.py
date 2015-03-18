import xmltodict

from DensityMap import DensityMap
import Functions as func


if __name__ == "__main__":

    infile = open("input.xml", "r")
    input_ = xmltodict.parse(infile)["input"]
    infile.close()

    for radius_in in func.to_list(input_["radius_in"], float):
        for radius_out in func.to_list(input_["radius_out"], float):
            dataset = DensityMap(
                filename=input_["datafile"],
                coordsystem=input_["system"],
                unit=input_["unit"],
                inclinations=func.to_list(input_["inclination"], float),
                radius_in=radius_in,
                radius_out=radius_out,
                diskmass=float(input_["diskmass"]),
                H=float(input_["H0"]),
                kappa=float(input_["kappa"]),
            )
            for star in func.to_list(input_["star"]):
                dataset.add_star(star)
            dataset.make_lightcurve(
                n_angle=int(input_["azimuthsteps"]),
                n_radius=int(input_["radiussteps"]),
                unit=input_["unit"]["angle"],
                show=True,
                save=True,
                outfolder=input_["outfolder"],
            )
