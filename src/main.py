import xmltodict

from DensityMap import DensityMap
import Functions as func


if __name__ == "__main__":

    infile = open("input.xml", "r")
    input_ = xmltodict.parse(infile)["input"]
    infile.close()

    for radius_out in func.to_list(input_["radius_out"], float):
        dataset = DensityMap(
            filename=input_["datafile"],
            inclinations=func.to_list(input_["inclination"], float),
            radius_in=float(input_["radius_in"]),
            radius_out=radius_out,
            kappa=float(input_["kappa"]),
        )
        for star in func.to_list(input_["star"]):
            dataset.add_star(star)
        dataset.make_lightcurve(
            n_angle=int(input_["azimuthsteps"]),
            n_radius=int(input_["radiussteps"]),
            unit=input_["unit"]["angle"],
            show=True,
        )
        # dataset.set_r0()
        # dataset.set_H()
        # dataset.set_sigma()
        # dataset.get_density_profile()

    # dataset.writeto(filename)
