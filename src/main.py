from DensityMap import DensityMap


if __name__ == "__main__":
    """Everything under this should just be considered a test block for now."""

    radius_star = .75
    radius_in = 1.0
    radius_out = 3.5
    n_radius = 10
    n_angle = 6
    inclinations = [0, 30]
    unit = "deg"
    kappa = 1.
    H = 1.
    dz = .2
    ratio = .2
    filename = "../data/data_cropped.p"

    for radius_out in [3.0]:
        dataset = DensityMap(
            filename=filename,
            radius_star=radius_star,
            radius_in=radius_in,
            radius_out=radius_out,
        )
        dataset.make_lightcurve(
            inclinations=inclinations,
            theta=360.,
            n_angle=n_angle,
            n_radius=n_radius,
            unit=unit,
            dz=dz,
            ratio=ratio,
            show=True,
        )

    # dataset.writeto(filename)
