import sys
import numpy as np
    #  Numerical Python. For performing computations on arrays/matrices fast.
import cPickle as pickle
    # Can be used to save (dump) and load Python objects to files. Much
    # faster than reading and writing ASCII tables.
import time
    # Used to time parts of code to look for bottlenecks.
import scipy.interpolate
import matplotlib.pyplot as plt


def load(
    filename,
    radius_in=0,
    radius_out=np.inf,
    method="pickle",
    separator=" ",
):
    """Load a dataset to analyse from a file.

    Assumed to be on the form 'x,y,density' or 'x,y,z,density' which
    represents a point in cartesian space and the density at that point. If
    the z-coordinate is not given it will be assumed to be 0 for every point.

    filename: (string) Full pathname to file containing dataset.
    radius_in: (float) Crop points that are closer than this to origo.
    radius_out: (float) Crop points that are further than this to origo.
    method: (string) What kind of loading algorithm to use.
    separator: (string) If method='ascii' this is the separator
        between the values each line. Usually a space or comma. Ignored if
        method='pickle'.

    return: (float, array) Array of shape (N, 4), where N is the number of
        data points.
    """

    t_start = time.time()  # Just to time the loading, in case of large dataset.
    infile = open(filename, "r")

    if method == "pickle":
        data = pickle.load(infile)
        mask = (
            (np.linalg.norm(data[:, 0:2], axis=1) >= radius_in) *
            (np.linalg.norm(data[:, 0:2], axis=1) <= radius_out)
        )
        data = data[np.where(mask)]

    elif method == "ascii":
        data = []
        for line in infile:
            line = [float(value) for value in line.split(separator)]
            if len(line) >= 3:
                if radius_in <= np.linalg.norm(line[0:2]) <= radius_out:
                    data.append()
        data = np.array(data)

    infile.close()
    t_end = time.time()  # End of timer.
    print "Loading took %f seconds." % (t_end - t_start)

    if data.shape[1] < 4:
        # Add the z-dimension if it is not already there.
        z = np.zeros((data.shape[0], 1))
        data = np.hstack((data[:, 0:2], z, data[:, 2, None]))
    return data


def writeto(data, filename, method="pickle", separator=" "):
    """Write a dataset to a file for later use.

    Assumed to be on the form 'x,y,z,density' which represents a point in
    cartesian space and the density at that point.

    filename: (string) Full pathname to outfile for writing data.
    method: (string) What kind of writing algorithm to use. Recommended to
        use 'pickle' if it will be loaded by this program later (faster) and
        'ascii' for an other purpose.
    separator: (string) If method='ascii' this is the separator between the
        values each line. Usually a space or comma. Ignored if method='pickle'.
    """

    t_start = time.time()  # Just to time the writing, in case of large dataset.

    if method == "pickle":
        outfile = open(filename, "wb")
        pickle.dump(data, outfile)

    elif method == "ascii":
        outfile = open(filename, "w")
        for line in data:
            outfile.write("%f%s%f%s%f%s%f\n" % (
                line[0], separator,
                line[1], separator,
                line[2], separator,
                line[3],
            ))

    outfile.close()
    t_end = time.time()  # End of timer.
    print "Writing took %f seconds." % (t_end - t_start)


def rotate(data, angle_x=0, angle_y=0, angle_z=0, unit="deg"):
    """Rotate entire dataset by an angle around any axis.

    data: (float, array) The dataset to be rotated. Array of shape (N, 4),
        where N is the number of datapoints.
    angle_x: (float) Angle to rotate around x-axis (roll).
    angle_y: (float) Angle to rotate around y-axis (pitch).
    angle_z: (float) Angle to rotate around z-axis (yaw).
    unit: (string) What unit angles are given in.
        'rad', 'deg', 'arcmin' or 'arcsec'.

    return: (float, array) The input data with coordinates rotated.
    """

    if unit == "rad":
        factor = 1.
    elif unit == "deg":
        factor = np.pi / 180.
    elif unit == "arcmin":
        factor = np.pi / 180. * 60
    elif unit == "arcsec":
        factor = np.pi / 180. * 3600
    angle_x *= factor
    angle_y *= factor
    angle_z *= factor

    R_x = np.matrix([
        [             1,                0,                0],
        [             0,  np.cos(angle_x), -np.sin(angle_x)],
        [             0,  np.sin(angle_x),  np.cos(angle_x)],
    ])
    R_y = np.matrix([
        [ np.cos(angle_y),              0,  np.sin(angle_y)],
        [               0,              1,                0],
        [-np.sin(angle_y),              0,  np.cos(angle_y)],
    ])
    R_z = np.matrix([
        [ np.cos(angle_z), -np.sin(angle_z),              0],
        [ np.sin(angle_z),  np.cos(angle_z),              0],
        [               0,                0,              1],
    ])
    rotation_matrix = R_z * R_y * R_x

    coords_in = data[:, 0:3]
    coords_out = (rotation_matrix * coords_in.transpose()).transpose()
    data[:, 0:3] = coords_out
    return data


def add_3d_points(data, H, n_layers=None, dz=None, thickness=None):
    """Expands a 2D-disk into 3rd (z) dimension assuming a simple model.

    z-coordinate must be 0 for every point in data. If the dataset is already
    3-dimensional then you NEED not and MUST not use this method. You must
    use this method BEFORE any rotation around x- or y-axis, NOT AFTER.
    Rotations around z-axis does not matter, as they do not alter z-data.

    Note that this method multiplies the size of the dataset provided, so
    the returned data can be very large.

    data: (float, array) The dataset to be rotated. Array of shape (N, 4),
        where N is the number of datapoints.
    H: (float) A physical parameter. We assume it to be a constant.
    n_layers: (int) How many data layers to add to EACH side of the disk.
        The total number of layers in z-direction becomes 2*n_layers+1.
    dz: (float) Distance between each layer.
    thickness: (float) Total thickness of disk to provide data for. This
        argument can be given instead of either n_layers or dz. If both
        n_layers and dz is provided then this is ignored.

    return: (float, array) The dataset with 2*n_layers*N more points added.
    """

    if n_layers is None:
        n_layers = int(round(0.5 * thickness / dz))
    elif dz is None:
        dz = 0.5 * thickness /n_layers

    N = data.shape[0]
    data_over = np.zeros((N*n_layers, 4))

    for k in xrange(1, n_layers+1):
        data_over[(k-1)*N : k*N, 0:2] = data[:, 0:2]
        data_over[(k-1)*N : k*N, 2] = k * dz
        data_over[(k-1)*N : k*N, 3] = data[:, 3] * np.exp(- k * dz / H)
    data_under = data_over.copy()
    data_under[:, 2] *= -1

    return np.vstack((data, data_over, data_under))


def points_to_image(data, sylinder=True):
    """Convert a set of datapoints into a datacube through linear interpolation.

    TODO: Since the coordinates are not included in the resulting image they
    should also be returned somehow as metadata to the image.

    data: (float, array) The dataset to be rotated. Array of shape (N, 4),
        where N is the number of datapoints.
    sylinder: (bool) If True, will only convert the relevant sylinder. If
        False, will convert the entire input dataset.

    return: (float, array) Array (datacube) of shape (n, n) of the density map.
    """

    N = data.shape[0]

    if sylinder:
        # n_x = round(N / ((radius_out-radius_in) / (2*radius_star*n_layers)))
        # y = round(N * (2*radius_star) / (radius_out-radius_in) / n_layers)
        n_x = 7
        n_y = 3  # Or the datacube gets too big! :(
        grid_x, grid_y, grid_z = np.mgrid[
            radius_in     : radius_out  : n_x*1j,
            - radius_star : radius_star : n_y*1j,
            - radius_star : radius_star : n_y*1j,
        ]
    else:
        # This is not really needed. Maybe for informative plots for the user.
        n = round(np.sqrt(N) * 4)
        grid_x, grid_y, grid_z = np.mgrid[
            - radius_out : radius_out : n*1j,
            - radius_out : radius_out : n*1j,
            -.5*thickness:.5*thickness: (2*n_layers+1)*1j,
        ]
        # grid_x, grid_y, grid_z = np.mgrid[
            # data[:, 0].min() : data[:, 0].max() : n*1j,
            # data[:, 1].min() : data[:, 1].max() : n*1j,
            # data[:, 2].min() : data[:, 2].max():  n*1j,
    datacube = scipy.interpolate.griddata(
        data[:, 0:3],
        data[:, 3],
        (grid_x, grid_y, grid_z),
        fill_value=0,
        method='linear',
    )
    return datacube


def get_sylinder(data):
    """TODO: Write docstring."""
    mask = (
        (data[:, 0] > 0) *
        (np.linalg.norm(data[:, 1:3], axis=1) <= radius_star)
    )
    data_sylinder = data[np.where(mask)]
    return data_sylinder


def space_sylinder(data, n_steps=None, dr=None):
    """TODO: Write docstring."""

    if n_steps is None:
        n_steps = int(round((radius_out-radius_in) / dr))
    elif dr is None:
        dr = (radius_out-radius_in) / n_steps

    radiuses = np.linspace(radius_in, radius_out, n_steps+1)
        # TODO Should not be made here, since it is the same always.
    densities = np.zeros(n_steps)

    for i in xrange(n_steps):
        mask = (
            (data[:, 0] >  radiuses[i]) *
            (data[:, 0] <= radiuses[i+1])
        )
        if data[np.where(mask), 3].size == 0:
            densities[i] = 0
            print "Warning: No stars in bin (%g, %g]. Density set to 0." \
                (densities[i], densities[i+1])
        elif:
            densities[i] = data[np.where(mask), 3].mean()
    return densities


def integrate(densities):
    """TODO: Write docstring."""

    n_steps = densities.shape[0]  # Assume sylinder is oriented in x-direction.
    dr = (radius_out - radius_in) / n_steps

    intensity = 1.  # Or whatever the full intensity of the star is.
    kappa = 1  #TODO Should be defined as a constant outside this function.

    for i in xrange(n_steps):
        tau = kappa * densities[i] * dr
        intensity *= np.exp(-tau)

    return intensity


def make_lightcurve(data, n_angle=None, dtheta=None, theta=None, unit="deg", n_radius=None, dr=None):
    """TODO: Write docstring."""

    if n_angle is None:
        n_angle = int(round(float(theta) / dtheta))
    elif dtheta is None:
        dtheta = float(theta) / n_angle

    angles = np.linspace(0, theta-dtheta, n_angle)
    lightcurve = np.zeros(n_angle)

    for i, angle in enumerate(angles):
        print "%f / %f" % (angle, theta)
        lightcurve[i] = integrate(space_sylinder(
            get_sylinder(rotate(data, angle_z=angle, unit=unit)),
            n_steps=n_radius,
            dr=dr,
        ))
    print "%f / %f" % (theta, theta)

    plt.plot(angles, lightcurve)
    plt.show()


if __name__ == "__main__":
    """Everything under this should just be considered a test block for now."""

    radius_star = 1.0
    radius_in = 1.0
    radius_out = 4.0
    n_layers = 1
    thickness = 0.002
    # filename = "../data/data_tiny_3d.p"
    filename = "../data/data_cropped.p"

    data = load(filename, method="pickle", \
        radius_in=radius_in, radius_out=radius_out)
    data = add_3d_points(data, H=1, n_layers=n_layers, thickness=thickness)

    make_lightcurve(data, theta=360., n_angle=360, n_radius=30)

    # writeto(data, "../data/data_micro_3d.p")

    # plt.plot(
        # data[::1, 0],
        # data[::1, 1],
        # "r+",
    # )
    # plt.show()
    # data = rotate(data, angle_z=25)
    # plt.plot(
        # data[::1, 0],
        # data[::1, 1],
        # "r+",
    # )
    # plt.show()
