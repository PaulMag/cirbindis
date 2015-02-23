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


def points_to_image(data):
    """Convert a set of datapoints into an image through interpolation.

    TODO: Only makes 2D image at the moment. Extend to 3D (not difficult,
    but resource heavy).

    TODO: Only convert a certain slice of the dataset, to limit computation
    time and memory used.

    TODO: Since the coordinates are not included in the resulting image they
    should also be returned somehow as metadata to the image.

    data: (float, array) The dataset to be rotated. Array of shape (N, 4),
        where N is the number of datapoints.

    return: (float, array) Array (image) of shape (n, n) of the density map.
    """

    N = data.shape[0]
    n = round(np.sqrt(N) * 4)
        # Resulting image will use same amount of memory as the datapoints.

    grid_x, grid_y, grid_z = np.mgrid[
        - radius_out : radius_out : n*1j,
        - radius_out : radius_out : n*1j,
        -.5*thickness:.5*thickness: (2*n_layers+1)*1j,
    ]
    grid_density = scipy.interpolate.griddata(
        data[:, 0:3],
        data[:, 3],
        (grid_x, grid_y, grid_z),
        method='linear',
    )
    return grid_density



if __name__ == "__main__":
    """Everything under this should just be considered a test block for now."""

    radius_star = 0.01
    radius_in = 0.01
    radius_out = 0.03
    n_layers = 1
    thickness = 0.2
    # filename = "../data/data_tiny_3d.p"
    filename = "../data/data_micro_3d.p"

    data = load(filename, method="pickle", \
        radius_in=radius_in, radius_out=radius_out)
    # writeto(data, "../data/data_micro.tab", method="ascii")
    # data = add_3d_points(data, H=1, n_layers=1, dz=0.1)
    # writeto(data, "../data/data_micro_3d.p")

    img = points_to_image(data)
    print img.shape
    plt.imshow(img[:, :, 0], interpolation="nearest", origin="lower")
    plt.show()

    import sys; sys.exit(0)

    # plt.plot(
        # data[::1, 0],
        # data[::1, 1],
        # "r+",
    # )
    # plt.show()
    # data = rotate(data, angle_x=25)
    # plt.plot(
        # data[::1, 0],
        # data[::1, 1],
        # "r+",
    # )
    # plt.show()
