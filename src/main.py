import numpy as np
    #  Numerical Python. For performing computations on arrays/matrices fast.
import cPickle as pickle
    # Can be used to save (dump) and load Python objects to files. Much
    # faster than reading and writing ASCII tables.
import time
    # Used to time parts of code to look for bottlenecks.


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
            outfile.write("%f%s%f%s%f\n" % (
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

    data: (float, array) The dataset to be rotated. Array of shape (N, 4), where N is the number of datapoints.
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


if __name__ == "__main__":

    radius_star = 0.3
    radius_in = 0.5
    radius_out = 3.0
    filename = "../data/data_cropped.p"

    data = load(filename, method="pickle", \
        radius_in=radius_in, radius_out=radius_out)

    import matplotlib.pyplot as plt
    plt.plot(
        data[::1, 0],
        data[::1, 1],
        "r+",
    )
    plt.show()
    data = rotate(data, angle_x=25)
    plt.plot(
        data[::1, 0],
        data[::1, 1],
        "r+",
    )
    plt.show()
