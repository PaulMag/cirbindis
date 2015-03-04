import sys
    # For doing meta things like receiving command-line arguments and exiting
    # the program.
import numpy as np
    #  Numerical Python. Contains mathematical functions for performing
    # computations on arrays/matrices fast.
import cPickle as pickle
    # Can be used to save (dump) and load Python objects to files. Much
    # faster than reading and writing ASCII tables.
import time
    # Used to time parts of code to look for bottlenecks.
import matplotlib.pyplot as plt
    # For plotting results.



class DensityMap:


    def __init__(self,
        data=None,
        filename=None,
        radius_star=None,
        radius_in=0,
        radius_out=np.inf,
    ):

        if data is not None:
            if data.shape[1] == 4
                self.data = data
            elif data.shape[1] == 3:
                self.data = np.vstack((
                    data[:, 0:2],
                    np.zeros(data.shape[0]),
                    data[:, 2, None],
                ))

        elif filename is not None:
            self.data = self.load(filename)

        self.data_rotated = None
        self.radius_star = radius_star
        self.radius_in = radius_in
        self.radius_out = radius_out


    def load(
        filename,
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
                (np.linalg.norm(data[:, 0:2], axis=1) >= self.radius_in) *
                (np.linalg.norm(data[:, 0:2], axis=1) <= self.radius_out)
            )
            data = data[np.where(mask)]

        elif method == "ascii":
            data = []
            for line in infile:
                line = [float(value) for value in line.split(separator)]
                if len(line) >= 3:
                    if (self.radius_in <=
                        np.linalg.norm(line[0:2]) <=
                        self.radius_out
                    ):
                        data.append()
            data = np.array(data)

        infile.close()
        t_end = time.time()  # End of timer.
        print "Loading took %f seconds." % (t_end - t_start)

        if data.shape[1] < 4:
            # Add the z-dimension if it is not already there.
            z = np.zeros((data.shape[0], 1))
            data = np.hstack((data[:, 0:2], z, data[:, 2, None]))
        self.data = data


    def rotate(self, angle_z=0, angle_y=0, angle_x=None, unit="deg"):
        """Rotate entire dataset by an angle around any axis.
        TODO: Update this docstring.

        data: (float, array) The dataset to be rotated. Array of shape (N, 4),
            where N is the number of datapoints.
        angle_z: (float) Angle to rotate around z-axis. This is the rotational
            axis for the disk. It can be gradually increased to simulate the
            orbital motion of the system..
        angle_y: (float) Angle to rotate around y-axis. The inclination between
            the disk and the field of view. This angle should always be the same
            for one analysis if the disk is not wobbling.
        angle_x: The x-axis is the line of sight. Rotations around this axis
            would have no effect on the received flux, therefore this angle is
            ignored.
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
        angle_y *= factor
        angle_z *= factor

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
        rotation_matrix = R_y * R_z

        coords_in = self.data[:, :~0]
        coords_out = (rotation_matrix * coords_in.transpose()).transpose()
        self.data_rotated = np.hstack((coords_out, self.data[:, ~0, None]))


    def add_3d_points(self, H, n_layers=None, dz=None, thickness=None, ratio=None):
        """Expands a 2D-disk into 3rd (z) dimension assuming a simple model.
        TODO: Update this docstring.

        z-coordinate must be 0 for every point in data. If the dataset is already
        3-dimensional then you NEED not and MUST not use this method. You must
        use this method BEFORE any rotation around x- or y-axis, NOT AFTER.
        Rotations around z-axis does not matter, as they do not alter z-data.

        Note that this method multiplies the size of the dataset provided, so
        the returned data can be very large.


        H: (float) A physical parameter. We assume it to be a constant.
        n_layers: (int) How many data layers to add to EACH side of the disk.
            The total number of layers in z-direction becomes 2*n_layers+1.
        dz: (float) Distance between each layer.
        thickness: (float) Total thickness of disk to provide data for. This
            argument can be given instead of either n_layers or dz. If both
            n_layers and dz are provided then thickness is ignored.
        ratio: (float) Automatically choose a thickness so the outer layers have
            a density of ratio*density_z0. If more than 1 of n_layers, dz,
            thickness are provided then ratio is ignored.

        return: (float, array) The dataset with 2*n_layers*N more points added.
        """

        if thickness is None:
            thickness = - np.log(ratio) * H
        if n_layers is None:
            n_layers = int(round(0.5 * thickness / dz))
        elif dz is None:
            dz = 0.5 * thickness / n_layers

        N = self.data.shape[0]
        data_over = np.zeros((N*n_layers, 4))

        for k in xrange(1, n_layers+1):
            data_over[(k-1)*N : k*N, 0:2] = self.data[:, 0:2]
            data_over[(k-1)*N : k*N, 2] = k * dz
            data_over[(k-1)*N : k*N, 3] = self.data[:, 3] * np.exp(- k * dz / H)
        data_under = data_over.copy()
        data_under[:, 2] *= -1

        print (
            "%d layers added on each side of the disk. "
            "Size of dataset is increased from %d to %d points."
            % (n_layers, self.data.shape[0], (2*n_layers+1)*self.data.shape[0])
        )
        self.data =  np.vstack((self.data, data_over, data_under))


    def distance(self, p1, p2=None):
        """Returns the distances from a set of points to a line.
        TODO: Update this docstring.

        points: (float, array) Coordinates represented by an array of shape
            (N, d) where N is the number of points and d is the number of
            dimensions (2 or 3).
        p1: (float, array) A point in space which defines a line with p2.
        p2: (float, array) A point in space which defines a line with p1. If
            p2 is not provided it is assumed that the line is parallell to
            the x-axis.

        return: (float, array) The shortest euclidian distances between
            points and the line (p1, p2).
        """
        if p2 is None:
            p2 = p1.copy()
            p2[0] += 1.
        return np.linalg.norm(
            np.cross(p1 - self.data[:, :~0], p2 - self.data[:, :~0]) /
            np.linalg.norm(p2 - p1),
            axis=1,
        )


    def get_sylinder(self, radius_star):
        """Slice out a sylinder shape from a set of datapoints.
        TODO: Update this docstring.

        The sylinder is always centered on the x-axis and spans x=(0, inf].

        data: (float, array) The dataset to be sliced. Array of shape (N, 4),
            where N is the number of datapoints.
        radius_star: (float) Radius of the sylinder (a distance in the yz-plane).

        return: (float, array) The slice of the dataset contained in the sylinder.
        """
        mask = (
            (selfdata[:, 0] > 0) *
            (np.linalg.norm(self.data[:, 1:3], axis=1) <= radius_star)
        )
        data_sylinder = self.data[np.where(mask)]
        return data_sylinder


