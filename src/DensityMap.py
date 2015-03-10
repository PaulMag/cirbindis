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

kappa = 1.



class DensityMap:


    def __init__(self,
        data=None,
        filename=None,
        inclinations=None,
        r0=1.,
        radius_star=None,
        radius_in=0,
        radius_out=np.inf,
    ):

        self.data_rotated = None
        # If the inclination is a single number, put it in a list:
        try:
            iter(inclinations)
            self.inclinations = inclinations
        except TypeError:
            self.inclinations = [inclinations]
        self.r0 = r0
        self.radius_star = radius_star
        self.radius_in = radius_in
        self.radius_out = radius_out

        if data is not None:
            if data.shape[1] == 4:
                self.data = data
            elif data.shape[1] == 3:
                self.data = np.vstack((
                    data[:, 0:2],
                    np.zeros(data.shape[0]),
                    data[:, 2, None],
                ))

        elif filename is not None:
            self.load(filename)


    def load(self,
        filename,
        method="pickle",
        separator=" ",
    ):
        """Load a dataset to analyse from a file.
        TODO: Update this docstring.

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

        t_start = time.time()
            # Just to time the loading, in case of large dataset.
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


    def writeto(self, filename, method="pickle", separator=" "):
        """Write a dataset to a file for later use.
        TODO: Update this docstring.

        Assumed to be on the form 'x,y,z,density' which represents a point in
        cartesian space and the density at that point.

        filename: (string) Full pathname to outfile for writing data.
        method: (string) What kind of writing algorithm to use. Recommended to
            use 'pickle' if it will be loaded by this program later (faster) and
            'ascii' for an other purpose.
        separator: (string) If method='ascii' this is the separator between the
            values each line. Usually a space or comma. Ignored if method='pickle'.
        """

        t_start = time.time()
            # Just to time the writing, in case of large dataset.

        if method == "pickle":
            outfile = open(filename, "wb")
            pickle.dump(self.data, outfile)

        elif method == "ascii":
            outfile = open(filename, "w")
            for line in self.data:
                outfile.write("%f%s%f%s%f%s%f\n" % (
                    line[0], separator,
                    line[1], separator,
                    line[2], separator,
                    line[3],
                ))

        outfile.close()
        t_end = time.time()  # End of timer.
        print "Writing took %f seconds." % (t_end - t_start)


    def set_H(self, H0, H_power):
        self.H0 = H0
        self.H_power = H_power

    def get_H(self, r):
        return self.H0 * (r / self.r0)**self.H_power

    def set_sigma(self, sigma0, sigma_power):
        self.sigma0 = sigma0
        self.sigma_power = sigma_power

    def get_sigma(self, r):
        return self.sigma0 * (r / self.r0)**self.sigma_power


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
        coords_out = \
            np.asarray(rotation_matrix * coords_in.transpose()).transpose()
        self.data_rotated = np.hstack((coords_out, self.data[:, ~0, None]))


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


    def get_sylinder(self):
        """Slice out a sylinder shape from a set of datapoints.
        TODO: Update this docstring.

        The sylinder is always centered on the x-axis and spans x=(0, inf].

        data: (float, array) The dataset to be sliced. Array of shape (N, 4),
            where N is the number of datapoints.
        radius_star: (float) Radius of the sylinder (a distance in the yz-plane).

        return: (float, array) The slice of the dataset contained in the sylinder.
        """
        try:
            mask = (
                (self.data_rotated[:, 0] > 0) *
                (   np.linalg.norm(self.data_rotated[:, 1:3], axis=1) <=
                    self.radius_star
                )
            )
            data_sylinder = self.data_rotated[np.where(mask)]
        except:
            mask = (
                (self.data[:, 0] > 0) *
                (   np.linalg.norm(self.data[:, 1:3], axis=1) <=
                    self.radius_star
                )
            )
            data_sylinder = self.data[np.where(mask)]
        return data_sylinder


    def make_lightcurve(self,
        inclinations=None,
        H=1.,
        n_angle=None,
        dtheta=None,
        theta=None,
        unit="deg",
        n_radius=None,
        dr=None,
        save=False,
        show=False,
    ):
        """Makes a lightcurve by calling the other methods for each orientation
        of the dataset. Sort of a main method.

        TODO: Complete docstring.
        """

        # Default setting for theta is a full revolution:
        if theta is None:
            if unit == "rad":
                theta = 2*np.pi
            elif unit == "deg":
                theta = 360.
            elif unit == "arcmin":
                theta = 360. * 60
            elif unit == "arcsec":
                theta = 360. * 3600

        if inclinations is None:
            inclinations = self.inclinations
        # If the inclination is a single number, put it in a list:
        try:
            iter(inclinations)
        except TypeError:
            inclinations = [inclinations]

        if n_angle is None:
            n_angle = int(round(float(theta) / dtheta))
        elif dtheta is None:
            dtheta = float(theta) / n_angle

        angles = np.linspace(0, theta-dtheta, n_angle)
        lightcurve = np.zeros((len(inclinations), n_angle))

        for i, angle in enumerate(angles):
            print "%f / %f" % (angle, theta)
            self.rotate(
                angle_z=angle,
                unit=unit,
            )
            sylinder = DensityMap(
                self.get_sylinder(),
                radius_star=self.radius_star,
                radius_in=self.radius_in,
                radius_out=self.radius_out,
            )
            sylinder.data = sylinder.get_sylinder()
            for j, inclination in enumerate(inclinations):
                densities, drs = space_sylinder(
                    sylinder.data,
                    self.radius_star,
                    inclination=inclination,
                    n_steps=n_radius,
                    dr=dr,
                    radius_in=self.radius_in,
                    radius_out=self.radius_out,
                )
                lightcurve[j, i] = integrate(densities, drs)
        print "%f / %f" % (theta, theta)

        lightcurve /= lightcurve.mean(axis=1)[:, None]

        for j, inclination in enumerate(inclinations):

            title = (
                "H=%g, kappa=%g, "
                "r_star=%g, r_in=%g, r_out=%g, dr=%g, "
                "dtheta=%g%s, inclination=%g%s"
                % ( H,
                    kappa,
                    self.radius_star,
                    self.radius_in,
                    self.radius_out,
                    (self.radius_out - self.radius_in) / n_radius,
                    float(theta) / n_angle,
                    unit,
                    inclination,
                    unit,
                )
            )
            xlabel = "rotation angle [%s]" % unit
            ylabel = "observed flux"

            if save:
                outfile = open("../results/%s.csv" % title.replace(", ", "__"), "w")
                outfile.write(title + "\n")
                outfile.write(xlabel + "\n")
                outfile.write(ylabel + "\n")
                for angle, flux in zip(angles, lightcurve[j]):
                    outfile.write("%f,%f\n" % (angle, flux))
                outfile.close()

            if show:
                plt.plot(
                    angles,
                    lightcurve[j],
                    label="inc=%g" % inclinations[j],
                )

        if show:
            plt.title(title)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.legend()
            plt.show()



def space_sylinder(
    data,
    radius_star,
    inclination=0,
    unit="deg",
    H=1.,
    n_steps=None,
    dr=None,
    radius_in=None,
    radius_out=None,
):
    """Bin a (sylinder shaped) set of datapoints into a set of mean densities.

    The sylinder is first sorted along the x-axis and is then cut along the
    x-axis like a loaf of bread. The mean density is then computed from each
    slice of the sylinder/bread.

    This method is the one using most time in this program.

    TODO: Decide on how to organize the other arguments.
    data: (float, array) The dataset to be binned. Array of shape (N, 4),
        where N is the number of datapoints.

    return: (float, array), (float, array) A list of mean densities and the
        corresponding list of delta radiuses for each bin. Both are arrays of
        length n_step. The arrays are order FROM inside of disk TO oustide
        of disk.
    """

    if unit == "rad":
        factor = 1.
    elif unit == "deg":
        factor = np.pi / 180.
    elif unit == "arcmin":
        factor = np.pi / 180. * 60
    elif unit == "arcsec":
        factor = np.pi / 180. * 3600
    inclination *= factor

    if n_steps is None:
        n_steps = int(round((radius_out-radius_in) / dr))
    dpoints = int(round(data.shape[0] / float(n_steps)))
        # How many datapoints to include in each bin.

    densities = np.zeros(n_steps)
    drs = np.zeros(n_steps)

    data = data[np.argsort(data[:, 0])]

    y0 = 0
    for i in xrange(n_steps):
        start = i*dpoints
        if i == n_steps-1:
            # If it is the last step, make sure the last few points are
            # included (in case there are some rounding problems).
            end = data.shape[0]
            drs[i] = data[end-1, 0] - data[start, 0]
            s = data[start:end].shape[0]
            drs[i] *= (s + 1.) / s
        else:
            end = (i+1)*dpoints
            drs[i] = data[end, 0] - data[start, 0]
        W = np.sqrt(
            radius_star**2 -
            (data[start:end, 1] - y0)**2
        ) / np.cos(inclination)
        z = (
            data[start:end, 0] * np.tan(inclination)
        )
        z1 = z - W
        z2 = z + W
        densities[i] = (
            data[start:end, 3] *
            H *
            (np.exp(- z1 / H) - np.exp(- z2 / H))
        ).sum() / (2 * np.sum(W))

    return densities, drs


def integrate(densities, drs):
    """Integrates the intensity through the layers of dust.

    densities: (float, array) List of (mean) densities through a field of
        view of the disk, ordered from inside to outside.
    drs: (float, array) List of the corresponding dr to each density
        measurement.

    return: (float) Perceived intensity outside the disk
    """

    intensity = 1.  # Or whatever the full intensity of the star is.

    for density, dr in zip(densities, drs):
        tau = kappa * density * dr
        intensity *= np.exp(-tau)

    return intensity
