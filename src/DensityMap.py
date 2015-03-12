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

from Star import Star
import Functions as func



class DensityMap:


    def __init__(self,
        data=None,
        filename=None,
        inclinations=None,
        radius_in=0,
        radius_out=np.inf,
        kappa=10.
    ):

        self.data_rotated = None
        # If the inclination is a single number, put it in a list:
        try:
            iter(inclinations)
            self.inclinations = inclinations
        except TypeError:
            self.inclinations = [inclinations]
        self.stars = []
        self.radius_in = radius_in
        self.radius_out = radius_out
        self.kappa = kappa  # [cm^2 / g]
            # Between 5 and 100 according to Boubier et al. 1999.

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


    def add_star(self, d=None, position=None, radius=None, intensity=None):
        """Make a Star instance and store it in a list of stars for the disk.

        d: (dictionairy) Must contain position, radius and intensity and can
            be provided instead of giving these other arguments individually.
        position: (float, array-like) Coordinates of the star.
        radius: (float) Radius of the star.
        intensity: (float) Intensity of the star.
        """
        self.stars.append(Star(
            d=d, position=position, radius=radius, intensity=intensity
        ))


    def load(self,
        filename,
        method=None,
        separator=" ",
    ):
        """Load a dataset to analyse from a file.

        Assumed to be on the form 'x,y,density' or 'x,y,z,density' which
        represents a point in cartesian space and the density at that point.
        If the z-coordinate is not given it will be assumed to be 0 for
        every point. Resulting data is an array of shape (N, 4), where N is
        the number of data points.

        filename: (string) Full pathname to file containing dataset.
        method: (string) What kind of loading algorithm to use. Can be
            'ascii' or 'pickle', If none is given, will try to automatically
            find out by looking at file ending.
        separator: (string) If method='ascii' this is the separator
            between the values each line. Usually a space or comma. Ignored if
            method='pickle'.
        """

        t_start = time.time()
            # Just to time the loading, in case of large dataset.
        infile = open(filename, "r")

        if method is None:
            if filename.endswith(".p") or filename.endswith(".pickle"):
                method = "pickle"
            else:
                method = "ascii"

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
        """Write self.data to a file for later use.

        Assumed to be on the form 'x,y,z,density' which represents a point in
        cartesian space and the density at that point.

        filename: (string) Full pathname to outfile for writing data.
        method: (string) What kind of writing algorithm to use. Recommended to
            use 'pickle' if it will be loaded by this program later (faster) and
            'ascii' for any other purpose.
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


    def set_r0(self, r0=1.49597871e13):
        self.r0=r0  # [centimeters]

    def set_H(self, H0=0.03, H_power=1/4.):
        self.H0 = H0
        self.H_power = H_power

    def get_H(self, r):
        return r * self.H0 * (r / self.r0)**self.H_power  # [centimeters]

    def set_sigma(self, sigma0=1700., sigma_power=-3/2.):
        self.sigma0 = sigma0
        self.sigma_power = sigma_power

    def get_sigma(self, r):
        return self.sigma0 * (r / self.r0)**self.sigma_power  # [g / cm^2]


    def rotate(self, angle_z=0, angle_y=0, angle_x=None, unit="deg"):
        """Rotate entire dataset by an angle around any axis.

        The original data is not changed. Rotated version of data stored in
        self.data_rotated. self.stars are also rotated.

        angle_z: (float) Angle to rotate around z-axis. This is the rotational
            axis for the disk. It can be gradually increased to simulate the
            orbital rotation of the system.
        angle_y: (float) Angle to rotate around y-axis. The inclination between
            the disk and the field of view. This angle should always be the same
            for one analysis if the disk is not wobbling.
        angle_x: The x-axis is the line of sight. Rotations around this axis
            would have no effect on the received flux, therefore this angle is
            ignored.
        unit: (string) What unit angles are given in.
            'rad', 'deg', 'arcmin' or 'arcsec'.
        """

        # Transform angles into radians:
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

        # Make rotation matrix:
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

        # Rotate the disk:
        coords_in = self.data[:, :~0]
        coords_out = \
            np.asarray(rotation_matrix * coords_in.transpose()).transpose()
        self.data_rotated = np.hstack((coords_out, self.data[:, ~0, None]))

        # Rotate the stars:
        for star in self.stars:
            star.position_rotated = np.asarray(
                rotation_matrix * star.position[:, None]
            ).transpose()[0]


    def distance(self, p1, p2=None):
        """Returns the distances from the (rotated) datapoints to a line.

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
            np.cross(
                p1 - self.data_rotated[:, :~0],
                p2 - self.data_rotated[:, :~0],
            ) /
            np.linalg.norm(p2 - p1),
            axis=1,
        )


    def get_sylinder(self, starno=None, star=None):
        """Slice out a sylinder shape from the data based on a star.

        The sylinder is always oriented along the x-axis. Its size and
        position is determined by a star. It spans from the position of the
        surface of the star until x=inf.

        starno: (int) Index to get star from self.stars. Ignored if star is
            given.
        star: (Star instance) A star to base the sylinder on.

        return: (float, array) The slice of the dataset contained in the
        sylinder.
        """

        if star is None:
            star = self.stars[starno]

        mask = (
            (self.data_rotated[:, 0] > 0) *
            (self.distance(star.position) <= star.radius)
        )
        data_sylinder = self.data_rotated[np.where(mask)]
        return data_sylinder


    def get_density_profile(self, sigma=1, skip=2000, show=True):
        """Returns and displays the density profile of self.data.
        TODO: Finish this docstring.
        """

        from scipy.ndimage.filters import gaussian_filter1d

        radiuses = np.linalg.norm(self.data[:, 0:2], axis=1)
        indices_sorted = np.argsort(radiuses)
        radiuses = radiuses[indices_sorted][::skip]
        densities = gaussian_filter1d(
            self.data[:, 3][indices_sorted],
            sigma=sigma,
            mode="nearest",
        )[::skip]

        if show:
            plt.plot(radiuses, densities, "b+")
            plt.xlabel("radius")
            plt.ylabel("density")
            plt.show()

        return radiuses, densities


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
        inclinations = func.to_list(inclinations)

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
            for k, star in enumerate(self.stars):
                sylinder = Sylinder(
                    star=star,
                    data=self.get_sylinder(star=star),
                    radius_in=self.radius_in,
                    radius_out=self.radius_out,
                    kappa=self.kappa
                )
                for j, inclination in enumerate(inclinations):
                    sylinder.space_sylinder(
                        inclination=inclination,
                        unit=unit,
                        n_steps=n_radius,
                        dr=dr,
                    )
                    lightcurve[j, i] += sylinder.integrate()
        print "%f / %f" % (theta, theta)

        lightcurve /= lightcurve.mean(axis=1)[:, None]

        for j, inclination in enumerate(inclinations):

            title = (
                "H=%g, kappa=%g, "
                "r_star=%g, r_in=%g, r_out=%g, dr=%g, "
                "dtheta=%g%s, inclination=%g%s"
                % ( H,
                    self.kappa,
                    #TODO One for each star.
                    self.stars[0].radius,
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



class Sylinder(DensityMap):


    def __init__(self,
        star,
        data,
        inclinations=None,
        radius_in=0,
        radius_out=np.inf,
        kappa=10.
    ):

        # If the inclination is a single number, put it in a list:
        try:
            iter(inclinations)
            self.inclinations = inclinations
        except TypeError:
            self.inclinations = [inclinations]
        self.star = star
        self.radius_in = radius_in
        self.radius_out = radius_out
        self.kappa = kappa  # [cm^2 / g]


        if data.shape[1] == 4:
            self.data = data
        elif data.shape[1] == 3:
            self.data = np.vstack((
                data[:, 0:2],
                np.zeros(data.shape[0]),
                data[:, 2, None],
            ))


    def space_sylinder(self,
        inclination=0,
        unit="deg",
        H=1.,
        n_steps=None,
        dr=None,
    ):
        """Bin this sylinder's datapoints into a set of mean densities.

        The sylinder is first sorted along the x-axis and is then cut along
        the x-axis like a loaf of bread. Each point is integrated
        analytically through its projected density from the bottom to the
        top of the sylinder. The mean density is then computed from each
        slice of the sylinder/bread.

        This is stored temporarily as self.densities and self.drs: (float,
        array), (float, array) A list of mean densities and the
        corresponding list of delta radiuses for each bin. Both are arrays
        of length n_step. The arrays are order FROM inside of disk TO
        oustide of disk.

        inclination: (float) The angle to incline the line of sight on the
            sylinder.
        deg: (string) Unit of the angle.
        H: (float) Thickness of the disk. Necessary for integral.
        n_steps: (int) How many slices to divide the sylinder in. Affects
            accuracy of integral.
        dr: (float) The width of each sylinder section. Ignored if n_steps
            is provided.

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
            n_steps = int(round((self.radius_out - self.radius_in) / dr))
        dpoints = int(round(self.data.shape[0] / float(n_steps)))
            # How many datapoints to include in each bin.

        densities = np.zeros(n_steps)
        drs = np.zeros(n_steps)

        data = self.data[np.argsort(self.data[:, 0])]

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
                self.star.radius**2 -
                (data[start:end, 1] - self.star.position[1])**2
            ) / np.cos(inclination)
            z = (
                (data[start:end, 0] - self.star.position[0]) *
                np.tan(inclination)
            )
            z1 = z - W
            z2 = z + W
            densities[i] = (
                data[start:end, 3] *
                H *
                (np.exp(- z1 / H) - np.exp(- z2 / H))
            ).sum() / (2 * np.sum(W))

        self.densities = densities
        self.drs = drs


    def integrate(self):
        """Integrates the intensity through the layers of dust.

        Assumes that space_sylinder has just been called and used its results.

        return: (float) Perceived intensity outside the disk
        """

        intensity = self.star.intensity

        for density, dr in zip(self.densities, self.drs):
            tau = self.kappa * density * dr
            intensity *= np.exp(-tau)

        return intensity
