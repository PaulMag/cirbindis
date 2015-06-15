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
from scipy import integrate
from scipy import special
import astropy.units as u
import textwrap
    # String manipulations.

from Star import Star
import Functions as func



class DensityMap:


    def __init__(self,
        data=None,
        filename=None,
        dataname=None,
        coordsystem="cartesian",
        outfolder=None,
        unit=None,
        inclinations=None,
        radius_in=0,
        radius_out=np.inf,
        diskmass=.01,
        diskradius=1000.,
        H=1.,
        kappa=10.
    ):

        self.data_rotated = None
        # If the inclination is a single number, put it in a list:
        try:
            iter(inclinations)
            self.inclinations = inclinations
        except TypeError:
            self.inclinations = [inclinations]
        if dataname is None or dataname == "":
            self.dataname = filename.split("/")[~0]
        else:
            self.dataname = dataname
        self.outfolder = outfolder
        self.unit = unit
        self.stars = []
        self.radius_in = radius_in
        self.radius_out = radius_out
        self.H = H
        self.kappa = kappa  # [cm^2 / g]
            # Between 5 and 100 according to Boubier et al. 1999.

        if data is not None:
            self.data = data
        elif filename is not None:
            self.load(filename)

        if coordsystem == "cartesian":
            pass
        elif coordsystem == "polar":
            x, y = func.pol2cart(self.data[:, 0], self.data[:, 1])
            self.data[:, 0], self.data[:, 1] = x, y
        else:
            raise KeyError("Coordinate system must be 'cartesian' or 'polar'.")

        self.set_physical_units(diskmass, diskradius)


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
            self.data = data[np.where(mask)]

        elif method == "ascii":
            data = []
            for line in infile:
                line = line.rstrip().split(separator)
                if len(line) >= 3:
                    line = [float(value) for value in line]
                    if (self.radius_in <=
                        np.linalg.norm(line[0:2]) <=
                        self.radius_out
                    ):
                        data.append(line)
            self.data = np.array(data)

        infile.close()
        t_end = time.time()  # End of timer.
        print "Loading took %f seconds." % (t_end - t_start)


    def set_physical_units(self, mass_total, r_total):
        """Convert density to physical units related to the total mass and
        size of the disk.
        """

        r_in = u.Quantity(self.radius_in, self.unit["distance"])
        r_out = u.Quantity(self.radius_out, self.unit["distance"])
        r_total = u.Quantity(r_total, self.unit["distance"])
        H = u.Quantity(self.H, self.unit["distance"])
        mass_total = u.Quantity(mass_total, self.unit["mass"])

        mass_central = (
            (r_out**0.5 - r_in**0.5) /
            (r_total**0.5 - r_in**0.5) *
            mass_total
        )
        rho_central = (
            mass_central / (2*np.pi * (r_out**2 - r_in**2) * 2*H)
        ).to(u.Unit(self.unit["mass"]) / u.Unit(self.unit["distance"])**3).value

        # Simple static alternative:
        # rho_central = (
            # mass_total.value / (
                # np.pi * u.Quantity(50, "AU").to(
                    # u.Unit(self.unit["distance"])
                # ).value**2 * 2*self.H
            # )
        # )

        # self.data[:, ~0] /= self.data[:, ~0].mean()
            # Normalize before scaling, or not?
        self.data[:, ~0] *= rho_central


    def writeto(self, filename, method=None, separator=" "):
        """Write self.data to a file for later use.

        Assumed to be on the form 'x,y,z,density' which represents a point in
        cartesian space and the density at that point.

        filename: (string) Full pathname to outfile for writing data.
        method: (string) What kind of writing algorithm to use. Recommended to
            use 'pickle' if it will be loaded by this program later (faster)
            and 'ascii' for any other purpose. If none is given, will try to
            automatically find out by looking at filename ending.
        separator: (string) If method='ascii' this is the separator between the
            values each line. Usually a space or comma. Ignored if
            method='pickle'.
        """

        t_start = time.time()
            # Just to time the writing, in case of large dataset.

        if method is None:
            if filename.endswith(".p") or filename.endswith(".pickle"):
                method = "pickle"
            else:
                method = "ascii"

        if method == "pickle":
            outfile = open(filename, "wb")
            pickle.dump(self.data, outfile)

        elif method == "ascii":
            outfile = open(filename, "w")
            for line in self.data:
                outfile.write("%f%s%f%s%f\n" % (
                    line[0], separator,
                    line[1], separator,
                    line[2],
                ))

        outfile.close()
        t_end = time.time()  # End of timer.
        print "Writing took %f seconds." % (t_end - t_start)


    def rotate(self, angle_z=0, unit="deg"):
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
        angle_z *= factor

        # Make rotation matrix:
        rotation_matrix = np.matrix([
            [ np.cos(angle_z), -np.sin(angle_z)],
            [ np.sin(angle_z),  np.cos(angle_z)],
        ])

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

        return np.abs(
            np.cross(
                p1 - self.data_rotated[:, :~0],
                p2 - self.data_rotated[:, :~0],
            ) /
            np.linalg.norm(p2 - p1),
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
        H=None,
        n_angle=None,
        dtheta=None,
        theta=None,
        unit="deg",
        n_radius=None,
        dr=None,
        save=False,
        show=False,
        outfolder=None,
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

        if H is None:
            H = self.H

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
                    unit=self.unit,
                    radius_in=self.radius_in,
                    radius_out=self.radius_out,
                    kappa=self.kappa
                )
                for j, inclination in enumerate(inclinations):
                    sylinder.space_sylinder(
                        inclination=inclination,
                        unit=unit,
                        H=H,
                        n_steps=n_radius,
                        dr=dr,
                    )
                    lightcurve[j, i] += sylinder.integrate()
        print "%f / %f" % (theta, theta)

        lightcurve /= lightcurve.mean(axis=1)[:, None]

        for j, inclination in enumerate(inclinations):

            starradius = ""
            starflux = ""
            for star in self.stars:
                starradius += "%g-" % star.radius
                starflux += "%g-" % star.intensity
            starradius = starradius.rstrip("-")
            starflux = starflux.rstrip("-")
            header = (
                "%s, H=%g, kappa=%g, "
                "r_star=%s, flux_star=%s, r_in=%g, r_out=%g, dr=%g, "
                "dtheta=%g%s, inc=%g%s"
                % ( self.dataname,
                    H,
                    self.kappa,
                    starradius,
                    starflux,
                    self.radius_in,
                    self.radius_out,
                    (self.radius_out - self.radius_in) / n_radius,
                    float(theta) / n_angle,
                    unit,
                    inclination,
                    unit,
                )
            )
            if save:
                outname = (
                    "%s__H=%g__"
                    "r_in=%g__r_out=%g__"
                    "inc=%02g"
                    % ( self.dataname,
                        H,
                        self.radius_in,
                        self.radius_out,
                        inclination,
                    )
                )
                if outfolder is None:
                    outfolder = self.outfolder
                func.make_folder(outfolder)
                outfile = open("%s/%s.csv" % (outfolder, outname), "w")
                outfile.write("#" + header + "\n")
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
            plt.title("\n".join(textwrap.wrap(header.split(", inc")[0], 70)))
            plt.xlabel("rotational angle [degree]")
            plt.ylabel("normalized intensity")
            plt.legend(loc="best")
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
            plt.tight_layout()
            plt.show()



class Sylinder(DensityMap):


    def __init__(self,
        star,
        data,
        unit=None,
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
        self.unit = unit
        self.radius_in = radius_in
        self.radius_out = radius_out
        self.kappa = kappa  # [cm^2 / g]

        self.data = data


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

        g = np.sqrt(2) * H  # Constant used several times in calculations.

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
                np.sum(
                    # \int_z1^z2 \rho_0 * e^{- z^2 / (2*H^2)} dz
                    g * data[start:end, ~0] * 0.5 * np.sqrt(np.pi) *
                    (special.erf(z2 / g) - special.erf(z1 / g))
                ) / (2. * np.sum(W))
            )

        self.densities = densities
        self.drs = drs


    def integrate(self):
        """Integrates the intensity through the layers of dust.

        Assumes that space_sylinder has just been called and used its results.

        return: (float) Perceived intensity outside the disk
        """

        kappa = self.kappa * u.Unit("cm2/gram").to(
            u.Unit(self.unit["distance"])**2 / u.Unit(self.unit["mass"])
        )
        intensity = self.star.intensity

        for density, dr in zip(self.densities, self.drs):
            tau = kappa * density * dr
            intensity *= np.exp(-tau)

        return intensity * (u.Unit(self.unit["intensity"])).to("erg / (cm2 s)")
