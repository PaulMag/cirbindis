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
import re
    # Regular expressions for data reading.
import matplotlib.pyplot as plt
from matplotlib import ticker
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
        H0=1.,
        R0=1.,
        H_power=1.,
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
        print (
            "Loading dataset '%s' from file '%s'..."
            % (self.dataname, filename)
        )
        self.outfolder = outfolder
        self.unit = unit
        self.stars = []
        self.radius_in = radius_in
        self.radius_out = radius_out
        self.diskmass = diskmass
        self.diskradius = diskradius
        self.H0 = H0
        self.R0 = R0
        self.H_power = H_power
        self.kappa = kappa  # [cm^2 / g]
            # Between 5 and 100 according to Bouvier et al. 1999.

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
                line = re.split("[, ]+", line.strip())
                    # Split on comma/space and strip spaces/newline.
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

        # Set a finite radius_out if there was no limit:
        if np.isinf(self.radius_out):
            self.radius_out = np.linalg.norm(self.data[:, 0:2], axis=1).max()


    def set_physical_units(self, mass_total=None, r_total=None):
        """Convert density to physical units related to the total mass and
        size of the disk.
        """

        if mass_total is None:
            mass_total = self.diskmass
        if r_total is None:
            r_total = self.diskradius

        r_in = u.Quantity(self.radius_in, self.unit["distance"])
        print "r_in =", r_in
        r_out = u.Quantity(self.radius_out, self.unit["distance"])
        print "r_out =", r_out
        r_total = u.Quantity(r_total, self.unit["distance"])
        print "r_total =", r_total
        H0 = u.Quantity(self.H0, self.unit["distance"])
        print "H0 =", H0
        mass_total = u.Quantity(mass_total, self.unit["mass"])
        print "mass_total =", mass_total

        mass_central = (
            (r_out**0.5 - r_in**0.5) /
            (r_total**0.5 - r_in**0.5) *
            mass_total
        )
        print "mass_central =", u.Quantity(mass_central, self.unit["mass"])
        rho_central = (
            mass_central / (np.pi * (r_out**2 - r_in**2) * np.sqrt(np.pi/2)*H0)
        ).to(u.Unit(self.unit["mass"]) / u.Unit(self.unit["distance"])**3).value
        print "rho_central =", u.Quantity(
            rho_central,
            u.Unit(self.unit["mass"]) / u.Unit(self.unit["distance"])**3,
        )
        rho_central_CGM = (
            mass_central / (np.pi * (r_out**2 - r_in**2) * np.sqrt(np.pi/2)*H0)
        ).to(u.Unit("g") / u.Unit("cm")**3).value
        print "rho_central =", u.Quantity(
            rho_central_CGM,
            u.Unit("g") / u.Unit("cm")**3,
        )

        self.data[:, ~0] /= self.data[:, ~0].mean()
            # Normalize before scaling, or not? Yes.
        self.data[:, ~0] *= rho_central

        r_innercavity = 2.0  # This variable is currently hardcoded.
        mask_innercavity = np.linalg.norm(self.data[:, 0:2], axis=1) <= r_innercavity
        data_innercavity = self.data[np.where(mask_innercavity)]
        rho_innercavity = data_innercavity[:, ~0].mean()
        mass_innercavity = rho_innercavity * (np.pi * (r_innercavity**2 - self.radius_in**2) * np.sqrt(np.pi/2)*self.H0)
        print "rho_innercavity r=[r_in,%g] =" % r_innercavity, u.Quantity(
            rho_innercavity,
            u.Unit(self.unit["mass"]) / u.Unit(self.unit["distance"])**3,
        )
        print "rho_innercavity r=[r_in,%g] =" % r_innercavity, u.Quantity(
            rho_innercavity,
            u.Unit(self.unit["mass"]) / u.Unit(self.unit["distance"])**3,
        ).to(u.Unit("g") / u.Unit("cm")**3)
        # ).to(u.Unit("g") / u.Unit("cm")**3).value)
        print "mass_innercavity r=[r_in,%g] =" % r_innercavity, u.Quantity(
            mass_innercavity,
            u.Unit(self.unit["mass"]),
        )


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
        angle_z *= + factor
            # Multiply by positive factor so that the disk is rotated
            # counter-clockwise. Then the line of sight (observer) is rotated
            # clockwise.

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
            (self.distance(star.position_rotated) <= star.radius)
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
        H0=None,
        R0=None,
        H_power=None,
        n_angle=None,
        dtheta=None,
        theta=None,
        unit="deg",
        n_radius=None,
        dr=None,
        lcurve_show=False,
        lcurve_savefig=False,
        lcurve_savecsv=False,
        dprofile_show=False,
        dprofile_savefig=False,
        normalizations=["stellar"],
        short_title=True,
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

        if H0 is None:
            H0 = self.H0
        if R0 is None:
            R0 = self.R0
        if H_power is None:
            H_power = self.H_power

        if n_angle is None:
            n_angle = int(round(float(theta) / dtheta))
        elif dtheta is None:
            dtheta = float(theta) / n_angle

        angles = np.linspace(0, theta-dtheta, n_angle)
        lightcurve = np.zeros((len(inclinations), n_angle))

        # Density profile:
        if dprofile_show or dprofile_savefig:
            fig_dprof = plt.figure(figsize=(12,6))
            fig_dprof.suptitle("%s" % self.dataname)
            axes_dprof = []
            nplots = [
                int(round(np.sqrt(len(inclinations)))),  # No of columns.
                int(np.ceil(np.sqrt(len(inclinations)))),  # No of rows.
                len(inclinations),  # Total no of density profile plots.
            ]
            radius_max = None
            plotcolors = ("b", "g", "r", "c", "m", "y")
            plotlinestyles = ('-', '--', '-.', ':')

        for i, angle in enumerate(angles):
            print "%6.2f / %g" % (angle, theta)
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
                        H0=H0,
                        R0=R0,
                        H_power=H_power,
                        n_steps=n_radius,
                        dr=dr,
                    )
                    lightcurve[j, i] += sylinder.integrate()

                    # Density profile:
                    if (dprofile_show or dprofile_savefig) and (k < 4):
                        # Only do for maximum 4 stars.
                        if (i == 0) and (k == 0):  # Only initialize axes once:
                            axes_dprof.append(fig_dprof.add_subplot(
                                nplots[1],
                                nplots[0],
                                j+1,
                            ))
                        if k == 0:
                            axes_dprof[j].plot(
                                sylinder.radiuses,
                                sylinder.densities *
                                    u.Unit(
                                        u.Unit(self.unit["mass"]) /
                                        u.Unit(self.unit["distance"])**3
                                    ).to("gram/cm3"),
                                "%s%s" % (plotcolors[i], plotlinestyles[k]),
                                label="%3g" % angle,
                            )
                        else:
                            axes_dprof[j].plot(
                                sylinder.radiuses,
                                sylinder.densities *
                                    u.Unit(
                                        u.Unit(self.unit["mass"]) /
                                        u.Unit(self.unit["distance"])**3
                                    ).to("gram/cm3"),
                                "%s%s" % (plotcolors[i], plotlinestyles[k]),
                            )
                        for l, density in enumerate(sylinder.densities):
                            if density == 0:
                                if sylinder.radiuses[l] > radius_max:
                                    radius_max = sylinder.radiuses[l]
                                break
                        axes_dprof[j].set_title("inc=%2g" % inclination)

        # Density profile:
        if dprofile_show or dprofile_savefig:
            if radius_max is None:
                radius_max = self.radius_out
            for j in range(nplots[2]):  # All subplots.
                axes_dprof[j].set_xlim([self.radius_in, radius_max])
                if False:  # Hardcoded switch.
                    axes_dprof[j].set_ylim([1e-30, 1e-12])
                try:
                    axes_dprof[j].set_yscale("log")  # Does not work if all 0s.
                except:
                    pass
                axes_dprof[j].yaxis.set_major_formatter( \
                    ticker.FormatStrFormatter('%.1e'))
            for j in range(nplots[0]):  # Bottom row.
                axes_dprof[~j].set_xlabel("radius [a]")
            for j in range(nplots[2] - nplots[0]):  # All except bottom row.
                axes_dprof[j].set_xticklabels([])
            for j in range(0, nplots[2], nplots[0]):  # Left coumn.
                axes_dprof[j].set_ylabel("density [g/cm^3]")
            axes_dprof[nplots[0]-1].legend(  # Only top right.
                title="v.angle [deg]=",
                loc="best",
            )

        print "%6.2f / %g" % (theta, theta)

        if "all" in normalizations:
            normalizations = ["stellar", "max", "mean"]

        for normalization in normalizations:

            if "stellar" in normalization or "unobscured" in normalization:
                unobscured_flux = 0.
                for star in self.stars:
                    unobscured_flux += star.intensity
                lightcurve /= unobscured_flux
            elif "max" in normalization:
                for j, maxflux in enumerate(lightcurve.max(axis=1)):
                    if maxflux > 0:
                        lightcurve[j] /= maxflux
                    # else lightcurve[j] is all zeros, so avoid dividing by 0
            elif "mean" in normalization:
                for j, maxflux in enumerate(lightcurve.max(axis=1)):
                    if maxflux > 0:
                        lightcurve[j] /= maxflux
                    # else lightcurve[j] is all zeros, so avoid dividing by 0

            if lcurve_show or lcurve_savefig:
                fig = plt.figure(figsize=(12,6))
                fig.gca().get_yaxis().get_major_formatter().set_useOffset(False)
                    # Always use absolute labels and not offsets.
                plt.minorticks_on()  # Turn default minorticks on for y-axis.
                ax = fig.add_subplot(1,1,1)
                ax.set_xlim([0, 360])
                ax.set_xticks(range(0, 360+1, 30))
                ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(3))
                    # A major tick every 30 deg and minor tick every 10 deg.

                # Set a specific tick step (hardcoded switch):
                if False:
                    step = 0.1
                    stepmin = lightcurve.min()
                    stepmax = lightcurve.max()
                    stepdiff = stepmax - stepmin
                    stepmin -= stepdiff  # Increase the range of the ticks just
                    stepmax += stepdiff  # to be sure.
                    stepmin = step * round(stepmin / step)  # Round off to the
                    stepmax = step * round(stepmax / step)  # nearest step.
                    ax.set_yticks(np.linspace(
                        stepmin,
                        stepmax,
                        int(round((stepmax - stepmin) / step)) + 1,
                    ))
                    # Set a specific label step (hardcoded switch):
                    if True:
                        # Automatic minor ticks between steps.
                        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
                    elif True:
                        # Manual label removal of certain steps.
                        steplabel = 2 * step
                        ylabels = []
                        for ytick in ax.get_yticks():
                            if (int(round(ytick/step)) %
                                int(round(steplabel/step)) == 0
                            ):
                                ylabels.append(ytick)
                            else:
                                ylabels.append("")
                        ax.set_yticklabels(ylabels)

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
                        H0,
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
                if lcurve_savecsv:
                    outname = (
                        "%s__H=%g__"
                        "r_in=%g__r_out=%g__"
                        "inc=%02g__%snorm"
                        % ( self.dataname,
                            H0,
                            self.radius_in,
                            self.radius_out,
                            inclination,
                            normalization,
                        )
                    )
                    if outfolder is None:
                        outfolder = self.outfolder
                    func.make_folder(outfolder)
                    func.make_folder(outfolder + "/csvtables")
                    outfile = open("%s/csvtables/%s.csv" \
                        % (outfolder, outname), "w")
                    outfile.write("#" + header + "\n")
                    for angle, flux in zip(angles, lightcurve[j]):
                        outfile.write("%f,%f\n" % (angle, flux))
                    outfile.close()

                if lcurve_show or lcurve_savefig:
                    ax.plot(
                        angles,
                        lightcurve[j],
                        label="%2g" % inclinations[j],
                    )

            if lcurve_show or lcurve_savefig:
                if short_title:
                    # Only use the name of the data in the title.
                    ax.set_title(self.dataname)
                    ax.set_position([0.10, 0.10, 0.80, 0.83])
                else:
                    # Use all metadata in the title.
                    ax.set_title(
                        "\n".join(textwrap.wrap(header.split(", inc")[0], 70))
                    )
                    ax.set_position([0.10, 0.10, 0.80, 0.80])
                ax.set_xlabel("viewing angle [degree]")
                ax.set_ylabel("flux [%s flux]" % normalization)
                ax.legend(
                    bbox_to_anchor=(1.11, 0.5),
                    title="inc [deg]=",
                    loc="right",
                    borderaxespad=0.,
                )

            if lcurve_savefig or dprofile_savefig:
                if outfolder is None:
                    outfolder = self.outfolder
                func.make_folder(outfolder)
                func.make_folder(outfolder + "/plots")
                outname = (
                    "%s__H=%g__"
                    "r_in=%g__r_out=%g__"
                    % ( self.dataname,
                        H0,
                        self.radius_in,
                        self.radius_out,
                    )
                )

            if lcurve_savefig:
                fig.savefig("%s/plots/%s%snorm.png" \
                    % (outfolder, outname, normalization))
            if dprofile_savefig:
                fig_dprof.savefig("%s/plots/%sdprofiles.png" \
                    % (outfolder, outname))

        if lcurve_show or dprofile_show:
            raw_input("Press <enter> to view the plots: ")
            if lcurve_show:
                fig.show()
            if dprofile_show:
                fig_dprof.show()
            raw_input("Press <enter> to close plots and exit: ")



class Sylinder(DensityMap):


    def __init__(self,
        star,
        data,
        unit=None,
        inclinations=None,
        radius_in=0,
        radius_out=np.inf,
        kappa=10.,
        ngrid=7,
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
        inclination=90.0,
        unit="deg",
        H0=0.1,
        R0=1.0,
        H_power=0,
        n_steps=None,
        dr=None,
        ngridz=12,
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
        H0: (float) Thickness of the disk. Necessary for integral.
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
        inclination = np.pi/2 - inclination
            # Convert from standard inclination definition to what is used
            # in these calculations.

        if ngridz is None:
            ngridz = int(max(min(10. * self.star.radius / H, 10), 2000) + 0.5)

        if n_steps is None:
            n_steps = int(round((self.radius_out - self.radius_in) / dr))
        dpoints = int(round(self.data.shape[0] / float(n_steps)))
            # How many datapoints to include in each bin.

        densitygrids = np.zeros((n_steps, ngridz))
        densities = np.zeros(n_steps)
        drs = np.zeros(n_steps)
        radiuses = np.zeros(n_steps)

        z1_grid = np.linspace(
            -self.star.radius,
            self.star.radius - 2. * self.star.radius / ngridz,
            ngridz,
        )
        z2_grid = z1_grid + 2. * self.star.radius / ngridz
        z1_grid /= np.cos(inclination)
        z2_grid /= np.cos(inclination)

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
            radiuses[i] = self.radius_in + np.sum(drs[:i]) + drs[i]/2.

            H = H0 * (radiuses[i] / R0)**H_power  # H in current bin.
            g = np.sqrt(2) * H  # Used several times in calculations.

            W = np.sqrt(
                self.star.radius**2 -
                (data[start:end, 1] - self.star.position_rotated[1])**2
            ) / np.cos(inclination)
                # In the docs sin(phi) is used since the opposite inclination
                # definition is explained.
            z = (
                (data[start:end, 0] - self.star.position_rotated[0]) *
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
            ) * H0 / H

            z1b = np.maximum(z[:,None] + z1_grid, z1[:,None])
            z2b = np.minimum(z[:,None] + z2_grid, z2[:,None])
            mask = z2b - z1b < 0
            z2b[mask] = z1b[mask]
            densitygrids[i, :] = (
                np.sum(
                    g * data[start:end, ~0, None] * 0.5 * np.sqrt(np.pi) *
                    (special.erf(z2b / g) - special.erf(z1b / g)),
                    axis=0,
                ) / (2. * np.sum(z2b - z1b, axis=0))
            )
        self.densitygrids = densitygrids
        self.densities = densities
        self.drs = drs
        self.radiuses = radiuses


    def integrate(self):
        """Integrates the flux through the layers of dust.

        Assumes that space_sylinder has just been called and used its results.

        return: (float) Perceived flux outside the disk
        """

        kappa = self.kappa * u.Unit("cm2/gram").to(
            u.Unit(self.unit["distance"])**2 / u.Unit(self.unit["mass"])
        )

        fluxgrid = np.zeros(self.densitygrids.shape[1])
        fluxR = np.sqrt(self.star.intensity / np.pi)  # The flux "radius".
        fluxRs = np.linspace(-fluxR, fluxR, fluxgrid.size+1)

        def fluxW(fluxr):
            return np.sqrt(fluxR**2 - fluxr**2)

        # print
        # print self.star.intensity
        # print fluxRs
        for i in xrange(fluxgrid.size):
            fluxgrid[i] = integrate.quad(fluxW, fluxRs[i], fluxRs[i+1])[0] * 2
        # print

        # print fluxgrid
        for densitygrid, dr in zip(self.densitygrids, self.drs):
            tau = kappa * densitygrid * dr
            fluxgrid *= np.exp(-tau)
            # print densitygrid
            # print fluxgrid

        return fluxgrid.sum()
