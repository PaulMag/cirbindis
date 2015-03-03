import sys
import numpy as np
import matplotlib.pyplot as plt


def load(filename, separator=","):
    """Load a dataset to plot from a file.

    Assumed to be on the form 'angle,flux'.
    The first 3 lines are the title, xlabel and ylabel.

    filename: (string) or (string, list) Full pathname to file containing
        dataset or a list of pathnames.
    separator: (string) This is the separator between the values each line.
        Usually a space or comma.
    """

    if isinstance(filename, basestring):
        filenames = [filename]
    else:
        try:
            iter(filename)
            filenames = filename
        except TypeError:
            raise TypeError(
                "'filename' must be a string or sequence of strings."
            )


    infiles = [open(filename, "r") for filename in filenames]

    titles = [infile.readline() for infile in infiles]
    xlabels = [infile.readline() for infile in infiles]
    ylabels = [infile.readline() for infile in infiles]
    title = titles[0]
    xlabel = xlabels[0]
    ylabel = ylabels[0]

    for infile in infiles:
        data = []
        for line in infile:
            line = [float(value) for value in line.split(separator)]
            data.append(line)
        angles, fluxes = np.array(data).T
        infile.close()
        plt.plot(angles, fluxes)

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()


if __name__ == "__main__":

    load(filename=sys.argv[1:])
