import sys
import numpy as np
import matplotlib.pyplot as plt


def load(filename, separator=","):
    """Load a dataset to plot from a file.

    Assumed to be on the form 'angle,flux'.
    The first 3 lines are the title, xlabel and ylabel.

    filename: (string) Full pathname to file containing dataset.
    separator: (string) This is the separator between the values each line.
        Usually a space or comma.
    """

    infile = open(filename, "r")

    title = infile.readline()
    xlabel = infile.readline()
    ylabel = infile.readline()

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

    load(filename=sys.argv[1])
