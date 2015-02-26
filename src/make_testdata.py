"""This script can be used to make an artificial density map with a smooth
sinusoidal shape. This can be used to test the analysis program.
"""

import numpy as np
import cPickle as pickle


Nx = int(6e2 + 1)
Ny = int(6e2 + 1)
strip_x = np.linspace(-4, 4, Nx)
strip_y = np.linspace(-4, 4, Ny)

x = np.zeros(Nx*Ny)
y = np.zeros(Nx*Ny)
z = np.zeros(Nx*Ny)

for i in xrange(Ny):
    x[i*Nx : (i+1)*Nx] = strip_x
for i in xrange(Nx):
    y[i :: Nx] = strip_y

def density(x, y):
    return (1. + np.sin(2*np.arctan2(y, x)))

def add_density(x, y):
    return np.vstack((x, y, z, density(x, y))).transpose()

data = add_density(x, y)

outfile = open("../data/testdata.p", "wb")
pickle.dump(data, outfile)
outfile.close()
