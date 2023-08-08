# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 14:54:29 2023

@author: abd__
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

variables = ['B', 'E', 'J', 'ji']
coordinates = ['x', 'y', 'z']
species = ['01', '02', '03']
molecules = ['H+', 'He++', 'H2O+']
planes = ['xy', 'xz', 'yz']

# Define the dimensions of the space
Nx = 180
Ny = 200
Nz = 430

# Define the limits of the space
xmin = -1.0e7   # [m]
xmax = 8.0e6    # [m]
ymin = -1.0e7   # [m]
ymax = 1.0e7    # [m]
zmin = -1.5e7   # [m]
zmax = 2.8e7    # [m]

# Generate the grid points in each direction
x = np.linspace(xmin, xmax, Nx)
y = np.linspace(ymin, ymax, Ny)

filename = 'D:/abd__/Downloads D/Bflip/Bflip/run0034/out/Amitis_comet_017200.h5'

f1 = h5py.File(filename,'r+')

# Particle Grid
X, Y = np.meshgrid(x, y, indexing='ij')

rho = f1['rho01'][:,:,0]
#rho = gaussian_filter(rho, 1)

plt.pcolormesh(X, Y, rho, shading='auto', cmap='jet')
plt.colorbar(label='Density', location='left')

# Define the subsampling factor
n = 10

# Generate the coordinates for all the grid points
X, Y = np.meshgrid(x[::n], y[::n], indexing='ij')
    
# Generate some example electric field data in the xy plane
vec_x = f1['Bx'][::n,::n,0]
vec_y = f1['By'][::n,::n,0]

# Calculate the magnitude of the electric field vectors
magnitude = np.sqrt(vec_x**2 + vec_y**2)

# Normalize the vectors to have the same length
normalized_x = vec_x / magnitude
normalized_y = vec_y / magnitude

# Plot the vector field in the xy plane with scaled vectors and colored based on magnitude
#plt.quiver(X, Y, normalized_x, normalized_y, ec = 'black'
           #, color='w', scale=20, lw = 0.1, width = 0.008)
plt.text(0.02e7,0.75e7,'No Gaussian KDE', color='black', fontsize=12)
plt.show()