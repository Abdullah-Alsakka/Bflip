# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 09:45:08 2023

@author: abd__
"""

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.ndimage import gaussian_filter
import h5py
import os

# Import file
filename = 'D:/abd__/Downloads D/Bflip/Bflip/run0034/out/Amitis_comet_006400.h5'
f1 = h5py.File(filename, 'r+')

# Generate some example data (replace this with your actual data)
# Let's assume you have arrays x, y, z representing the coordinates of the dots
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
x  = f1['rho01'][:,0,0]
y = f1['rho01'][0,0,:]
z = f1['rho01'][0,:,0]

# Generate an array of densities corresponding to each dot (replace this with your actual density data)
density = np.random.rand(100)

# Define the density
density = z
#f1['rho01']

# Normalize the density values to use them as sizes and colors for the dots
norm_density = (density - np.min(density)) / (np.max(density) - np.min(density))

# Create the 3D figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Define the colormap for coloring the dots based on density
color_map = plt.cm.get_cmap('jet')

# Create the scatter plot
scatter = ax.scatter(x, y, z, c=density, s=100*norm_density, cmap=color_map)

# Add a colorbar to the plot to show the density scale
cbar = fig.colorbar(scatter)
cbar.set_label('Density')

# Set axis labels (you can customize them as needed)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Show the plot
plt.show()
