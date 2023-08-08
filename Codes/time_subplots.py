# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 21:33:46 2023

@author: abd__
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import h5py

path = "D:/abd__/Downloads D/Bflip/Bflip/run0034/out"
file_list = os.listdir(path)

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
z = np.linspace(zmin, zmax, Nz)
    
# Create a figure with rows=2 and columns=4 (adjust as needed)
num_rows, num_cols = 2, 4
fig, axs = plt.subplots(num_rows, num_cols, figsize=(10,10), sharex=True, sharey=True)
        
# Call the plot_bow_shock function for each time slice and place the colormesh plots in the subplots
X, Y = np.meshgrid(x, z, indexing='ij')

times = np.linspace(0,50,8).reshape([2,4])
for i in range(num_rows):
    for j in range(num_cols):
        time_index = int(times[i, j])  # Calculate the correct time index
        filename = 'D:/abd__/Downloads D/Bflip/Bflip/run0034/out/'+file_list[time_index]
        f1 = h5py.File(filename,'r+')
        rho = gaussian_filter(f1['rho01'][:,int(Ny/2),:], 1)
        im = axs[i, j].pcolormesh(X, Y, rho, cmap='jet')
        axs[i, j].text(-0.4e7, 2.4e7, f'H+ Density\nTime: {time_index*160/50:.2f} sec',
                       bbox=dict(facecolor='white', edgecolor='black'
                                 , boxstyle='round,pad=0.3'),
                       size = 10)
    
# Adjust the spacing between subplots
plt.tight_layout()

# Display the plot
plt.show()