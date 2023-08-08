# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 10:27:47 2023

@author: abd__
"""

"""
>>> CUSTOMIZE <<<
"""

plane = 'xz'   # Choose from xy, yz, xz
var = 'B'   # Choose from B, E, J, ji
species = '01'   # Choose from 01, 02, 03 corresponding to H+, He++, H20+

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import h5py
import os

path = "D:/abd__/Downloads D/Bflip/Bflip/run0034/out"
file_list = sorted(os.listdir(path))  # Sort the file list

# Define the dimensions of the space
Nx = 180
Ny = 200
Nz = 430

molecules = ['H+', 'He++', 'H2O+']

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

x_axis = locals()[plane[0]]
y_axis = locals()[plane[1]]

filename = 'D:/abd__/Downloads D/Bflip/Bflip/run0034/out/Amitis_comet_008400.h5'

f1 = h5py.File(filename,'r+')

plt.figure(figsize=(3.25,10))

# Particle Grid
X, Y = np.meshgrid(locals()[plane[0]], locals()[plane[1]], indexing='ij')

if plane == 'xy':
    rho = f1['rho'+species][:,:,int(Nz/2)]    
elif plane == 'xz':
    rho = f1['rho'+species][:,int(Ny/2),:]
elif plane == 'yz':
    rho = f1['rho'+species][int(Nx/2),:,:]
    
rho = gaussian_filter(rho, 1)

plt.pcolormesh(X, Y, rho, shading='auto', cmap='jet')
cbar = plt.colorbar(location='top')
cbar.set_label('$n_{H^+} [m^{-3}]$', fontsize=15)
cbar.ax.tick_params(labelsize=15)
'''
# Define the subsampling factor
n = 10

# Generate the coordinates for all the grid points
X, Y = np.meshgrid(locals()[plane[0]][::n], locals()[plane[1]][::n], indexing='ij')
    
# Generate some example electric field data in the xy plane

if plane == 'xy':
    if var == 'ji':
        vec_x = f1[var+plane[0]+species][::n,::n,int(Nz/2)]
        vec_y = f1[var+plane[1]+species][::n,::n,int(Nz/2)]                    
    else:
        vec_x = f1[var+plane[0]][::n,::n,int(Nz/2)]
        vec_y = f1[var+plane[1]][::n,::n,int(Nz/2)]
elif plane == 'xz':
    if var == 'ji':
        vec_x = f1[var+plane[0]+species][::n,int(Ny/2),::n]
        vec_y = f1[var+plane[1]+species][::n,int(Ny/2),::n]                    
    else:
        vec_x = f1[var+plane[0]][::n,int(Ny/2),::n]
        vec_y = f1[var+plane[1]][::n,int(Ny/2),::n]                    
elif plane == 'yz':
    if var == 'ji':
        vec_x = f1[var+plane[0]+species][int(Nx/2),::n,::n]
        vec_y = f1[var+plane[1]+species][int(Nx/2),::n,::n]                    
    else:
        vec_x = f1[var+plane[0]][int(Nx/2),::n,::n]
        vec_y = f1[var+plane[1]][int(Nx/2),::n,::n]                    

# Calculate the magnitude of the electric field vectors
magnitude = np.sqrt(vec_x**2 + vec_y**2)

# Normalize the vectors to have the same length
normalized_x = vec_x / magnitude
normalized_y = vec_y / magnitude

# Plot the vector field in the xy plane with scaled vectors and colored based on magnitude
plt.quiver(X, Y, normalized_x, normalized_y, ec = 'black'
           , color='w', scale=20, lw = 0.1, width = 0.006)
'''
V = np.array([[0,1]])
origin = np.array([[0.5e7],[-1e7]])
plt.quiver(*origin , V[:,0], V[:,1], color = 'w', scale = 6, width = 0.015)
plt.text(+0.58e7, -0.8e7, '$E$', color='w', size = 20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.axis('equal')
plt.xlim(xmin,xmax)
plt.ylim(zmin,zmax)
plt.xlabel('$x$ [m]', fontsize = 15)
plt.ylabel('$z$ [m]', fontsize = 15)
#plt.title(f'{var} Field in the {plane} Plane for {molecules[species.index(species)]}'
#         ,fontsize=20)

# plt.savefig('D:/abd__/Downloads D/Bflip/Reports/Report 1/Figure_1.png')
plt.show()