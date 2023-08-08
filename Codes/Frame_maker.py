# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 01:23:01 2023

@author: Abdullah Alsakka
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import os

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

for plane in planes:
    for sp in species:
        for var in variables:
            os.mkdir(f'D:/abd__/Downloads D/Bflip/Plots/{var}_{plane}_{sp}')
            print(f'We are looking at the vector field of {var} at species {sp} in plane {plane}')
            
            for i in range(len(file_list)):
                filename = 'D:/abd__/Downloads D/Bflip/Bflip/run0034/out/'+file_list[i]
                
                f1 = h5py.File(filename,'r+')
                
                # Particle Grid
                X, Y = np.meshgrid(locals()[plane[0]], locals()[plane[1]], indexing='ij')
            
                if plane == 'xy':
                    rho = f1['rho'+sp][:,:,int(Nz/2)]    
                elif plane == 'xz':
                    rho = f1['rho'+sp][:,int(Ny/2),:]
                elif plane == 'yz':
                    rho = f1['rho'+sp][int(Nx/2),:,:]
                
                rho = gaussian_filter(rho, 1)
            
                plt.pcolormesh(X, Y, rho, shading='auto', cmap='jet')
                plt.colorbar(label='Density', location='left')
                
                # Define the subsampling factor
                n = 10
                
                # Generate the coordinates for all the grid points
                X, Y = np.meshgrid(locals()[plane[0]][::n], locals()[plane[1]][::n], indexing='ij')
                    
                # Generate some example electric field data in the xy plane
                
                if plane == 'xy':
                    if var == 'ji':
                        vec_x = f1[var+plane[0]+sp][::n,::n,int(Nz/2)]
                        vec_y = f1[var+plane[1]+sp][::n,::n,int(Nz/2)]                    
                    else:
                        vec_x = f1[var+plane[0]][::n,::n,int(Nz/2)]
                        vec_y = f1[var+plane[1]][::n,::n,int(Nz/2)]
                elif plane == 'xz':
                    if var == 'ji':
                        vec_x = f1[var+plane[0]+sp][::n,int(Ny/2),::n]
                        vec_y = f1[var+plane[1]+sp][::n,int(Ny/2),::n]                    
                    else:
                        vec_x = f1[var+plane[0]][::n,int(Ny/2),::n]
                        vec_y = f1[var+plane[1]][::n,int(Ny/2),::n]                    
                elif plane == 'yz':
                    if var == 'ji':
                        vec_x = f1[var+plane[0]+sp][int(Nx/2),::n,::n]
                        vec_y = f1[var+plane[1]+sp][int(Nx/2),::n,::n]                    
                    else:
                        vec_x = f1[var+plane[0]][int(Nx/2),::n,::n]
                        vec_y = f1[var+plane[1]][int(Nx/2),::n,::n]                    
                
                # Calculate the magnitude of the electric field vectors
                magnitude = np.sqrt(vec_x**2 + vec_y**2)
                
                # Normalize the vectors to have the same length
                normalized_x = vec_x / magnitude
                normalized_y = vec_y / magnitude
                
                # Define the colormap for coloring the vectors based on magnitude
                #color_map = plt.cm.get_cmap('inferno')
                
                # Plot the vector field in the xy plane with scaled vectors and colored based on magnitude
                plt.quiver(X, Y, normalized_x, normalized_y, ec = 'black'
                           , color='w', scale=20, lw = 0.1, width = 0.006)
                #plt.xlabel('X')
                #plt.ylabel('Y')
                plt.title(f'{var} Field in the {plane} Plane for {molecules[species.index(sp)]}')
                #plt.colorbar(label='Magnitude')
                
                plt.savefig(f'D:/abd__/Downloads D/Bflip/Plots/{var}_{plane}_{sp}/plot-{i}.png')
                plt.show()
                #plt.close()