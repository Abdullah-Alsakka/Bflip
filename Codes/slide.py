# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 17:15:24 2023

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
from matplotlib.widgets import Slider, Button
from scipy.ndimage import gaussian_filter
import h5py
import os

path = "D:/abd__/Downloads D/Bflip/Bflip/run0034/out"
file_list = sorted(os.listdir(path))  # Sort the file list

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

x_axis = locals()[plane[0]]
y_axis = locals()[plane[1]]

# Create the figure and the line that we will manipulate
fig, ax = plt.subplots()

# adjust the main plot to make room for the sliders
fig.subplots_adjust(left=0.25, bottom=0.25)

# Initialize the colorbar
colorbar = None

# Make a horizontal slider to control the time.
axtime = fig.add_axes([0.25, 0.1, 0.65, 0.03])
time_slider = Slider(
    ax=axtime,
    label='Time',
    valmin=0,
    valmax=len(file_list)-1,
    valinit=0,
    valstep=1,
    closedmax=True,
)

# The function to be called anytime the time slider's value changes
def update(val):
    index = int(time_slider.val)  # Get the index of the selected file
    filename = os.path.join(path, file_list[index])
    f1 = h5py.File(filename, 'r+')
    
    if plane == 'xy':
        rho = f1['rho'+species][:,:,int(Nz/2)]    
    elif plane == 'xz':
        rho = f1['rho'+species][:,int(Ny/2),:]
    elif plane == 'yz':
        rho = f1['rho'+species][int(Nx/2),:,:]
        
    rho = gaussian_filter(rho, 1)
    
    ax.clear()
    
    X, Y = np.meshgrid(x_axis, y_axis, indexing='ij')
    ax.pcolormesh(X, Y, rho, shading='auto', cmap='jet')
    
    ''' 
    # Define the subsampling factor
    n = 10
    
    # Generate the coordinates for all the grid points
    X, Y = np.meshgrid(x_axis[::n]
                       , y_axis[::n]
                       , indexing='ij')
    
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
    ax.quiver(X, Y, normalized_x, normalized_y, ec='black'
              , color='w', scale=20, lw=0.1, width=0.006)
    '''
    '''     
    # Observers
    
    f2 = np.load('D:/abd__/Downloads D/Bflip/Bflip/run0034/observer/Observer_comet_020000_G000.npz')
    
    if plane == 'xz':
        ax.plot(f2['observers'][:,0], f2['observers'][:,2], 'r.', ms = 1)
    elif plane == 'yz':
        ax.plot(f2['observers'][:,1], f2['observers'][:,2], 'r.', ms = 1)
        
    #ax.set_xlim(-3e6, 3e6)
    ax.set_ylim(min(f2['observers'][:,2]), max(f2['observers'][:,2]))
    '''
    fig.canvas.draw_idle()

# Register the update function with the time slider
time_slider.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')

def reset(event):
    time_slider.reset()
button.on_clicked(reset)

plt.show()
