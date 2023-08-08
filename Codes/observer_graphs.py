# -*- coding: utf-8 -*-
"""
Created on Sat Aug 5 13:54:55 2023
@author: abd__
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import os

f1 = loadmat('D:/abd__/Downloads D/Bflip/Bflip/run0034/observer/Observer_comet_010000_G000.mat')
f2 = np.load('D:/abd__/Downloads D/Bflip/Bflip/run0034/observer/Observer_comet_020000_G000.npz')

obsDict = [
    'time',
    'pBx', 'pBy', 'pBz', 'pBdx', 'pBdy', 'pBdz',
    'pEx', 'pEy', 'pEz',
    'pJx', 'pJy', 'pJz',
    'prho01', 'pjix01', 'pjiy01', 'pjiz01',
    'prho02', 'pjix02', 'pjiy02', 'pjiz02',
    'prho03', 'pjix03', 'pjiy03', 'pjiz03'
]

# Indices of specific values you want to plot (here, the first 3)
specific_values = [1, 2, 3, 7, 8, 9, 10, 13, 17, 21]

i = 43

# Create a figure with the 'sharex' option
fig, axs = plt.subplots(len(specific_values), 1, figsize=(12, 12), sharex=True)

for k, j in enumerate(specific_values):
    ax = axs[k]
    ax.plot(
        [f1['prob'][i-1, 0:-1, 0], f2['prob'][i-1, 0:-1, 0]],
        [f1['prob'][i-1, 0:-1, j], f2['prob'][i-1, 0:-1, j]],
        'r.', ms=0.3
    )
    ax.set_ylabel(obsDict[j])

    # Set the y-axis tick label format to plain without the scale (1e-8)
    ax.yaxis.set_major_formatter('{:.0e}'.format)

# Set x-axis label for the last subplot
axs[-1].set_xlabel('Time')

# Make x-axis tick labels visible only for the last subplot
plt.setp(axs[-1].get_xticklabels(), visible=True)

# Adjust spacing between subplots
plt.tight_layout()

# Show the plot
plt.show()
