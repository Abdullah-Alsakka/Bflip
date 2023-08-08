# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 03:53:30 2023

@author: abd__
"""

import imageio.v2 as imageio

variables = ['B']
coordinates = ['x', 'y', 'z']
species = ['01', '02', '03']
molecules = ['H+', 'He++', 'H2O+']
planes = ['xy', 'xz', 'yz']

frame_duration = 0.1

for plane in planes:
    for sp in species:
        for var in variables:
            with imageio.get_writer(f'D:/abd__/Downloads D/Bflip/Plots/{var}_{plane}_{sp}_CC/GIF.gif', mode='i') as writer:
                for i in range(51):
                    image = imageio.imread(f'D:/abd__/Downloads D/Bflip/Plots/{var}_{plane}_{sp}_CC/plot-{i}.png')
                    
                    for _ in range(int(frame_duration / 0.04)):
                        writer.append_data(image)