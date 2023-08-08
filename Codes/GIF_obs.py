# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 14:07:06 2023

@author: abd__
"""

import numpy as np
import imageio.v2 as imageio

obsDict = [
                'time',
                'pBx' ,
                'pBy',
                'pBz',
                'pBdx',
                'pBdy',
                'pBdz',
                'pEx',
                'pEy',
                'pEz',
                'pJx',
                'pJy',
                'pJz',
                'prho01', 
                'pjix01',
                'pjiy01', 
                'pjiz01',
                'prho02',
                'pjix02',
                'pjiy02',
                'pjiz02',
                'prho03',
                'pjix03',
                'pjiy03',
                'pjiz03'
                    ]

frame_duration = 0.2

for var in obsDict[1:-1]:
    with imageio.get_writer(f'D:/abd__/Downloads D/Bflip/Plots/Observers/gif_{var}.gif', mode='i') as writer:

        order = np.arange(60,30,-1).tolist() + np.arange(1,31).tolist()        
        for obs in order:
                
                image = imageio.imread(f'D:/abd__/Downloads D/Bflip/Plots/Observers/Observer {obs}/obs{obs}-{var}_vs_Time.png')
                    
                for _ in range(int(frame_duration / 0.04)):
                    writer.append_data(image)