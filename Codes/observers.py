# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 16:18:54 2023

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

for i in range(np.shape(f1['prob'])[0]):
    os.mkdir(f'D:/abd__/Downloads D/Bflip/Plots/Observers/Observer {i+1}')
    for j in range(1,np.shape(f1['prob'])[2]):
        plt.figure(figsize=(12,4))
        plt.plot([f1['prob'][i,0:-1,0], f2['prob'][i,0:-1,0]]
                 ,[f1['prob'][i,0:-1,j], f2['prob'][i,0:-1,j]]
                 , 'r.', ms = 0.5)
        plt.xlabel('Time')
        plt.ylabel(obsDict[j])
        plt.title(f'Observer {i+1}')
        plt.savefig(f'D:/abd__/Downloads D/Bflip/Plots/Observers/Observer {i+1}/obs{i+1}-{obsDict[j]}_vs_Time.png')
        plt.show()