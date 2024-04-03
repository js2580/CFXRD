# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 16:35:39 2023

@author: js2580
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('test.dat')
intensity = data[:,1]
pixel = data[:,0]

pixelSize = 0.075
SD = 114.4234

twotheta = np.arctan(pixel*pixelSize/SD)
twoTheta = twotheta * 180 / np.pi

plt.figure()
plt.plot(twoTheta,intensity)

np.savetxt('twoThea-Intensity.dat', np.c_[twoTheta,intensity], fmt='%10.5f')

