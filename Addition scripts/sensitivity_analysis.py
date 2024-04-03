# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 11:09:24 2023

@author: js2580
"""
import numpy as np
import matplotlib.pyplot as plt

# from matplotlib import rc
# plt.rcParams["font.family"] = "Times New Roman"

import matplotlib.style
import matplotlib as mpl
mpl.rcParams['font.size'] = 16

# Set the global label size
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['axes.labelsize'] = 14  # Set the y-axis label size



d0 =  3.46
d = np.linspace(3.46, 3.47, 1000)
delta_d = d - d0

strain = ((delta_d)/d0)*100

plt.figure()
plt.plot(delta_d, strain, 'b')
plt.xlabel('$\Delta d_{002} - Å$')
plt.ylabel('Strain - %')
plt.ylim([0, 0.30])
plt.xlim([0, 0.01])
# plt.grid(False)#
plt.show()




d0 =  2.07
d = np.linspace(2.07, 2.08, 1000)
delta_d = d - d0

strain = ((delta_d)/d0)*100

plt.figure()
plt.plot(delta_d, strain, 'b')
plt.xlabel('$\Delta d_{100} - Å$')
plt.ylabel('Strain - %')
plt.ylim([0, 0.50])
plt.xlim([0, 0.01])
# plt.grid(False)


