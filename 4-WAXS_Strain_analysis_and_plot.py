# -*- coding: utf-8 -*-
"""
Created on Sat Aug 20 11:49:46 2022

@author: js2580
"""

# ** Import necessary libraries
# Plot style
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
import matplotlib.markers as markers
from matplotlib.ticker import ScalarFormatter
from matplotlib import rc
plt.rcParams["font.family"] = "Times New Roman"
import matplotlib.style
import matplotlib as mpl
mpl.rcParams['font.size'] = 14

class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
             self.format = r'$\mathdefault{%s}$' % self.format

import numpy as np
import pandas as pd

from CFXRD.CFXRD import *

import seaborn as sns

# ** Define directory paths
ParentDir = 'C:\\Users\\js2580\\OneDrive - University of Bath\\Phd\\Research\\Synchrotron\\XRD\\MM30528 - doping\\WAXS analysis\\Analysis'
Dir = ParentDir + '\\data\\lattice_strain\\\\'


# ** Define scan set numbers
set5 = [389632,389633,389634] # 0% Humpback bridge sample

scanNo = set5

suffix = '_002_75FWHM_100_120FWHM_fullRing_800bins'

dat = pd.read_csv(Dir +  str(scanNo[0]) +  suffix + '.csv') 
columns = dat.columns

# %%
metaData = CFXRD() 
metaData.setDetectorParams(pixelSize= 0.075, SD= 127.032, BeamEnergy = 15.2) # in mm and energy in keV 
# SD = 127.168 #mm
# wavelength = 0.8156856 #Angstrom 15.2keV
# pixelSize = 0.075 #mm   #Eiger 4M

# ** Combine all set into a single grid array arranging in term of motors positions
gridDat = CFXRD.combineDataintoGridarray(Dir, scanNo, filetype = suffix + '.csv')

Cat = gridDat[:,:,columns.get_loc('Cat')].copy()


#### 002 Peak ####
data = gridDat[:,:,columns.get_loc('xc_CF_002')]
twoTheta, d002 =   metaData.Qspacing_to_Dspacing(data)
d0 = 3.4945454545454546 #optimised from FEA 150N
epsilon002  = CFXRD.Epsilon_cal(metaData, d002, d0, dtype='d-spacing')
# epsilon002 = metaData.get_neighbor_average(epsilon002, Cat)


#### 100 Peak ####
data = gridDat[:,:,columns.get_loc('xc_CF_100')]
twoTheta, d100 =   metaData.Qspacing_to_Dspacing(data)
d0 = 2.0793636363636367 #optimised from FEA 150N
epsilon100  = CFXRD.Epsilon_cal(metaData, d100, d0, dtype='d-spacing')
# epsilon100 = metaData.get_neighbor_average(epsilon100, Cat)



# %%

fig, ax = plt.subplots()
# Create a custom marker
marker = markers.MarkerStyle('x')

# Overlay the marker on the image
# for i in range(Cat.shape[0]):
#     for j in range(Cat.shape[1]):
#         if (Cat[i, j] == 'Fibre') and (np.isnan(epsilon002[i,j])):
#             ax.plot(j, i, marker=marker, color='k', markersize=6, linewidth=12)

# for i in range(Cat.shape[0]):
#     for j in range(Cat.shape[1]):
#         if gridDat[:,:,columns.get_loc('err_xc_CF_002')][i,j] > 0.000555:
#             ax.plot(j, i, marker=marker, color='k', markersize=6, linewidth=12)



plt.axis('off')
plt.set_cmap('jet') # plt.set_cmap('rainbow') # plt.set_cmap('tab20c') 

plot = plt.imshow(epsilon002, aspect=('equal')) 
plt.clim(-9e-3, 0)

cbar = fig.colorbar(plot, format=OOMFormatter(-3, mathText=False))
plt.colorbar(label = 'Relative radial {002} lattice strain, (mm/mm)', format ='%.2e')
plt.title('{002}')
plt.tight_layout()


###############################################################################
fig, ax = plt.subplots()
# Create a custom marker
marker = markers.MarkerStyle('x')

# Overlay the marker on the image
for i in range(Cat.shape[0]):
    for j in range(Cat.shape[1]):
        if (Cat[i, j] == 'Fibre') and (np.isnan(epsilon100[i,j])):
            ax.plot(j, i, marker=marker, color='k', markersize=6, linewidth=12)

# for i in range(Cat.shape[0]):
#     for j in range(Cat.shape[1]):
#         if gridDat[:,:,columns.get_loc('err_xc_CF_100')][i,j] > 0.000469:
#             ax.plot(j, i, marker=marker, color='k', markersize=6, linewidth=12)

# formatter = ScalarFormatter(useMathText=True)
# formatter.set_scientific(True)
# formatter.set_powerlimits((-3, 3))  # Adjust the power limits as needed
# plt.gca().yaxis.set_major_formatter(formatter)

plt.axis('off')
plt.set_cmap('jet') # plt.set_cmap('rainbow') # plt.set_cmap('tab20c') 

plot = plt.imshow(epsilon100, aspect=('equal')) 
plt.clim(-3e-3, 2e-3)

cbar = fig.colorbar(plot, format=OOMFormatter(-3, mathText=False))
plt.colorbar(label = 'Relative axial {100} lattice strain, (mm/mm)', format ='%.2e')
plt.title('{100}')
plt.tight_layout()

###############################################################################
#%%
q002 = gridDat[:,:,columns.get_loc('xc_CF_002')]
err002 = gridDat[:,:,columns.get_loc('err_xc_CF_002')]
err002_avg = np.nanmean(err002)
percentile002 = np.nanpercentile(err002, 90)
print(f'err002_avg: {err002_avg}')
print(f'percentile002: {percentile002}')
# print(f'strain error: {percentile002/3.484}')
# print(f'strain error: {2*np.pi/1.8**2*percentile002/3.484}')
print(f'strain error: {np.nanpercentile(1.8*err002/q002**2, 90) * 100}%')
print(f'strain range: {np.nanmax(epsilon002) - np.nanmin(epsilon002)}')

print('####################')
q100 = gridDat[:,:,columns.get_loc('xc_CF_100')]
err100 = gridDat[:,:,columns.get_loc('err_xc_CF_100')]
err100_avg = np.nanmean(err100)
percentile100 = np.nanpercentile(err100, 90)
print(f'err100_avg: {err100_avg}')
print(f'percentile100: {percentile100}')
# print(f'strain error: {percentile100/2.0758}')
# print(f'strain error: {2*np.pi/3.02**2*percentile100/2.0758}')
print(f'strain error: {np.nanpercentile(3.0*err100/q100**2, 90) * 100}%')
# print(f'strain distribution: {np.nanstd(epsilon100)}')
print(f'strain range: {np.nanmax(epsilon100) - np.nanmin(epsilon100)}')


#%%
print('####################')
print(f'Total points: {np.count_nonzero(np.isnan(err002.astype(float)))}')


plt.show()


#%%
# print('####################')
# err_angle = gridDat[:,:,columns.get_loc('err_angle2')]
# errangle_avg = np.nanmean(err_angle)
# percentileangle = np.nanpercentile(err_angle, 90)
# print(f'errAngle_avg: {errangle_avg}')
# print(f'percentileAngle: {percentileangle}')
# # print(f'strain error: {percentile002/3.484}')
# # print(f'strain error: {2*np.pi/1.8**2*percentile002/3.484}')
