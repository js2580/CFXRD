# -*- coding: utf-8 -*-
"""
Created on Sat Aug 20 11:49:46 2022

@author: js2580
"""

# ** Import necessary libraries
# Plot style
import matplotlib.pyplot as plt
from matplotlib import rc
plt.rcParams["font.family"] = "Times New Roman"
import matplotlib as mpl
mpl.rcParams['font.size'] = 14


import numpy as np
import pandas as pd

from CFXRD.CFXRD import *

import seaborn as sns

# ** Define directory paths
ParentDir = '.\\'
Dir = ParentDir + '\\data\\lattice_strain\\\\'


# ** Define scan set numbers
setScan = [389632,389633,389634] # 0% Humpback bridge sample

scanNo = setScan

suffix = '_002_75FWHM_100_120FWHM_fullRing_800bins'

dat = pd.read_csv(Dir +  str(scanNo[0]) +  suffix + '.csv') 
columns = dat.columns

# %%
metaData = CFXRD() 

# Define experimental parameters
pixelSize = 0.075 # Pixel size
SD = 127.032 # Sample to detector distance (mm)
BeamEnergy = 15.2 # Beam energy (keV)
metaData.setDetectorParams(pixelSize = pixelSize, SD = SD, BeamEnergy = BeamEnergy) 

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



###############################################################################

CFXRD.Mapping_Plot(epsilon002, category=Cat, cbarTitle='{002}', cbarLabel='Relative radial {002} lattice strain, (mm/mm)', 
            cbarMax=2E-3, cbarMin=-9e-3, Marker = 'ON', label = 'ON')

###############################################################################
CFXRD.Mapping_Plot(epsilon100, category=Cat, cbarTitle='{100}', cbarLabel='Relative radial {100} lattice strain, (mm/mm)', 
            cbarMax=-3E-3, cbarMin=2e-3, Marker = 'ON', label = 'ON')

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


