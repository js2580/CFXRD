# -*- coding: utf-8 -*-
"""
Created on Sat Aug 20 11:49:46 2022

@author: js2580
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
import matplotlib.markers as markers


import pandas as pd

from CFXRD.CFXRD import *

from matplotlib import rc
plt.rcParams["font.family"] = "Times New Roman"


import matplotlib.style
import matplotlib as mpl
mpl.rcParams['font.size'] = 16

# Set the global label size
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['axes.labelsize'] = 14  # Set the y-axis label size
mpl.rcParams['legend.fontsize'] = 14


import seaborn as sns



# %%
ParentDir = '.\\'
Dir = ParentDir + '\\data\\lattice_strain\\\\'

# 5%
set1 = [389563,389564,389565]
set2 = [389573,389581,389582]
set3 = [389599,389600,389601]
set4 = [389614,389615,389616]

# 0%
set5 = [389632,389633,389634]
set6 = [389639,389640,389641]
set7 = [389646,389647,389648]
set8 = [389653,389654,389655]

setP = [389719,'dump']

setDry = [389730]


scanNo = set7
# set6 
# scanNo = [389633]
if scanNo == set5:
    # start = 32 #central
    # length = 25 #central
    load = '0'
if scanNo == set6:
    # start = 25 #central
    # length = 25 #central
    start = 35
    length = 5
    
    load = '75'
elif scanNo == set7:
    # start = 30 #central
    # length = 25 #central
    start = 38
    length = 5
    
    
    load = '150'


suffix = '_002_75FWHM_100_120FWHM_fullRing_800bins'
dat = pd.read_csv(Dir +  str(scanNo[0]) +  suffix + '.csv') 
columns = dat.columns


metaData = CFXRD() 
metaData.setDetectorParams(pixelSize= 0.075, SD= 127.032, BeamEnergy = 15.2) # in mm and energy in keV 
# SD = 127.168 #mm
# wavelength = 0.8156856 #Angstrom 15.2keV
# pixelSize = 0.075 #mm   #Eiger 4M

gridDat = CFXRD.combineDataintoGridarray(Dir, scanNo, filetype = suffix + '.csv')

Cat = gridDat[:,:,columns.get_loc('Cat')].copy()

#%% 
data = gridDat[:,:,columns.get_loc('xc_CF_002')]
twoTheta, d002 =   metaData.Qspacing_to_Dspacing(data)

if scanNo == set6:
    # 75 N
    DVC_mean = -0.0011998445303358664 * 100 #in percentage
    FEA_mean = -0.004542 * 100 #in percentage
    # d0_range = np.linspace(3.473,3.476,100)
    d0_range = np.linspace(3.4750,3.490,100)

elif scanNo == set7:
    # 150 N
    DVC_mean = -0.002377637686147884 * 100 #in percentage
    FEA_mean = -0.006557 * 100 #in percentage
    # d0_range = np.linspace(3.479,3.481,100)
    d0_range = np.linspace(3.490,3.495,100)

model_mean = FEA_mean
mean = []
for d0 in d0_range:
    epsilon002  = CFXRD.Epsilon_cal(metaData, d002, d0, dtype='d-spacing')
    epsilon002 = metaData.get_neighbor_average(epsilon002, Cat)
    epsilon002 = epsilon002[start:start+length,:]  * 100 #in percentage
    mean.append(np.nanmean(epsilon002))
#%%
plt.figure()
plt.plot(d0_range, mean, color='b', ls='-', lw='2', label='XRD: mean strain')
plt.axhline(model_mean, color='r', ls='--', lw='2', label='DVC: mean strain')
plt.xlabel('$d^0_{002} - Å$')
plt.ylabel('Strain - %')
plt.legend()
plt.xlim([np.min(d0_range), np.max(d0_range)])

d0 = d0_range[np.argmin(np.abs(np.array(mean) - model_mean))] 

print(f'Best d0_002: {d0}')

# d0 = 3.4744242424242424 # same as 75N
# d0 = 3.4936078707599414 #adjust according to Newport distance  12/11/2023
# d0 = 3.466373457440805 #new dry fibre 800bins 29/11/2023
# d0 = 3.4859090909090913 #optimised from FEA 75N
d0 = 3.4945454545454546 #optimised from FEA 150N

epsilon002  = CFXRD.Epsilon_cal(metaData, d002, d0, dtype='d-spacing')
epsilon002 = metaData.get_neighbor_average(epsilon002, Cat)

#Crop Epsilon data to centre part only
epsilon002 = epsilon002[start:start+length,:]

#%%
################################################################################
data = gridDat[:,:,columns.get_loc('xc_CF_100')]
twoTheta, d100 =   metaData.Qspacing_to_Dspacing(data)

if scanNo == set6:
    # 75 N
    DVC_mean = -0.0006340986591155148 * 100 #in percentage
    FEA_mean = -0.000086 * 100 #in percentage
    # d0_range = np.linspace(2.077,2.078,100)
    d0_range = np.linspace(2.076,2.078,100)

elif scanNo == set7:
    # 150 N
    DVC_mean = -0.0017011298574350308 * 100 #in percentage
    FEA_mean = -0.000153 * 100 #in percentage
    # d0_range = np.linspace(2.082,2.083,100)
    d0_range = np.linspace(2.079,2.083,100)


model_mean = FEA_mean
mean = []
for d0 in d0_range:
    epsilon100  = CFXRD.Epsilon_cal(metaData, d100, d0, dtype='d-spacing')
    epsilon100 = metaData.get_neighbor_average(epsilon100, Cat)
    epsilon100 = epsilon100[start:start+length,:] * 100 #in percentage
    mean.append(np.nanmean(epsilon100))
#%%
plt.figure()
plt.plot(d0_range, mean, color='b', ls='-', lw='2', label='XRD: mean strain')
plt.axhline(model_mean, color='r', ls='--', lw='2', label='DVC: mean strain')
plt.xlabel('$d^0_{100} - Å$')
plt.ylabel('Strain - %')
plt.legend()
plt.xlim([np.min(d0_range), np.max(d0_range)])


d0 = d0_range[np.argmin(np.abs(np.array(mean) - model_mean))] 

print(f'Best d0_100: {d0}')

# d0 = 2.077424242424242 # same as 75N
# d0 = 2.0702549548363174 # #adjust according to Newport distance  12/11/2023
# d0 = 2.081591800198599 #new dry fibre 800bins 29/11/2023
# d0 = 2.0762828282828285 # 75N
# d0 = 2.0793636363636367 # 150N

epsilon100  = CFXRD.Epsilon_cal(metaData, d100, d0, dtype='d-spacing')
epsilon100 = metaData.get_neighbor_average(epsilon100, Cat)

#Crop Epsilon data to centre part only
epsilon100 = epsilon100[start:start+length,:]
# %%

fig, ax = plt.subplots()

plt.axis('off')
plt.set_cmap('jet') # plt.set_cmap('rainbow') # plt.set_cmap('tab20c') 

plt.imshow(epsilon002, aspect=('equal')) 
# plt.clim(-4e-3, 4e-3)
plt.colorbar(label = 'Relative radial {002} lattice strain, (mm/mm)', format ='%.2e')
# plt.title('{002}')
plt.tight_layout()

x = epsilon002
plt.figure()
sns.set(style="whitegrid")
hist = sns.histplot(data=x.flatten(), bins=30, facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5, stat="probability")
plt.xlabel("Radial strain")
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


###############################################################################
fig, ax = plt.subplots()

plt.axis('off')
plt.set_cmap('jet') # plt.set_cmap('rainbow') # plt.set_cmap('tab20c') 

plt.imshow(epsilon100, aspect=('equal')) 
# plt.clim(-3e-3, 3e-3)
plt.colorbar(label = 'Relative axial {100} lattice strain, (mm/mm)', format ='%.2e')
# plt.title('{100}')
plt.tight_layout()

x = epsilon100
plt.figure()
sns.set(style="whitegrid")
sns.histplot(data=x.flatten(), bins=30, facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5, stat="probability")
plt.xlabel("Axial strain")
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))



np.savetxt(f'{load}N_lattice_strain.txt', (epsilon002.flatten(), epsilon100.flatten()))
###############################################################################
#%%
err002 = gridDat[:,:,columns.get_loc('err_xc_CF_002')]
err002_avg = np.nanmean(err002)
percentile002 = np.nanpercentile(err002, 90)
print(f'err002_avg: {err002_avg}')
print(f'percentile002: {percentile002}')
# print(f'strain error: {percentile002/3.484}')
print(f'strain error: {2*np.pi/1.8**2*percentile002/3.484}')
print(f'strain range: {np.nanmax(epsilon002) - np.nanmin(epsilon002)}')
print(f'strain std: {np.nanstd(epsilon002)}')
print(f'002 mean: {np.nanmean(epsilon002)}')
print('####################')

err100 = gridDat[:,:,columns.get_loc('err_xc_CF_100')]
err100_avg = np.nanmean(err100)
percentile100 = np.nanpercentile(err100, 90)
print(f'err100_avg: {err100_avg}')
print(f'percentile100: {percentile100}')
# print(f'strain error: {percentile100/2.0758}')
print(f'strain error: {2*np.pi/3.02**2*percentile100/2.0758}')
print(f'strain range: {np.nanmax(epsilon100) - np.nanmin(epsilon100)}')
print(f'strain std: {np.nanstd(epsilon100)}')
print(f'100 mean: {np.nanmean(epsilon100)}')

plt.show()


print(f'{np.nanstd(epsilon002)}')