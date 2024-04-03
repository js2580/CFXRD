# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 11:12:59 2023

@author: js2580
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import lmfit

#import NEXSUS/h5 file reader
import h5py 

ParentDir = 'C:\\Users\\js2580\\OneDrive - University of Bath\\Phd\Research\\Synchrotron\\XRD\\WAXS analysis\\Analysis\\'
nexus_dir = ParentDir + '\\data\\caking\dry fibre\\_eigerScan'
scanNo = 389730

pixelSize= 0.075
SD= 127.032 
BeamEnergy = 15.2
wavelength = 6.62607015E-34 * 3E8 / (BeamEnergy * 1.60218E-16 ) * 10**10

# Open and Read data from NEXUS file 
dataset = h5py.File( nexus_dir + str(scanNo)  + '_caking.nxs' , 'r')
Intensity = np.asarray(dataset['processed/result/data'])
# dat2 = np.asarray(dataset['processed/intermediate/4-Azimuthal Integration/data']) # Full integration
Q = np.asarray(dataset['processed/result/q'])
Angle = np.asarray(dataset['processed/result/azimuthal angle (degrees)'])

# Close NEXUS file
dataset.close()

d = 2*np.pi/Q
twotheta = 2*np.arcsin(wavelength/(2*d)) * 180 / np.pi


array_with_nan = np.where(Intensity[0] == 0, np.nan, Intensity[0])
average_with_nan = np.nanmean(array_with_nan, axis=0)
average = np.where(np.isnan(average_with_nan), 0, average_with_nan)

plt.figure()
plt.plot(Q, average)