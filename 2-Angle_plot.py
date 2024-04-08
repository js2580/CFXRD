# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 09:26:28 2023

@author: js2580
"""

# ** Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.signal import find_peaks
import h5py 
from tqdm import tqdm

from CFXRD.CFXRD import *

# ** Define parent directory path
ParentDir = '.\\'
# ** Define input directory path
InputDir = ParentDir + '\\data\\fibre_angle\\'

# ** Define column names
ColumnNames = ['ScanNo', 'PointNo', 'Xmotor', 'Ymotor', 'Cat',\
           'angle1', 'err_angle1','angle_redchi1', 'FWHM_angle1',\
           'angle2', 'err_angle2', 'angle_redchi2', 'FWHM_angle2',\
           'xc_CF_002', 'err_xc_CF_002', 'chi_CF_002',\
           'xc_CF_100', 'err_xc_CF_100', 'chi_CF_100', \
           'xc_CF_100*', 'err_xc_CF_100*', 'chi_CF_100*',\
           'xc_CF_110', 'err_xc_CF_110', 'chi_CF_110',\
           'xc_CF_004', 'err_xc_CF_004', 'chi_CF_004',\
           ]
    

# ** Define scan set numbers
setScan = [389632,389633,389634]

full_set = setScan

fig = plt.figure()

plt.axis('equal')
plt.axis('off')
fig.set_size_inches(10, 10)

for scanNo in full_set:
     # Call functions
     metaData = FibrePlot()
     df = pd.read_csv(InputDir +  str(scanNo) +  '.csv')
     
     metaData.motorPosition(df['Xmotor'], df['Ymotor'], df['Cat'], df['angle1'], df['angle2'], df['angle_redchi1'], df['angle_redchi2'])
     metaData.Plot() # Plot fibre orientation
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    