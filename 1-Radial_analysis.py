# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 10:11:32 2023

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

# ** Define output directory path
outDir = ParentDir + '\\data\\fibre_angle\\'
# ** Define cake NEXUS file input directory path
nexus_dir = ParentDir + '\\data\\caking\\_eigerScan'
# ** Define motor input directory path
motor_dir = ParentDir + '\\data\\motor\\'
# ** Define motor names
motorX = 'gtbX3'
motorY = 'gtbY2'

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

"""
# Multiple scanning points or Single scanning point
# singleScan: True --> ignore motor positions (NAN value)
# singleScan: True --> read motor positions (True value)
"""
singleScan = False

# ** Define scan set numbers
setScan = [389632,389633,389634] # 0% Humpback bridge sample

# ** Define multiple scan numbers
scanlist = setScan 

# ! Loop through each scanning set
for scanNo in scanlist:
    # Call functions
    metaData = CFXRD() 
    # Define experimental parameters
    metaData.setDetectorParams(pixelSize = 0.075, SD= 127.032, BeamEnergy = 15.2) # in mm and energy in keV
    
    # Open and Read data from NEXUS file 
    dataset = h5py.File( nexus_dir + str(scanNo)  + '_caking.nxs' , 'r')
    Intensity = np.asarray(dataset['processed/result/data'])
    Angle = np.asarray(dataset['processed/result/azimuthal angle (degrees)'])
    Q = np.asarray(dataset['processed/result/q'])
    
    # Close NEXUS file
    dataset.close()
    
    metaData.readCakeData(Intensity = Intensity, AzimuthalAngle = Angle, RadialDist = Q)
    # Define q range for 002
    minQ = 1.5
    maxQ = 2.0
    # Define q range for 100
    # minQ = 2.8
    # maxQ = 3.25
    metaData.cake_RadialIntegration(radialMax = maxQ, radialMin = minQ, angleMax = 450, angleMin = 90)
    
    metaData.GenResultArray(ColumnNames = ColumnNames)
    
    # ! Loop through scanning points in each scanning set (z-axis)
    for pointNo in tqdm(range(metaData.TotalSlice)):
    # for pointNo in [1543]:
    # for pointNo in range(810,820):

        x = metaData.AzimuthalAngle
        y = metaData.Intensity_Radial[pointNo]
    
        ############# @@@@@@@@@@@@@ #############
        # Plot original data
        # plt.figure()
        # plt.plot(x, y, '*-b')
        ############# @@@@@@@@@@@@@ #############
        
        # Differentiate intensity to get sharp slope drop 
        ydiff = np.diff(y)
        ydiffdiff = np.diff(ydiff)
        
        # ** Filter noise or dead pixels 
        # Todo: these values need to be adjusted to specific data
        idx = np.ones(metaData.AzimuthalBin-1, dtype=int)
        idx[ydiff >= 200 ] = 0
        idx[ydiff <= -200] = 0
        idx[ydiff == 0] = 0
        idx = np.append(idx, [0])
        idx[y == 0 ] = 0
        
        # for filter in range(0,1):
        #     y = y[idx == 1]
        #     x = x[idx == 1]
            
        #     idx = np.ones(len(y)-1, dtype=int)
        #     ydiff = np.diff(y)
        #     idx[ydiff >= 1] = 0
        #     idx[ydiff <= -1] = 0
        #     idx[ydiff == 0] = 0
        #     idx = np.append(idx, [0])
            
        #####################################
        # Filter outliner statistically#
            
        # Q1, Q3 = np.quantile(y[idx==1], [0.25, 0.75])
        # IQR = Q3 - Q1
        # lower = Q1 -  1.5 * IQR
        # upper = Q3 +  1.5 * IQR
        
        mean = np.mean(y[idx==1])
        # sigma = np.std(y[idx==1]) 
        # lower = mean - 4 * sigma
        # upper = mean + 4 * sigma
        
        # idx[y > upper ] = 0
        # idx[y < lower ] = 0
        #####################################
        
        # ** Automatic find peak index
        # peak_pos = find_peaks(y[idx==1], distance=30, height=(80, ), prominence=(10, ), width=(5, ))[0]
        peak_pos = find_peaks(y[idx==1], distance=20, height=(80, ), prominence=(10, ), width=(5, ))[0]


        ############# @@@@@@@@@@@@@ #############
        ##Plot filter peaks
        plt.figure()
        plt.plot(x[idx==1], y[idx==1], 'g.')    
        plt.xlabel('Azimuthal angle (degrees)')
        plt.ylabel('Intensity counts')
        plt.plot(x[idx==1][peak_pos], y[idx==1][peak_pos], 'rx')
       ############# @@@@@@@@@@@@@ #############
       
        # if len(peak_pos) == 2 and mean >= 40 and np.max(y[idx==1][peak_pos])>=80:
        if len(peak_pos) == 2 and np.max(y[idx==1][peak_pos])>=80:
            cat = 'Fibre'            
            #### Peak No 1 ####
            # ** Define fitting range
            peak = x[idx==1][peak_pos[0]] # degrees
            offset = 40 # degrees
            # LowerLim = np.absolute(x[idx==1] - (peak - offset)).argmin()
            # UpperLim = np.absolute(x[idx==1] - (peak + offset)).argmin()
            LowerLim = (peak - offset)
            UpperLim = (peak + offset)

            # ** Define fitting functions
            spec = {
                'x': x[idx==1],
                'y': y[idx==1],
                # 'x': x,
                # 'y': y,
                'model': [
                            {'type': 'PseudoVoigtModel',
                                   'params': {
                                        # 'center'    : peak,
                                        # 'sigma'     : (x[idx==1][UpperLim] - x[idx==1][LowerLim])/2,
                                        # 'height'    : y[idx==1][LowerLim:UpperLim].max(),
                                        'fraction'  : 1
                                       },
                                   'help':{
                                       # 'center'  : {'min': 190,'max': 220},
                                       # 'sigma'   : {'min': 0},
                                       # 'amplitude': {'max' : }
                                       # 'fractioin' : {'min' : 1}
                                       }
                            },
                        ]
    
                }

            FittingOutput = metaData.PeakModelGen(spec = spec, LowerLim = LowerLim, UpperLim = UpperLim)
            # FittingOutput.fit_report() # !uncomment to show fitting report 
            # metaData.PeakResidualPlot(FittingOutput) # ! uncomment to show residual plot
            
            # ** Retrieve and assign values to be stored 
            xc_angle1 = FittingOutput.params['m0_center'].value
            xcstd_angle1 = FittingOutput.params['m0_center'].stderr
            xcchi_angle1 = FittingOutput.redchi
            fwhm_angle1  = FittingOutput.params['m0_fwhm'].value
            
            
            metaData.ResultArray.loc[pointNo, 'angle1'] = xc_angle1
            metaData.ResultArray.loc[pointNo, 'err_angle1'] = xcstd_angle1
            metaData.ResultArray.loc[pointNo, 'angle_redchi1'] = xcchi_angle1
            metaData.ResultArray.loc[pointNo, 'FWHM_angle1'] = fwhm_angle1
            
            ###### Peak No 2 ######
            # ** Define fitting range
            peak = x[idx==1][peak_pos[1]] # degrees
            offset = 40 # degrees
            LowerLim = (peak - offset)
            UpperLim = (peak + offset)

            # ** Define fitting functions
            spec = {
                'x': x[idx==1],
                'y': y[idx==1],
                # 'x': x,
                # 'y': y,
                'model': [
                            {'type': 'PseudoVoigtModel',
                                    'params': {
                                        # 'center'    : peak,
                                        # 'sigma'     : (x[idx==1][UpperLim] - x[idx==1][LowerLim])/2,
                                        # 'height'    : y[idx==1][LowerLim:UpperLim].max(),
                                        'fraction'  : 1
                                        },
                                    'help':{
                                        # 'center'  : {'min': 190,'max': 220},
                                        # 'sigma'   : {'min': 0},
                                        # 'amplitude': {'max' : }
                                        # 'fractioin' : {'min' : 1}
                                        }
                            },
                        ]
        
                }
        
            FittingOutput = metaData.PeakModelGen(spec = spec, LowerLim = LowerLim, UpperLim = UpperLim)
            # FittingOutput.fit_report() # !uncomment to show fitting report 
            # metaData.PeakResidualPlot(FittingOutput) # ! uncomment to show residual plot
            
            # ** Retrieve and assign values to be stored 
            xc_angle2 = FittingOutput.params['m0_center'].value
            xcstd_angle2 = FittingOutput.params['m0_center'].stderr
            xcchi_angle2 = FittingOutput.redchi
            fwhm_angle2  = FittingOutput.params['m0_fwhm'].value
            
            
            metaData.ResultArray.loc[pointNo, 'angle2'] = xc_angle2
            metaData.ResultArray.loc[pointNo, 'err_angle2'] = xcstd_angle2
            metaData.ResultArray.loc[pointNo, 'angle_redchi2'] = xcchi_angle2
            metaData.ResultArray.loc[pointNo, 'FWHM_angle2'] = fwhm_angle2
                
        elif mean > 15:
            cat = 'Resin'
        else:
            cat = 'Off'
        if singleScan == True:
            pass
            x_pos = np.nan
            y_pos = np.nan
        else:
            x_pos, y_pos = metaData.motor_position(motor_dir, motorX, motorY, scanNo, pointNo)

        
        metaData.ResultArray.loc[pointNo, 'ScanNo'] = scanNo
        metaData.ResultArray.loc[pointNo, 'PointNo'] = pointNo
        metaData.ResultArray.loc[pointNo, 'Xmotor'] = x_pos
        metaData.ResultArray.loc[pointNo, 'Ymotor'] = y_pos
        metaData.ResultArray.loc[pointNo, 'Cat'] = cat
        
    metaData.ResultArray.to_csv(outDir + str(scanNo) + '.csv', na_rep='NaN', index=False)
    
        
        
        
        
        
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    