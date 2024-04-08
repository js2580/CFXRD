# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 13:03:29 2023

@author: js2580
"""


# ** Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import lmfit
import h5py 
from tqdm import tqdm

import time

from CFXRD.CFXRD import *

# ** Define directory paths
ParentDir = '.\\'

outDir = ParentDir + '\\data\\lattice_strain\\'
output_suffix = '_eta.csv'
nexus_dir = ParentDir + '\\data\\caking\\_eigerScan'
angle_dir = ParentDir + '\\data\\fibre_angle\\'


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
setScan = [389632,389633,389634] # 0% Humpback bridge sample

scanlist = setScan

fit002 = True
fit100 = True
fit100s = False
fit110 = False

# ! Loop through each scanning set
for scanNo in scanlist:
    metaData = CFXRD() 
    metaData.setDetectorParams(pixelSize= 0.075, SD= 127.032, BeamEnergy = 15.2) # in mm and energy in keV
    
    # Open and Read data from NEXUS file 
    dataset = h5py.File( nexus_dir + str(scanNo)  + '_caking.nxs' , 'r')
    Intensity = np.asarray(dataset['processed/result/data'])
    # dat2 = np.asarray(dataset['processed/intermediate/4-Azimuthal Integration/data']) # Full integration
    Q = np.asarray(dataset['processed/result/q'])
    Angle = np.asarray(dataset['processed/result/azimuthal angle (degrees)'])
    
    # Close NEXUS file
    dataset.close()
    
    # ** Call function to correct for two theta shift
    # correct q-spacing from sample displacement
    twoTheta, d = metaData.Qspacing_to_Dspacing(Q)
    # new 2theta
    eta = metaData.twoThetaShift(disp=1, R=127.032, twoTheta=twoTheta)
    Q = metaData.twoTheta_to_Qspacing(twoTheta-eta)
    
    metaData.readCakeData(Intensity = Intensity, AzimuthalAngle = Angle, RadialDist = Q)
    
    minQ = Q.min()
    maxQ = Q.max()
    metaData.cake_AzimuthalIntegration(radialMax = maxQ, radialMin = minQ, angleMax = 450, angleMin = 90)
    
    metaData.GenResultArray(ColumnNames = ColumnNames)
    
    # ! Check if there is an exisiting file
    if os.path.isfile(outDir + str(scanNo) + output_suffix) == True:
        #Assigned exisisting date to Result array
        metaData.ResultArray = pd.read_csv(outDir + str(scanNo) + output_suffix)
    else:
        # if there is no existing file
        # read fibre orientatin angle from angle output file
        df = pd.read_csv(angle_dir +  str(scanNo) +  '.csv')
        #Assigned to the first few columns
        metaData.ResultArray.iloc[:,:len(df.columns)] = df 
        
        
    # Read from excel file
    angle1 = metaData.ResultArray['angle1']
    angle2 = metaData.ResultArray['angle2']

    
    angle_chi1 = metaData.ResultArray['angle_redchi1']
    angle_chi2 = metaData.ResultArray['angle_redchi2']
    
    fwhm_angle1 = metaData.ResultArray['FWHM_angle1']    
    fwhm_angle2 = metaData.ResultArray['FWHM_angle2']    

    
    cat = metaData.ResultArray['Cat']
    
    t1 = time.perf_counter()
    
    # ! Loop through scanning points in each scanning set (z-axis)
    for pointNo in tqdm(range(metaData.TotalSlice)):
    # for pointNo in [1312]:
    # for pointNo in range(700,705):
        print(f'im No: {pointNo}')
        if cat[pointNo] == 'Fibre':
            if np.absolute(1 - angle_chi1[pointNo]) < np.absolute(1 - angle_chi2[pointNo]):
                angle = angle1[pointNo] 
                fwhm = fwhm_angle1[pointNo]
            else:
                angle = angle2[pointNo] 
                fwhm = fwhm_angle2[pointNo]
                
            index_A, index_invA, index_B, index_invB = metaData.Anisotropic_selective_angle(PeakAngle = angle, fwhmA = 75, fwhmB = 120) # Pristine

            allA = metaData.Integrate_selective_angle_intensity(metaData.Intensity_Azimuthal[pointNo][index_A + index_invA])
            allB = metaData.Integrate_selective_angle_intensity(metaData.Intensity_Azimuthal[pointNo][index_B + index_invB])

            
            # plt.figure()
            # for l in [79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]:
            #   plt.plot(Q,metaData.Intensity[pointNo][l])

            ############# @@@@@@@@@@@@@ #############
            
            # plt.figure()
            # plt.plot(Q, allA, label= 'All A')
            # plt.plot(Pixel, metaData.Integrate_selective_angle_intensity(metaData.Intensity[pointNo][index_A + index_invA]), label= 'All A')
                       
            # plt.plot(Q, allB, label= 'All B')
            # plt.plot(Q, metaData.Integrate_selective_angle_intensity(metaData.Intensity[pointNo][index_B + index_invB]), label= 'All B')
            # allB = metaData.Integrate_selective_angle_intensity(metaData.Intensity[pointNo][index_B + index_invB])
            # plt.legend()
            
            # plt.figure()
            # plt.plot(Pixel, subAll, '-x', label = 'abs subtract A - B')
            # # plt.plot(Pixel, allA, '-x', label = 'all A')
            # plt.legend()
            
            # ######################################
            '''
            Perfrom curve fitting: #############  CF 002 peak  #############
            '''
            if fit002 == True:
                x = Q
                y = allA
                
                ############# @@@@@@@@@@@@@ #############
                # plt.figure()
                # plt.plot(x, y)
                #######################################
                p002 = 1.8 #d002 = 3.48
                LowerLim = 1.2
                UpperLim = 2.2
                # LowerLim = 270
                # UpperLim = 470
                
                spec = {
                        'x': x,
                        'y': y,
                        'model': [
                                    {'type': 'PseudoVoigtModel',
                                            'params': {
                                                'center'  : p002,
                                            #      # 'sigma'   : 50
                                                },
                                            'help':{
                                                'center'  : {'min': 1.6},
                                                # 'sigma'   : {'min': 0.5},
                                                }
                                    },
                                    {'type': 'PseudoVoigtModel',
                                            'params': {
                                                'center'  : 1.3,
                                                    'sigma'   : 1,
                                                },
                                                'help':{
                                                    'center'  : {'max': 1.6},
                                                # 'sigma'   : {'min': 30},
                                                }
                                    },
                                    
                            ]
                        }
                    
                FittingOutput = metaData.PeakModelGen(spec = spec, LowerLim = LowerLim, UpperLim = UpperLim )
                # FittingOutput.fit_report() # ! uncomment to show fitting report 
                # metaData.PeakResidualPlot(FittingOutput) # ! uncomment to show residual plot
                
                # ** Retrieve and assign values to be stored 
                # Peak centre
                metaData.ResultArray.loc[pointNo, 'xc_CF_002'] = FittingOutput.params['m0_center'].value
                metaData.ResultArray.loc[pointNo, 'err_xc_CF_002'] = FittingOutput.params['m0_center'].stderr
                metaData.ResultArray.loc[pointNo, 'chi_CF_002'] = FittingOutput.redchi # FittingOutput.chisqr
                # FWHM 
                if spec["model"][0]["type"] == 'GaussianModel':    
                    FWHM = 2 * FittingOutput.params['m0_sigma'].value * np.sqrt(2 * np.log(2)) # Gaussian
                elif spec["model"][0]["type"] == 'PseudoVoigtModel':
                    FWHM = FittingOutput.params['m0_sigma'].value * np.sqrt(2 * np.log(2)) # PseudoVoigtModel
                metaData.ResultArray.loc[pointNo, 'FWHM_002'] = FWHM
            # ######################################
            
            '''
            Perfrom curve fitting: #############  CF 100 peak  #############
            '''
            if fit100 == True:
                x = Q
                y = allB
                
                ############# @@@@@@@@@@@@@ #############
                # plt.figure()
                # plt.plot(x, y)
                #######################################
                
                p100 = 3.0
                # # LowerLim = 2.4
                # # UpperLim = 3.5
                LowerLim = 2.7
                UpperLim = 3.3
                
                spec = {
                        'x': x,
                        'y': y,
                        'model': [
                                    {'type': 'PseudoVoigtModel',
                                            'params': {
                                                'center'  : p100,
                                                'sigma'   : 0.1
                                                },
                                                'help':{
                                                'center'  : {'max': 3.2},
                                                'sigma'   : {'max': 0.4},
                                                }
                                    },
                                    {'type': 'GaussianModel',
                                            'params': {
                                                'center'  : 2.8,
                                                'sigma'   : 0.55,
                                                },
                                            'help':{
                                                # 'center'  : {'min': 3.2},
                                                'sigma'   : {'min': 0.4},
                                                }
                                    },
                                    
                            ]
                        }

                FittingOutput = metaData.PeakModelGen(spec = spec, LowerLim = LowerLim, UpperLim = UpperLim )
                # FittingOutput.fit_report() # ! uncomment to show fitting report 
                # metaData.PeakResidualPlot(FittingOutput) # ! uncomment to show residual plot
                
                # ** Retrieve and assign values to be stored 
                # Peak Centre
                metaData.ResultArray.loc[pointNo, 'xc_CF_100'] = FittingOutput.params['m0_center'].value
                metaData.ResultArray.loc[pointNo, 'err_xc_CF_100'] = FittingOutput.params['m0_center'].stderr
                metaData.ResultArray.loc[pointNo, 'chi_CF_100'] = FittingOutput.redchi # FittingOutput.chisqr
                # FWHM 
                if spec["model"][0]["type"] == 'GaussianModel':    
                    FWHM = 2 * FittingOutput.params['m0_sigma'].value * np.sqrt(2 * np.log(2)) # Gaussian
                elif spec["model"][0]["type"] == 'PseudoVoigtModel':
                    FWHM = FittingOutput.params['m0_sigma'].value * np.sqrt(2 * np.log(2)) # PseudoVoigtModel
                metaData.ResultArray.loc[pointNo, 'FWHM_100'] = FWHM

                # ######################################
                
            '''
            Perfrom curve fitting: #############  CF 100* peak  #############
            '''
            if fit100s == True:
                
                x = Q
                y = allA
                
                ############# @@@@@@@@@@@@@ #############
                # plt.figure()
                # plt.plot(x, y)
                #######################################
                
                p100s = 3.0
                LowerLim = 2.7
                UpperLim = 3.4
                
                spec = {
                        'x': x,
                        'y': y,
                        'model': [
                                    {'type': 'PseudoVoigtModel',
                                            'params': {
                                                'center'  : p100s,
                                                # 'sigma'   : 30
                                                },
                                            'help':{
                                                # 'center'  : {'min': p100s - 10, 'max': p100s + 10},
                                                # 'sigma'   : {'min': 10, 'max': 80},
                                                }
                                    },
                                    {'type': 'PseudoVoigtModel',
                                            'params': {
                                                'center'  : 2.8,
                                                    'sigma'   : 0.55,
                                                },
                                            'help':{
                                                # 'center'  : {'max': 370},
                                                'sigma'   : {'min': 0.5},
                                                }
                                    },
                                    
                            ]
                        }
                
                FittingOutput = metaData.PeakModelGen(spec = spec, LowerLim = LowerLim, UpperLim = UpperLim )
                # FittingOutput.fit_report() # ! uncomment to show fitting report 
                # metaData.PeakResidualPlot(FittingOutput) # ! uncomment to show residual plot
                
                # ** Retrieve and assign values to be stored 
                metaData.ResultArray.loc[pointNo, 'xc_CF_100*'] = FittingOutput.params['m0_center'].value
                metaData.ResultArray.loc[pointNo, 'err_xc_CF_100*'] = FittingOutput.params['m0_center'].stderr
                metaData.ResultArray.loc[pointNo, 'chi_CF_100*'] = FittingOutput.redchi # FittingOutput.chisqr
                
                # ######################################
                
            '''
            Perfrom curve fitting: #############  CF 110 peak  #############
            '''
            if fit110 == True:
                x = metaData.RadialDist
                y = dat2[pointNo] #
                # y = np.sum(metaData.Intensity[pointNo], axis = 0)
                
                ############# @@@@@@@@@@@@@ #############
                # plt.figure()
                # plt.plot(x, y)
                #######################################
                
                p110 = 1480
                LowerLim = 1000
                UpperLim = 1530
                
                spec = {
                        'x': x,
                        'y': y,
                        'model': [
                                    {'type': 'PseudoVoigtModel',
                                            'params': {
                                                'center'  : p110,
                                                'sigma'   : 100
                                                },
                                            'help':{
                                                'center'  : {'min': 1300},
                                                # 'sigma'   : {'min': 10, 'max': 150},
                                                }
                                    },
                                    {'type': 'PseudoVoigtModel',
                                            'params': {
                                                'center'  : LowerLim,
                                                    'sigma'   : 250,
                                                },
                                            'help':{
                                                'center'  : {'max': 1300},
                                                # 'sigma'   : {'min': 150},
                                                }
                                    },
                                    
                            ]
                        }
                FittingOutput = metaData.PeakModelGen(spec = spec, LowerLim = LowerLim, UpperLim = UpperLim )
                # FittingOutput.fit_report() # ! uncomment to show fitting report 
                # metaData.PeakResidualPlot(FittingOutput) # ! uncomment to show residual plot
                
                # ** Retrieve and assign values to be stored 
                metaData.ResultArray.loc[pointNo, 'xc_CF_110'] = FittingOutput.params['m0_center'].value
                metaData.ResultArray.loc[pointNo, 'err_xc_CF_110'] = FittingOutput.params['m0_center'].stderr
                metaData.ResultArray.loc[pointNo, 'chi_CF_110'] = FittingOutput.redchi # FittingOutput.chisqr
        #######################################

    metaData.ResultArray.to_csv(outDir + str(scanNo) + output_suffix, na_rep='NaN', index=False)
    t2 = time.perf_counter()

    print(f'Finished in {t2-t1} seconds')
    
    
    
    
    
    
    
    
    
    
    



