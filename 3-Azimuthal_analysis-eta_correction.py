# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 13:03:29 2023

@author: js2580
"""


#import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import lmfit

#import NEXSUS/h5 file reader
import h5py 

#import progressive bar
from tqdm import tqdm

import time


from CFXRD.CFXRD import *


# Define directory paths
ParentDir = 'C:\\Users\\js2580\\OneDrive - University of Bath\\Phd\\Research\\Synchrotron\\XRD\\MM30528 - doping\\WAXS analysis\\Analysis'

outDir = ParentDir + '\\data\\lattice_strain\\'
output_suffix = '_eta.csv'
nexus_dir = ParentDir + '\\data\\caking\\_eigerScan'
angle_dir = ParentDir + '\\data\\fibre_angle\\'


# Define column names
ColumnNames = ['ScanNo', 'PointNo', 'Xmotor', 'Ymotor', 'Cat',\
           'angle1', 'err_angle1','angle_redchi1', 'FWHM_angle1',\
           'angle2', 'err_angle2', 'angle_redchi2', 'FWHM_angle2',\
           'xc_CF_002', 'err_xc_CF_002', 'chi_CF_002', 'FWHM_002', 'err_FWHM_002',\
           'xc_CF_100', 'err_xc_CF_100', 'chi_CF_100', 'FWHM_100', 'err_FWHM_100', \
           'xc_CF_100*', 'err_xc_CF_100*', 'chi_CF_100*',\
           'xc_CF_110', 'err_xc_CF_110', 'chi_CF_110',\
           'xc_CF_004', 'err_xc_CF_004', 'chi_CF_004',\
           'xc_CF_100_Full', 'err_xc_CF_100_Full', 'chi_CF_100_Full',\
           'xc_Al_111_L', 'err_xc_Al_111_L', \
           'xc_Al_200_L', 'err_xc_Al_200_L', \
           'xc_Al_111_T', 'err_xc_Al_111_T', \
           'xc_Al_200_T', 'err_xc_Al_200_T', \
           'xc_Al_220', 'err_xc_Al_220',\
           'xc_Al_111_Full', 'err_xc_Al_111_Full', \
           'xc_Al_200_Full', 'err_xc_Al_200_Full'
           ]

setDry = ['389730_backgroundsubtracted']
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


# scanlist = [389599,389600,389601]
# scanlist = [389633]
scanlist = set5
# scanlist = set5 + set6 + set7 + set8
# scanlist = setDry
doping = False

fit002 = True
fit100 = True
fit100s = False
fit110 = False
fitAl = False

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
    
    #Check if there is an exisiting file
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
    
    # Loop the stack (z-axis)
    for pointNo in tqdm(range(metaData.TotalSlice)):
    # for pointNo in [1312]:
    # for pointNo in range(700,705):
        #print(f'im No: {pointNo}')
        if cat[pointNo] == 'Fibre':
            if np.absolute(1 - angle_chi1[pointNo]) < np.absolute(1 - angle_chi2[pointNo]):
                angle = angle1[pointNo] 
                fwhm = fwhm_angle1[pointNo]
            else:
                angle = angle2[pointNo] 
                fwhm = fwhm_angle2[pointNo]
                
            # index_A, index_invA, index_B, index_invB = metaData.Anisotropic_selective_angle(PeakAngle = angle, fwhmA = 80, fwhmB = 120) # Pristine
            index_A, index_invA, index_B, index_invB = metaData.Anisotropic_selective_angle(PeakAngle = angle, fwhmA = 75, fwhmB = 90)
            
            
                        
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
                if doping == False or doping == True:
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
                    # FittingOutput.fit_report()
                    # metaData.PeakResidualPlot(FittingOutput)
                    
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
                
                if doping == False:
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
                    # LowerLim = 1
                    # UpperLim = 4
                    # spec = {
                    #         'x': x,
                    #         'y': y,
                    #         'model': [
                    #                     {'type': 'PseudoVoigtModel',
                    #                           'params': {
                    #                               'center'  : 3,
                    #                                 # 'sigma'   : 0.1
                    #                               },
                    #                             'help':{
                    #                                 'center'  : {'min': 2.5},
                    #                                 # 'sigma'   : {'max': 2},
                    #                                 }
                    #                     },
                    #                     {'type': 'PseudoVoigtModel',
                    #                             'params': {
                    #                                 'center'  : 1.4,
                    #                                 # 'sigma'   : 0.55,
                    #                                 },
                    #                             'help':{
                    #                                 'center'  : {'min':1.2, 'max': 1.6},
                    #                                 # 'sigma'   : {'max': 1},
                    #                                 }
                    #                     },
                    #                     {'type': 'PseudoVoigtModel',
                    #                             'params': {
                    #                                 'center'  : 2.5,
                    #                                 # 'sigma'   : 0.55,
                    #                                 },
                    #                             'help':{
                    #                                 # 'center'  : {'min': 3.2},
                    #                                 # 'sigma'   : {'min': 0.5},
                    #                                 }
                    #                     },
                                        
                                        
                    #             ]
                    #         }
                  
                elif doping == True:
                    p100 = 3.0
                    LowerLim = 2.6
                    UpperLim = 3.8
                    
                    spec = {
                            'x': x,
                            'y': y,
                            'model': [
                                        {'type': 'GaussianModel',
                                              'params': {
                                                  'center'  : p100,
                                                   'sigma'   : 25
                                                  },
                                                'help':{
                                                    'center'  : {'min':690, 'max':725},
                                                    'sigma'   : {'min': 10, 'max': 80},
                                                    }
                                        },
                                        
                                        {'type': 'PseudoVoigtModel',
                                                'params': {
                                                    'center'  : 725,
                                                      'sigma'   : 200,
                                                    },
                                                'help':{
                                                    # 'center'  : {'max': 370},
                                                    'sigma'   : {'min': 80},
                                                    }
                                        },
                                        # {'type': 'LinearModel',
                                                
                                        # },
                                        
                                        {'type': 'LorentzianModel',
                                                'params': {
                                                    'center'  : 635,
                                                       'sigma'   : 5,
                                                    },
                                                'help':{
                                                    'center'  : {'min':615, 'max': 650},
                                                    # 'sigma'   : {'min': 80},
                                                    }
                                        },
                                        
                                        {'type': 'LorentzianModel',
                                                'params': {
                                                    'center'  : 745,
                                                       'sigma'   : 14,
                                                    },
                                                'help':{
                                                    'center'  : {'min':725, 'max': 760},
                                                    'sigma'   : {'min':5,'max': 30},
                                                    }
                                        },
                                        
                                        
                                ]
                            }

                    
                FittingOutput = metaData.PeakModelGen(spec = spec, LowerLim = LowerLim, UpperLim = UpperLim )
                # FittingOutput.fit_report()
                # metaData.PeakResidualPlot(FittingOutput)
                
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
                
                if doping == True:
                    xc_AL111 = FittingOutput.params['m2_center'].value
                    xcstd_AL111 = FittingOutput.params['m2_center'].stderr
                    # xcchi_AL111 = FittingOutput.chisqr
                    xcchi_AL111 = FittingOutput.redchi
                    
                    metaData.ResultArray.loc[pointNo, 'xc_Al_111_T'] = xc_AL111
                    metaData.ResultArray.loc[pointNo, 'err_xc_Al_111_T'] = xcstd_AL111
                    
                    xc_AL200 = FittingOutput.params['m3_center'].value
                    xcstd_AL200 = FittingOutput.params['m3_center'].stderr
                    # xcchi_AL200 = FittingOutput.chisqr
                    xcchi_AL200 = FittingOutput.redchi
                    
                    metaData.ResultArray.loc[pointNo, 'xc_Al_200_T'] = xc_AL200
                    metaData.ResultArray.loc[pointNo, 'err_xc_Al_200_T'] = xcstd_AL200
                
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
                
                if doping == False:
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
                
                elif doping == True:
                    
                    x = metaData.RadialDist
                    y = allA
                    
                    
                    p100 = 711
                    LowerLim = 610
                    UpperLim = 800
                    
                    spec = {
                            'x': x,
                            'y': y,
                            'model': [
                                        {'type': 'PseudoVoigtModel',
                                              'params': {
                                                  'center'  : p100,
                                                   'sigma'   : 25
                                                  },
                                                'help':{
                                                    'center'  : {'min':690, 'max':725},
                                                    # 'sigma'   : {'min': 10, 'max': 80},
                                                    }
                                        },
                                        
                                        {'type': 'PseudoVoigtModel',
                                                'params': {
                                                    'center'  : 725,
                                                      'sigma'   : 200,
                                                    },
                                                'help':{
                                                    # 'center'  : {'max': 370},
                                                    'sigma'   : {'min': 80},
                                                    }
                                        },
                                        
                                        {'type': 'PseudoVoigtModel',
                                                'params': {
                                                    'center'  : 635,
                                                       'sigma'   : 5,
                                                    },
                                                'help':{
                                                    'center'  : {'min':615, 'max': 650},
                                                    # 'sigma'   : {'min': 80},
                                                    }
                                        },
                                        
                                        {'type': 'PseudoVoigtModel',
                                                'params': {
                                                    'center'  : 745,
                                                       'sigma'   : 14,
                                                    },
                                                'help':{
                                                    'center'  : {'min':725, 'max': 760},
                                                    'sigma'   : {'min':5,'max': 30},
                                                    }
                                        },
                                        
                                        
                                ]
                            }
                    
                    
                FittingOutput = metaData.PeakModelGen(spec = spec, LowerLim = LowerLim, UpperLim = UpperLim )
                # FittingOutput.fit_report()
                # metaData.PeakResidualPlot(FittingOutput)
                
                metaData.ResultArray.loc[pointNo, 'xc_CF_100*'] = FittingOutput.params['m0_center'].value
                metaData.ResultArray.loc[pointNo, 'err_xc_CF_100*'] = FittingOutput.params['m0_center'].stderr
                metaData.ResultArray.loc[pointNo, 'chi_CF_100*'] = FittingOutput.redchi # FittingOutput.chisqr
                
                if doping == True:
                    xc_AL111 = FittingOutput.params['m2_center'].value
                    xcstd_AL111 = FittingOutput.params['m2_center'].stderr
                    # xcchi_AL111 = FittingOutput.chisqr
                    xcchi_AL111 = FittingOutput.redchi
                    
                    metaData.ResultArray.loc[pointNo, 'xc_Al_111_L'] = xc_AL111
                    metaData.ResultArray.loc[pointNo, 'err_xc_Al_111_L'] = xcstd_AL111
                    
                    xc_AL200 = FittingOutput.params['m3_center'].value
                    xcstd_AL200 = FittingOutput.params['m3_center'].stderr
                    # xcchi_AL200 = FittingOutput.chisqr
                    xcchi_AL200 = FittingOutput.redchi
                    
                    metaData.ResultArray.loc[pointNo, 'xc_Al_200_L'] = xc_AL200
                    metaData.ResultArray.loc[pointNo, 'err_xc_Al_200_L'] = xcstd_AL200
                    
                    
                    
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
                
                if doping == False:
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
                    # FittingOutput.fit_report()
                    # metaData.PeakResidualPlot(FittingOutput)
                    
                    
                    metaData.ResultArray.loc[pointNo, 'xc_CF_110'] = FittingOutput.params['m0_center'].value
                    metaData.ResultArray.loc[pointNo, 'err_xc_CF_110'] = FittingOutput.params['m0_center'].stderr
                    metaData.ResultArray.loc[pointNo, 'chi_CF_110'] = FittingOutput.redchi # FittingOutput.chisqr
        #######################################
                              
            '''
            Perfrom curve fitting: #############  Al peaks  #############
            '''
            if fitAl == True:
                metalAngle = 20
                index_C = metaData.MetalDoping(metalAngle)
                
                # plt.figure()
                xc_Al_111 = []
                xcstd_Al_111 = []
                MidAngle = np.arange(metaData.AzimuthalAngle[0] + metalAngle/2, metaData.AzimuthalAngle[0] + metalAngle/2 + (len(index_C)) * metalAngle, metalAngle)
                tmp = []
                tmpstd = []
                for indexSet in index_C:
                # for indexSet in [ [70, 71]]:
                    x = metaData.RadialDist
                    y = np.sum(metaData.Intensity[pointNo][indexSet], axis = 0)
                    #MidAngle.append((metaData.AzimuthalAngle[indexSet[0]] + metaData.AzimuthalAngle[indexSet[-1]])/2)
                    
                    # plt.plot(x,y)
                    
                    LowerLim = 605
                    UpperLim = 800
                       
                    spec = {
                           'x': x,
                           'y': y,
                           'model': [
                                       {'type': 'PseudoVoigtModel',
                                             'params': {
                                                 'center'  : 711,
                                                  'sigma'   : 25
                                                 },
                                               'help':{
                                                   'center'  : {'min':690, 'max':725},
                                                    'sigma'   : {'min': 10, 'max': 80},
                                                   }
                                       },
                                       
                                       # {'type': 'PseudoVoigtModel',
                                       #         'params': {
                                       #             'center'  : 725,
                                       #               'sigma'   : 200,
                                       #             },
                                       #         'help':{
                                       #             # 'center'  : {'max': 370},
                                       #              'sigma'   : {'min': 50, 'max':300},
                                       #             }
                                       # },
                                       {'type': 'LinearModel',
                                       #         'params': {
                                       #             'center'  : 725,
                                       #               'sigma'   : 200,
                                       #             },
                                                # 'help':{
                                                    # 'center'  : {'max': 370},
                                                     # 'intercept'   : {'min': 0},
                                                    # }
                                        
                                               
                                       },
                                       
                                       {'type': 'LorentzianModel',
                                               'params': {
                                                   'center'  : 630,
                                                      'sigma'   : 10,
                                                   },
                                               'help':{
                                                   'center'  : {'min':620, 'max': 640},
                                                    'sigma'   : {'min': 1, 'max':20},
                                                   }
                                       },
                                       
                                       {'type': 'LorentzianModel',
                                               'params': {
                                                   'center'  : 745,
                                                      'sigma'   : 14,
                                                   },
                                               'help':{
                                                   'center'  : {'min':725, 'max': 760},
                                                   'sigma'   : {'max': 20},
                                                   }
                                       },
                                       
                                       
                               ]
                           }
                    FittingOutput = metaData.PeakModelGen(spec = spec, LowerLim = LowerLim, UpperLim = UpperLim )
                    # FittingOutput.fit_report()
                    # metaData.PeakResidualPlot(FittingOutput)
                    
                    xc_Al_111.append(FittingOutput.params['m2_center'].value)
                    xcstd_Al_111.append(FittingOutput.params['m2_center'].stderr)
                    tmp.append(FittingOutput.params['m0_center'].value)
                    tmpstd.append(FittingOutput.params['m0_center'].stderr)
                    
                plt.figure()    
                xcstd_Al_111 = np.array(xcstd_Al_111, dtype=np.float64)
                xc_Al_111 = np.array(xc_Al_111, dtype=np.float64)
                #MidAngle = np.array(MidAngle, dtype=np.float64)
                
                tmp = np.array(tmp, dtype=np.float64)
                tmpstd = np.array(tmpstd, dtype=np.float64)
                
            
                # output , idx_mask = utils.percentilesRemoverNan(xcstd_Al_111)
                # plt.errorbar(MidAngle[idx_mask==1],xc_Al_111[idx_mask==1], yerr=xcstd_Al_111[idx_mask==1], fmt='-x', capsize=10)
                # plt.xlabel('Azimuthal angle ($^\circ$)')
                # plt.ylabel('Peak centre in term of radial distance (pixel)')
                # fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
                # #ax.errorbar(MidAngle[idx_mask==1]*np.pi/180,xc_Al_111[idx_mask==1], yerr=len(xcstd_Al_111[idx_mask==1]) * [3], fmt='-o', capsize= 7)
                # ax.errorbar(MidAngle[idx_mask==1]*np.pi/180,xc_Al_111[idx_mask==1], yerr=xcstd_Al_111[idx_mask==1], fmt='-o', capsize= 1)
                # ax.set_theta_direction(-1)
                # ax.set_rlim(np.min(xc_Al_111[idx_mask==1]) - 0.5, np.max(xc_Al_111[idx_mask==1]) + 0.5)

                
                # fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
                # ax.errorbar(MidAngle[idx_mask==1]*np.pi/180,tmp[idx_mask==1], yerr=tmpstd[idx_mask==1], fmt='-x', capsize=10)
                # ax.set_theta_direction(-1)
                # ax.set_rlim(700,730)


    
    metaData.ResultArray.to_csv(outDir + str(scanNo) + output_suffix, na_rep='NaN', index=False)
    t2 = time.perf_counter()

    print(f'Finished in {t2-t1} seconds')
    
    
    
    
    
    
    
    
    
    
    



