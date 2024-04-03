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
ParentDir = 'C:\\Users\\js2580\\OneDrive - University of Bath\\Phd\Research\\Synchrotron\\XRD\\WAXS analysis\\Analysis\\'

outDir = ParentDir + '\\data\\lattice_strain\\'
output_suffix = '_Dry.csv'
nexus_dir = ParentDir + '\\data\\caking\\800\\_eigerScan'
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

setDry = ['389730']

scanNo = setDry[0]

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
    

pointNo = 0
# Read from excel file
# angle1 = metaData.ResultArray['angle1']
# angle2 = metaData.ResultArray['angle2']

angle = 174


index_A, index_invA, index_B, index_invB = metaData.Anisotropic_selective_angle(PeakAngle = angle, fwhmA = 75, fwhmB = 120) # Pristine

allA = metaData.Integrate_selective_angle_intensity(metaData.Intensity_Azimuthal[pointNo][index_A + index_invA])
allB = metaData.Integrate_selective_angle_intensity(metaData.Intensity_Azimuthal[pointNo][index_B + index_invB])


x = Q
y = allA
# y = allA - allB

############# @@@@@@@@@@@@@ #############
# plt.figure()
# plt.plot(x, y)

# p002 = 1.8
# LowerLim = 1.50 #optimised for allA
# UpperLim = 2.00 #optimised for allA
# # LowerLim = 1.20 #optimised for allA - allB
# # UpperLim = 2.70 #optimised for allA - allB

# spec = {
#         'x': x,
#         'y': y,
#         'model': [
#                     {'type': 'PseudoVoigtModel',
#                             'params': {
#                                 'center'  : p002,
#                           #      # 'sigma'   : 50
#                                 },
#                             'help':{
#                                 # 'center'  : {'min': 1.6},
#                                 # 'sigma'   : {'min': 0.5},
#                                 }
#                     },
#                     # {'type': 'GaussianModel',
#                     #         'params': {
#                     #             'center'  : p002,
#                     #       #      # 'sigma'   : 50
#                     #             },
#                     #         'help':{
#                     #             'center'  : {'min': 1.6},
#                     #             # 'sigma'   : {'min': 0.5},
#                     #             }
#                     # },
#                     # {'type': 'PseudoVoigtModel',
#                     #         'params': {
#                     #             'center'  : 2.0,
#                     #               # 'sigma'   : 1,
#                     #             },
#                     #             'help':{
#                     #             # 'center'  : {'max': 1.6},
#                     #             # 'sigma'   : {'min': 30},
#                     #             }
#                     # },
#                     {'type': 'LinearModel',
                                                  
#                                             },
                    
#             ]
#         }

# FittingOutput = metaData.PeakModelGen(spec = spec, LowerLim = LowerLim, UpperLim = UpperLim )
# FittingOutput.fit_report()
# metaData.PeakResidualPlot(FittingOutput)



x = Q
y = allB

############# @@@@@@@@@@@@@ #############
plt.figure()
plt.plot(x, y)
#######################################

p100 = 3.0
# LowerLim = 2.7
# UpperLim = 3.3
LowerLim = 2.88
UpperLim = 3.20


spec = {
        'x': x,
        'y': y,
        'model': [
                    {'type': 'PseudoVoigtModel',
                          'params': {
                              'center'  : p100,
                                # 'sigma'   : 0.1
                              },
                              'help':{
                                # 'center'  : {'max': 3.2},
                                # 'sigma'   : {'max': 0.4},
                                }
                    },
                    # {'type': 'GaussianModel',
                    #         'params': {
                    #             # 'center'  : 3,
                    #             # 'sigma'   : 0.55,
                    #             },
                    #         'help':{
                    #             # 'center'  : {'min': 3.2},
                    #             # 'sigma'   : {'min': 0.4},
                    #             }
                    # },
                    {'type': 'LinearModel',
                                                                     
                                                                }
                    
            ]
        }

FittingOutput = metaData.PeakModelGen(spec = spec, LowerLim = LowerLim, UpperLim = UpperLim )
FittingOutput.fit_report()
metaData.PeakResidualPlot(FittingOutput)
