# -*- coding: utf-8 -*-
"""
Created on Sun Aug 28 12:08:08 2022

@author: js2580
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.markers as markers
import pandas as pd
import math
from math import sin, cos
from lmfit import models
from lmfit import Parameters
from scipy import ndimage
import warnings

# import h5py 
# import matplotlib.pyplot as plt


from matplotlib import rc
plt.rcParams["font.family"] = "Times New Roman"

import matplotlib.style
import matplotlib as mpl
mpl.rcParams['font.size'] = 14

class CFXRD:
    """This is a class to analyse X-ray diffraction of carbon fibre, containing necessary functions for the analysis.
    """

    def setDetectorParams(self, pixelSize, SD, wavelength = None, BeamEnergy = None):
        """This function stores experimental parameters.

        Args:
            :pixelSize (real): Pixel size in mm
            :SD (real): Sample-to-Detector distnace in mm
            :wavelength (real): wavelegnth in angstorm
            :BeamEnergy (real): Beam energy in keV

        Raises:
            :TypeError: Only wavelength or Beam energy is required
        """
        self.pixelSize = pixelSize
        self.SD = SD
        if (wavelength != None) & (BeamEnergy == None):
            self.wavelength = wavelength
        elif (wavelength == None) & (BeamEnergy != None):
            self.wavelength = 6.62607015E-34 * 3E8 / (BeamEnergy * 1.60218E-16 ) * 10**10 #angstorm
        else:
             raise TypeError("Only wavelength or Beam energy is required")

    def readCakeData(self, Intensity, AzimuthalAngle, RadialDist):
        """This function read input caked data (polar coordinate: Azimuthal angle VS Radial distance)
        Caked data will be 3D array e.g. [No. slice, azimuthal angle, radial distance]

        Min, Max and data range are calculated

        Args:
            :Intensity (real): Intensity
            :AzimuthalAngle (real): Azimuthal angle
            :RadialDist (real): Radial distance
        """
        self.Intensity = Intensity
        self.TotalSlice = len(Intensity)
        
        self.AzimuthalAngle = AzimuthalAngle
        self.AzimuthalRange = [AzimuthalAngle.min(), AzimuthalAngle.max()]
        self.AzimuthalBin = len(AzimuthalAngle)
        self.AzimuthalPerBin = (AzimuthalAngle.max() - AzimuthalAngle.min()) / (self.AzimuthalBin - 1)

        self.RadialDist = RadialDist
        self.RadialRange = [RadialDist.min(), RadialDist.max()]
        self.RadialBin = len(RadialDist)
        self.RadialPerBin = (RadialDist.max() - RadialDist.min()) / (self.RadialBin - 1)

    def cake_RadialIntegration(self, radialMax, radialMin, angleMax = 450, angleMin = 90):
        """This function performs average radial integration from 2D cake data (Azimuthal angle VS Radial distance) into 1D data (Intensity vs Azimuthal angle)
        

        Args:
            :radialMax (real): Maximum value in radial distance
            :radialMin (real): Minimum value in radial distance
            :angleMax (real): Maximum azimuthal angle. Defaults to 450.
            :angleMin (real): Minimum azimuthal angle. Defaults to 90.
        """
        self.Intensity_Radial  = self.Intensity[:, 
                                           np.argmin(np.abs(self.AzimuthalAngle-angleMin)):np.argmin(np.abs(self.AzimuthalAngle-angleMax))+1 ,
                                            np.argmin(np.abs(self.RadialDist-radialMin)):np.argmin(np.abs(self.RadialDist-radialMax))+1]
        self.Intensity_Radial = np.nanmean(self.Intensity_Radial, axis=2)


    def cake_AzimuthalIntegration(self, radialMax, radialMin, angleMax = 450, angleMin = 90):
        """This function performs average azimuthal integration from 2D cake data (Azimuthal angle VS Radial distance) into 1D data (Intensity vs Radial distance, 2θ or q)

        Args:
            :radialMax (real): Maximum value in radial distance
            :radialMin (real): Minimum value in radial distance
            :angleMax (real): Maximum azimuthal angle. Defaults to 450.
            :angleMin (real): Minimum azimuthal angle. Defaults to 90.
        """
        self.Intensity_Azimuthal  = self.Intensity[:, 
                                        np.argmin(np.abs(self.AzimuthalAngle-angleMin)):np.argmin(np.abs(self.AzimuthalAngle-angleMax))+1 ,
                                            np.argmin(np.abs(self.RadialDist-radialMin)):np.argmin(np.abs(self.RadialDist-radialMax))+1]
    
    
    def GenResultArray(self, ColumnNames = None):
        """This function generate a pre-define table data to stored neccessary parameters

        Args:
            :ColumnNames (list): List of header names needs to be stored
        """
        NumColumn = len(ColumnNames)
        temp = np.zeros((self.TotalSlice, NumColumn), dtype=object)
        temp[:] = np.nan
        self.ResultArray = temp
        if ColumnNames  == None:
            ColumnNames = [None] * self.TotalSlice
        else:
            self.ColumnNames = ColumnNames
        self.ResultArray = pd.DataFrame(self.ResultArray, columns = ColumnNames)
    #
    def Anisotropic_selective_angle(self, PeakAngle, fwhmA, fwhmB):
        """This function determines index angles from the centre of the peak to cover the defined angles.        
        Index is adjusted when goes beyond the maximum index range.

        Args:
            :PeakAngle (real): Angle at the centre of the peak
            :fwhmA (real): Full-Width-Half-Maximum (002) or Angle width
            :fwhmB (real): Full-Width-Half-Maximum (100) or Angle width

        Returns:
            :indexes (list): List of indexes in first/second peaks of 002 and 100
        
        .. note::
            :index_A = 002 angle
            :index_invA = invert 002 (180 degrees offset or 180/anglePerbin)
            :index_B = 100 angle (90 degrees offset to 002)
            :index_invB = invert 100 (180 degrees offset or 180/anglePerbin)
        """
        selectIndex = np.absolute(self.AzimuthalAngle - PeakAngle).argmin()
        offset90 = np.absolute(self.AzimuthalAngle - (PeakAngle + 90)).argmin()
    
        # Change Full widht half maximum into bin width from the peak centre in term of index
        BinWidthA = math.ceil((fwhmA/2)/self.AzimuthalPerBin) 
        BinWidthB = math.ceil((fwhmB/2)/self.AzimuthalPerBin) 
        #

        #
        index_A = list()
        index_invA = list()
        index_B = list()
        index_invB = list()
        #
        # offset from the centre peak by BinWidth
        index_A = list(
            np.linspace( 
                round(selectIndex)-BinWidthA, 
                round(selectIndex)+BinWidthA, 
                BinWidthA*2 + 1, dtype=(int) ) 
            )
        #
        index_invA = list(
            np.linspace( 
                round(selectIndex + math.ceil(180/self.AzimuthalPerBin))-BinWidthA, 
                round(selectIndex + math.ceil(180/self.AzimuthalPerBin))+BinWidthA, 
                BinWidthA*2 + 1, dtype=(int) ) 
            )
        #
        index_B = list(
            np.linspace( 
                round(offset90)-BinWidthB, 
                round(offset90)+BinWidthB, 
                BinWidthB*2 + 1, dtype=(int) ) 
            )
        #
        index_invB = list(
            np.linspace( 
                round(offset90 + math.ceil(180/self.AzimuthalPerBin))-BinWidthB, 
                round(offset90 + math.ceil(180/self.AzimuthalPerBin))+BinWidthB, 
                BinWidthB*2 + 1, dtype=(int) ) 
            )
        def AdjustIndex(index):
            # Adjust index that goes beyond the maximum index range 
            # Example [359,360,361,362] --> [359, 0, 1 ,2] where 359 is the maximum
            if any(i >= self.AzimuthalBin for i in index):
                idx = np.array(list(map(int, [number >= self.AzimuthalBin for number in index])))
                index = np.array(index)
                index[idx == 1] = index[idx == 1] - self.AzimuthalBin
            else:
                pass
            #
            return list(index)
        #
        index_A = AdjustIndex(index_A)
        index_invA = AdjustIndex(index_invA)
        #
        index_B = AdjustIndex(index_B)
        index_invB = AdjustIndex(index_invB)
        #

        return index_A, index_invA, index_B, index_invB
    #
    def Filter_selective_angle(self, indexes, pointNo, lowThreshold, upThreshold):

        """Filter selective angle index where it consists of missing data due to the detector stripes (dead pixel)
            lower and upper Thresholds set acccording the missing data (excluding mask data at beam centre)

        Args:
            :indexes (integer): _description_
            :pointNo (integer): _description_
            :lowThreshold (real): _description_
            :upThreshold (real): _description_

        Returns:
            :indexes: indexes
        """
        pop_list = []
        for k, index in enumerate(indexes):
            argmin = min(np.where(self.Intensity[pointNo][index] != 0)[0])  #Determine the first non-zero data after beam centre mask
            ydiff = np.diff(self.Intensity[pointNo][index][argmin:])
            if min(ydiff) <= lowThreshold or max(ydiff) >= upThreshold:
                indexes.pop(np.where(np.array(indexes) == index)[0][0]) # Update new indexes
        return indexes 
    #
    def Integrate_selective_angle_intensity(self, array):
        """This function average the intensities over the selective angle or indexes
            At the point where intensity = 0, it is masked out from the averaging calculation
            It then return the NAN value back to zero for further curve fitting

        Args:
            :array (real): _description_

        Returns:
            :average: average value ignore NAN value
        """
        array_with_nan = np.where(array == 0, np.nan, array)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            average_with_nan = np.nanmean(array_with_nan, axis=0) #average NAN over masked area is expected runtime warning
        average = np.where(np.isnan(average_with_nan), 0, average_with_nan)
        
        return average
    #

    def PeakModelGen(self, spec, LowerLim = None, UpperLim = None):

        """This function generates composite fitting model, see more detail: https://chrisostrouchov.com/post/peak_fit_xrd_python/

        Args:
            :spec (dict): models specifications
            :LowerLim (real): Lower limit x-axis
            :UpperLim (real): Upper limit x-axis

        Raises:
            :NotImplementedError: Only 'GaussianModel', 'LorentzianModel', 'VoigtModel' are available

        Returns:
            :FittingOutput: Fitting output
        """
        # .. note::
        #     spec = {
        #     'x': x,
        #     'y': y,
        #     'model': [
        #                 {'type': 'PseudoVoigtModel',
        #                     'params': {
        #                         'center'  : 1.5,
        #                         'sigma'   : 2,
        #                         'fraction': 0.5
        #                         },
        #                     'help':{
        #                         'center'  : {'min': -3, 'max': 3},
        #                         'sigma'   : {'min': 1E-6, 'max': 6},
        #                         'fraction': {'min': 0, 'max': 0.5}
        #                         }
        #                 },
        #                 {'type': 'LinearModel',
        #                     'params': {
        #                         'slope'       : 2,
        #                         'intercept'   : 2,
        #                         },
        #                     'help':{
        #                         'slope'       : {'min': -3, 'max': 3},
        #                         'intercept'   : {'min': 1E-6, 'max': 6},
        #                         }
        #                 },
                        
        #         ]
        #     }
            
        # .. warning::   
        #     LowerLim and UpperLim is the actual value on the x-axis that you want to limit on
        
        composite_model = None
        params = None
        x = spec['x']
        y = spec['y']

        if (LowerLim != None) & (UpperLim != None):
            # Pre-defined indexMask with zeros
            indexMask = np.zeros(len(x), dtype=int)
            # if Datatype == 'AzimuthalData':
                # head = int(LowerLim/self.RadialPerBin)
                # tail = int(UpperLim/self.RadialPerBin)
                
            head = np.absolute(x - LowerLim).argmin()
            tail = np.absolute(x - UpperLim).argmin() 
                
            indexMask[head:tail+1] = 1
            x = x[indexMask == 1]
            y = y[indexMask == 1]
            # Update spec x and y with bounds
            spec['x'] = x 
            spec['y'] = y 
        #
        x_min = np.min(x)
        x_max = np.max(x)
        x_range = x_max - x_min
        y_max = np.max(y)
        self.spec = spec
        #
        for i, basis_func in enumerate(spec['model']):
            prefix = f'm{i}_'
            model = getattr(models, basis_func['type'])(prefix=prefix)  # get a model from specified function from lmlit
            if basis_func['type'] in ['GaussianModel', 'LorentzianModel', 'VoigtModel']: # for now VoigtModel has gamma constrained to sigma
                model.set_param_hint('sigma', min=1e-6, max=float(x_range))
                model.set_param_hint('center', min=float(x_min), max=float(x_max))
                model.set_param_hint('height', min=1e-6, max=float(1.1*y_max))
                model.set_param_hint('amplitude', min=1e-6)
                # default guess is horrible!! do not use guess()
                # default_params = {
                    # prefix+'center': x_min + x_range * random.random(),
                    # prefix+'height': y_max * random.random(),
                    # prefix+'sigma': x_range * random.random()
                    # }
                default_params = Parameters()
                default_params.add(f'{prefix}center', float(x[np.absolute(y - y_max).argmin()]))
                default_params.add(f'{prefix}height', float(y_max))
                default_params.add(f'{prefix}sigma', float((x_max - x_min)/2))
                
            elif basis_func['type'] in ['PseudoVoigtModel']: 
                model.set_param_hint('sigma', min=1e-6, max=float(x_range)) #max=x_range
                model.set_param_hint('center', min=x_min, max=float(x_max))
                model.set_param_hint('height', min=1e-6, max=float(1.1*y_max))
                model.set_param_hint('amplitude', min=1e-6)
                model.set_param_hint('fraction', min=1e-6, max=1.0)
                # default_params = {
                    # prefix+'center': float(x_min + x_range * random.random()),
                    # prefix+'height': float(y_max * random.random()),
                    # prefix+'sigma': float(x_range * random.random()),
                    # prefix+'fraction': 0.5
                    # }
                # fraction = 1 is full LorentzianModel and fraction = 0 is full GaussianModel
                # default_params = {
                #     prefix+'center': float(x[np.absolute(y - y_max).argmin()]),
                #     prefix+'height': float(y_max),
                #     prefix+'sigma' : float((x_max - x_min)/2),
                #     prefix+'fraction': 0.5
                #     }
                default_params = Parameters()
                default_params.add(f'{prefix}center', float(x[np.absolute(y - y_max).argmin()]))
                default_params.add(f'{prefix}height', float(y_max))
                default_params.add(f'{prefix}sigma', float((x_max - x_min)/2))
                default_params.add(f'{prefix}fraction', 0.5)
                    
            elif basis_func['type'] in ['LinearModel']: 
                model.set_param_hint('slope', min= float(-np.inf), max= float(np.inf)) #max=x_range
                model.set_param_hint('intercept', min= float(-np.inf), max= float(np.inf))
                default_params = Parameters()
                default_params.add(f'{prefix}slope', float((y[-1] - y[0]) / (x[-1] - x[0])))
                default_params.add(f'{prefix}intercept', float(y[0]))
                
            else:
                raise NotImplementedError(f'model {basis_func["type"]} not implemented yet. Only GaussianModel, LorentzianModel, VoigtModel, LinearModel are available')
            #
            if 'help' in basis_func:  # allow override of settings in parameter with specified values
                for param, options in basis_func['help'].items():
                    model.set_param_hint(param, **options)
            #
            model_params = model.make_params(**default_params, **basis_func.get('params', {}))
            #
            if params is None:
                params = model_params   # params for the first peak
            else:
                params.update(model_params) # update params for the following peak
                
            if composite_model is None: # model for the first peak
                composite_model = model
            else:
                composite_model = composite_model + model # combine models
        #
        FittingOutput = composite_model.fit(spec['y'], params, x= spec['x'], method = 'least_squares', fit_kws={'xtol': 2e-10, 'ftol': 2e-10}, max_nfev=1e20)
        """""
        .. code-block:: python
            #Test script

            #import matplotlib for testing
            import matplotlib.pyplot as plt
            import math

            def g(x, A, μ, σ):
                return A / (σ * math.sqrt(2 * math.pi)) * np.exp(-(x-μ)**2 / (2*σ**2))

            x = np.linspace(-3, 3, 1000)
            y =  g(x, 30, 0, 1) + (2*x +1.2)
            fig, ax = plt.subplots()
            ax.plot(x, y)

            spec = {
                    'x': x,
                    'y': y,
                    'model': [
                                {'type': 'PseudoVoigtModel',
                                    'params': {
                                        'center'  : 1.5,
                                        'sigma'   : 2,
                                        'fraction': 0.5
                                        },
                                    'help':{
                                        'center'  : {'min': -3, 'max': 3},
                                        'sigma'   : {'min': 1E-6, 'max': 6},
                                        'fraction': {'min': 0, 'max': 0.5}
                                        }
                                },
                                {'type': 'LinearModel',
                                    'params': {
                                        'slope'       : 2,
                                        'intercept'   : 2,
                                        },
                                    'help':{
                                        'slope'       : {'min': -3, 'max': 3},
                                        'intercept'   : {'min': 1E-6, 'max': 6},
                                        }
                                },
                                
                        ]
                    }

            model, params = PeakModelGen(spec)
            output = model.fit(spec['y'], params, x=spec['x'])
            output.plot()

            components = output.eval_components(x=spec['x'])

            for i, model in enumerate(spec['model']):
                prefix = f'm{i}_'
                ax.plot(spec['x'], components[prefix])
        """
        return FittingOutput
    
    
    def PeakResidualPlot(self, FittingOutput):
        """This function plot peak fitting residual

        Args:
            :FittingOutput (dict): Fitting output; requires PeakModelGen() first
        """
        FittingOutput.plot() # Function in lmfit library to plot residual graph

        ###############################################################
        # #Plot residual graph (alternative way) allows to change labels
        # best_fit = FittingOutput.best_fit
        # x = 2*np.pi/self.spec['x']  
        # raw = self.spec['y']
        # residual = FittingOutput.residual
        # # Create a figure with two subplots
        # fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

        # # Plot the main data and the fit on the first subplot
        # ax1.plot(x, raw, 'o', color = 'blue', label = 'Raw data')
        # ax1.plot(x, best_fit, '-', color = 'red', label='Best fit', linewidth = 3)
        # # ax1.plot(x, model(fit_result.params, x), 'r-', label='Fit')
        # ax1.set_ylabel('Intensity - Counts')
        # ax1.legend()

        # # Plot the residuals on the second subplot
        # ax2.axhline(0, color='black', linestyle='--', linewidth=0.8)
        # ax2.scatter(x, residual, marker='o', color='blue')
        # ax2.set_xlabel('q-spacing - Å$^{-1}$')
        # ax2.set_ylabel('Residuals')
        # ax2.grid(True)

        # # Adjust layout and show the plot
        # plt.tight_layout()
        ###############################################################
        
        fig, ax = plt.subplots()
        ax.scatter(self.spec['x'], self.spec['y'], s=4)
        components = FittingOutput.eval_components(x=self.spec['x'])
        tmp = ['raw data']

        for i, model in enumerate(self.spec['model']):
            prefix = f'm{i}_'
            ax.plot(self.spec['x'], components[prefix])
            tmp.append(prefix)
        ax.legend(tmp)

        plt.xlabel("q-spacing - Å$^{-1}$")
        plt.ylabel("Intensity - Counts")


        # components = self.output.eval_components(x = self.output.userkws['x'])
        # print(components)
        # fig, ax = plt.subplots()
        # ax.plot(self.output.userkws['x'], self.output.data)

        # for i, model in enumerate(self.spec['model']):
        #     prefix = f'm{i}_'
        #     ax.plot(self.output.userkws['x'], components[prefix])
        return

    def motor_position(self, motor_dir, motorX, motorY, scanNo, pointNo):
        """This function reads motor positions input file and are adjusted for visuallisation

        Args:
            :motor_dir (str): Motor directory
            :motorX (str): Motor name in x-axis, e.g. gtbX3
            :motorY (str): Motor name in y-axis, e.g. gtbY2
            :scanNo (int): Scan number (unique number)
            :pointNo (int): Point/Slice number

        Returns:
            :x_pos: X position
            :y_pos: Y position
        """

        #Read data file
        with open(motor_dir + str(scanNo) + ".dat", "r") as file:
            for l_no, line in enumerate(file):
                # search string
                try:
                    #try to replace \t (indentation) to space
                    line = line.replace("\t"," ")
                except:
                    pass
                if f'{motorY} {motorX}' in line:
                    headerline = l_no
                    break # don't look for next lines
        # Read into DataFrame
        data = pd.read_csv(motor_dir + str(scanNo) + ".dat", delimiter= r"\s+", skiprows = headerline, skip_blank_lines=False)
        
        
        x_pos = round(data[motorX].iloc[pointNo], 2)
        y_pos = round(data[motorY].iloc[pointNo], 2)
        
        #### To adjust motor position for visualisation 
        if (x_pos * 10) % 2 != 0.0:
            x_pos -= 0.1
        if (y_pos * 10) % 2 != 0.0:
            y_pos -= 0.1
        
        return x_pos, y_pos
    
    #
    
    def combineDataintoGridarray(Dir:str, scanNo:list, filetype = '.csv'):
        """This function combine multiples scans into a single array with the unique positions
            Motor position can be sorted 
            df = df.sort_values(['Ymotor', 'Xmotor'], ascending = [True, True])

        Args:
            :Dir (str): File directory
            :scanNo (list): Scan number
            :filetype (str): Extension file type. Defaults to '.csv'.

        Returns:
            :girdDat: Array grid data
        """
        #Read the first scan to retrieve columns name first
        dat = pd.read_csv(Dir +  str(scanNo[0]) +  filetype) 
        columns = dat.columns
        df = pd.DataFrame(columns = columns)
        #Iterate through scanlist no. and combine them together
        for scan in scanNo:
            dat = pd.read_csv(Dir +  str(scan) +  filetype)  
            df = pd.concat([df, dat])
        
        df['Ymotor'] = np.round(df['Ymotor'], 2)
        df['Xmotor'] = np.round(df['Xmotor'], 2)
        
        df = df.drop_duplicates(subset=['Ymotor', 'Xmotor'], keep='last')
        gb = df.groupby('Ymotor')    
        gb = [gb.get_group(x) for x in gb.groups]
        
        Ymotor = np.array(df['Ymotor'])
        Xmotor = np.array(df['Xmotor'])
        
        uniqueY = list(np.unique(Ymotor))
        uniqueX = list(np.unique(Xmotor))
        
        noRow = len(np.unique(Ymotor))
        noCol = len(uniqueX)
        gridSize = (noRow,noCol)
            
        gridDat = np.ndarray((noRow,noCol,len(columns)), dtype=object)
        gridDat[:] = np.nan
        
        for row in range(0, noRow):
            if len(gb[row]) != len(uniqueX):
                templist = list(gb[row]['Xmotor'])
                idx = []
                for num in templist:
                    idx.append(uniqueX.index(num))
                dummy = np.empty((len(uniqueX), len(columns)))
                dummy[:] = np.nan
                                
                subdf = pd.DataFrame(dummy, columns=columns)
                
                for i, index in enumerate(idx):
        
                    subdf.iloc[[index]] = gb[row].iloc[i]
                    
                subdf['Xmotor'] = uniqueX
                subdf['Ymotor'] = [uniqueY[row]] * len(uniqueX)
            
                gb[row] = subdf
        
        df = pd.concat(gb, axis=0)
        # Invert axis by chaning ascending
        df = df.sort_values(['Ymotor', 'Xmotor'], ascending = [True, True])
        
        
        for z in range(0,len(columns)):
            gridDat[:,:,z] = np.reshape(df.iloc[:,z].to_numpy(), gridSize)
            pass
        
        return gridDat
    
    def pixel_to_Dspacing(self, pixel):
        """Convert pixel to d-spacing

        Args:
            :pixel (real): pixel value

        Returns:
            :twoTheta: 2θ in degree
            :d: d-spacing in angstorm

        .. warning::
            Requires Pixel size, Wavelength and Sample to detector distance.
            Define these via setDetectorParams().
        """
        pixel = np.array(pixel, dtype=float)
        pixel[pixel == float(0)] = np.nan
        mask = np.isfinite(pixel, dtype=bool)
        
        d = np.zeros(np.shape(pixel))
        d[:] = np.nan
        
        twoTheta = np.zeros(np.shape(pixel))
        twoTheta[:] = np.nan
        
        twoTheta[mask] = np.arctan(pixel[mask]*self.pixelSize/self.SD)
        twoTheta[mask] = twoTheta[mask] * 180 / np.pi
        
        d[mask] = self.wavelength/(2*np.sin(np.radians(twoTheta[mask]/2)))
        
        return twoTheta, d
    
    def Qspacing_to_Dspacing(self, q):
        """Convert q-spacing to d-spacing

        Args:
            :q (real): q-spacing

        Returns:
            :twoTheta: 2θ in degree
            :d: d-spacing in angstorm

        .. warning::
            Requires Wavelength.
            Define this via setDetectorParams().
        """
        q = np.array(q, dtype=float)
        q[q == float(0)] = np.nan
        mask = np.isfinite(q, dtype=bool)
        
        d = np.zeros(np.shape(q))
        d[:] = np.nan
        
        twoTheta = np.zeros(np.shape(q))
        twoTheta[:] = np.nan


        d[mask] = (2*np.pi)/(q[mask])
        
        twoTheta[mask] = 2*np.arcsin(self.wavelength/(2*d[mask]))
        twoTheta[mask] = twoTheta[mask] * 180 / np.pi

        return twoTheta, d
    
    def twoTheta_to_Qspacing(self, twoTheta):
        """Convert 2θ to q-spacing

        Args:
            :twoTheta (real): 2θ in degree

        Returns:
            :q: q-spacing in angstorm^-1

        .. warning::
            Requires Wavelength.
            Define this via setDetectorParams().
        """
        #Input = 2theta in degrees
        twoTheta = np.array(twoTheta, dtype=float)
        twoTheta[twoTheta == float(0)] = np.nan
        mask = np.isfinite(twoTheta, dtype=bool)
        
        q = np.zeros(np.shape(twoTheta))
        q[:] = np.nan
        q[mask]= 4*np.pi/self.wavelength*np.sin(np.radians(twoTheta[mask]/2))

        return q
    
    def twoTheta_to_Dspacing(self, twoTheta):
        """Convert 2θ to d-spacing

        Args:
            :twoTheta (real): 2θ in degree

        Returns:
            :d: d-spacing in angstorm

        .. warning::
            Requires Wavelength.
            Define this via setDetectorParams().
        """
        #2theta in degrees
        twoTheta = np.array(twoTheta, dtype=float)
        twoTheta[twoTheta == float(0)] = np.nan
        mask = np.isfinite(twoTheta, dtype=bool)
        
        d = np.zeros(np.shape(twoTheta))
        d[:] = np.nan
        d[mask]= self.wavelength/ 2 / np.sin(np.radians(twoTheta[mask]/2))

        return d
    
    def Dspacing_to_twoTheta(self, d):
        """Convert d-spacing to 2θ

        Args:
            :d: d-spacing in angstorm

        Returns:
            :twoTheta (real): 2θ in degree

        .. warning::
            Requires Wavelength.
            Define this via setDetectorParams().
        """
        #Input = 2theta in degrees
        d = np.array(d, dtype=float)
        d[d == float(0)] = np.nan
        mask = np.isfinite(d, dtype=bool)
        
        twoTheta = np.zeros(np.shape(d))
        twoTheta[:] = np.nan
        twoTheta[mask]= np.degrees(np.arcsin(self.wavelength/ 2 / d))*2

        return twoTheta
    
    def lattice_strain_d_spacing(self, d, d_0):
        """Convert d-spacing to lattice strain

        Args:
            :d (real): Observed d-spacing
            :d_0 (real): Strain free d-spacing

        Returns:
            :epsilon (real): Lattice strain
        """

        epsilon = (d - d_0)/d_0
        
        return epsilon
    
    def lattice_strain_q_spacing(self, q, q_0):
        """Convert q-spacing to lattice strain

        Args:
            :q (real): Observed q-spacing
            :q_0 (real): Strain free q-spacing

        Returns:
            :epsilon (real): Lattice strain
        """

        epsilon = (q_0 - q)/q
        
        return epsilon
    
    def Epsilon_cal(self, deform:np.array, ref, dtype):
        """Calculate lattice strain
            Outliner is removed via Interquartile criterion

        Args:
            :deform (np.array): Strain values
            :ref (real): Strain free values
            :dtype (str): 'd-spacing', 'q-spacing'

        Raises:
            :ValueError: Invalid data type

        Returns:
            :epsilon (real): Lattice strain
        """
        dtypes = ['d-spacing', 'q-spacing']

        if dtype not in dtypes:
            raise ValueError("Invalid sim type. Expected one of: %s" % dtypes)
        if dtype == 'd-spacing':
            epsilon = self.lattice_strain_d_spacing(deform, ref)
        elif dtype == 'q-spacing':
            epsilon = self.lattice_strain_q_spacing(deform, ref)

        epsilon, _ = utils.InterquartileRemover(epsilon)
        # epsilon = utils.filter_nan_gaussian_conserving(epsilon, sigma=4)
        # epsilon = filter_nan_gaussian_conserving2(epsilon, sigma=1)
        # epsilon = utils.filter_nan_gaussian_david(epsilon, sigma=1)

        return epsilon

    def twoThetaShift(self, disp, twoTheta):
        """Specimen-displacement correction for powder X-ray diffraction in Debye-Scherrer geometry with a flat area detector, see https://doi.org/10.1107/S1600576722011360

        Args:
            :disp (real): Sample Displacement. Defaults to 3. in mm.
            :SD (real): Sample to Detector distance. Defaults to 127.032.
            :twoTheta (np.array): Range of 2θ in degree. Defaults to np.linspace(7, 43, 43-7+1).

        Returns:
            :eta: Correction angle
        """
        numerator = disp*np.sin(np.radians(twoTheta * 2))
        denominator = 2*(self.SD - disp * np.sin(np.radians(twoTheta))**2)
        eta = np.degrees(np.arctan(numerator/denominator))
        # new 2theta = 2theta - eta
        return eta
    
    def scherrer(self, twoTheta, fwhm, k):
        """Scherrer equation for crystal size calculation

        Args:
            twoTheta (real): Twotheta in degree
            fwhm (real): Full-width half maximum in degree
            k (real): Scherrer constant. 0.90 for graphitic layers stacked perpendicular to fibre axis. 1.84 for graphitic layer in the plane along fibre axis.
        """
        fwhm = np.array(fwhm, dtype=float)
        twoTheta = np.array(twoTheta, dtype=float)
        L = (k*self.wavelength)/(np.radians(fwhm)*np.cos(np.radians(twoTheta/2)))
        return L
        
    def scherrer_q(self, fwhm, k):
        """Scherrer equation for crystal size calculation in q-spacing

        Args:
            fwhm (real): Full-width half maximum in Å
            k (real): Scherrer constant. 0.90 for graphitic layers stacked perpendicular to fibre axis. 1.84 for graphitic layer in the plane along fibre axis.
            
        Returns:
        :L: Crystal size in Å
        """
        L = k * 2 * np.pi/fwhm
        return L

    def get_neighbor_average(self, Input:np.array, Cat): 
        """Allow value to propagate to NaN by averaing neighbor values
        """
        def find_neighbors(i, j):
            value_list = []
            sumvalue = 0
            div = 0
            k0 = [0]
            k1 = [-1, 0, 1]
            k2 = [-2, -1, 0, 1, 2]
            k = k1
            # i: row
            # j: column
            for offset_i in k: 
                for offset_j in k:
                    new_i = i + offset_i
                    new_j = j + offset_j
                    if (new_i >= 0 and new_j >= 0 and new_i < Input.shape[0] and new_j < Input.shape[1]):
                        value_list.append(Input[new_i][new_j])
            avg = np.nanmean(value_list)
            err_avg =  np.nanstd(value_list) 
           
            return avg, err_avg
        
        for m in range(Cat.shape[0]):
            for n in range(Cat.shape[1]):
                if (Cat[m, n] == 'Fibre') and (np.isnan(Input[m,n])):
                    Input[m,n], _ = find_neighbors(m,n)

        return Input
    
    def Mapping_Plot(Input: np.array, cbarTitle: str, category: np.array, cbarLabel: str,
                 cbarMax: float = None, cbarMin: float = None, Marker: str = 'OFF',
                 label: str = 'OFF', sci_exp: int = -3):
        """ This function plots heating mapping of the input 2D array like image pixels

        Args:
            Input (np.array): Input 2D array
            cbarTitle (str): Colour bar title
            category (np.array): Category to classify data points such Fibre or Off-sample positions
            cbarLabel (str): Colour bar label
            cbarMax (float): Colour bar max intensity value
            cbarMin (float): Colour bar min intensity value
            Marker (str, optional): Error marker point (x) Defaults to 'OFF'.
            label (str, optional): Pixel number label. Defaults to 'OFF'.

        Raises:
            ValueError: ON or OFF label plot
        """
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
                    
        # Use numpy nanmax/nanmin for default values if not provided
        if cbarMax is None:
            cbarMax = np.nanmax(Input)
        if cbarMin is None:
            cbarMin = np.nanmin(Input)
            
        Cat = category
        height, width = Input.shape

        fig, ax = plt.subplots()
        marker = markers.MarkerStyle('x')

        # Overlay the marker on the image
        if Marker.lower() == 'on':
            for i in range(height):
                for j in range(width):
                    if (Cat[i, j] == 'Fibre') and (np.isnan(Input[i,j])):
                        ax.plot(j, i, marker=marker, color='k', markersize=6, linewidth=12)
        elif Marker.lower() == 'off':
            pass
        else:
            raise ValueError('Please specify ON or OFF label plot')

        # plt.axis('off')
        plt.set_cmap('jet') # plt.set_cmap('rainbow') # plt.set_cmap('tab20c') 
        plot = plt.imshow(Input, aspect=('equal'), origin='lower') #
        plt.clim(cbarMin, cbarMax)

        cbar = fig.colorbar(plot, format=OOMFormatter(sci_exp, mathText=False))
        cbar.set_label(cbarLabel)
        plt.title(cbarTitle)
        plt.tight_layout()    

        if label.lower() == 'on':
            # Get the dimensions of the data
            # Iterate over each pixel and annotate
            index = 0 
            for i in range(height):
                for j in range(width):
                    plt.text(j, i, index, ha='center', fontsize=8, va='center', color='black')
                    index += 1
        elif label.lower() == 'off':
            pass
        else:
            raise ValueError('Please specify ON or OFF label plot')



#
class utils:
    def InterquartileRemover(input):
        # Change to array format
        input = np.array(input)
        
        q75, q25 = np.nanpercentile(input,[75,25])
        intr_qr = q75-q25
        
        uplim = q75+(1.5*intr_qr)
        lowlim = q25-(1.5*intr_qr)

        idx_mask = np.ones_like(input)

        idx_mask[input > uplim] = 0
        idx_mask[input < lowlim] = 0

        output = np.zeros_like(input)
        output[:] = np.nan
        output[idx_mask==1] = input[idx_mask==1]

        return output, idx_mask
    
    def filter_nan_gaussian_david(Input:np.array, sigma=1):
        """Allows intensity to leak into the nan area.
        According to Davids answer: https://stackoverflow.com/a/36307291/7128154
        """
        gauss = Input.copy()
        gauss[np.isnan(gauss)] = 0
        gauss = ndimage.gaussian_filter(
                gauss, sigma=sigma, mode='constant', cval=0)

        norm = np.ones(shape=Input.shape)
        norm[np.isnan(Input)] = 0
        norm = ndimage.gaussian_filter(
                norm, sigma=sigma, mode='constant', cval=0)

        # avoid RuntimeWarning: invalid value encountered in true_divide
        norm = np.where(norm==0, 1, norm)
        gauss = gauss/norm
        gauss[np.isnan(Input)] = np.nan

        return gauss
    
    
    def filter_nan_gaussian_conserving(Input:np.array, sigma):
        """Apply a gaussian filter to an array with nans.

        Intensity is only shifted between not-nan pixels and is hence conserved.
        The intensity redistribution with respect to each single point
        is done by the weights of available pixels according
        to a gaussian distribution.
        All nans in arr, stay nans in gauss.
        """
        nan_msk = np.isnan(Input)

        loss = np.zeros(Input.shape)
        loss[nan_msk] = 1
        loss = ndimage.gaussian_filter(
                loss, sigma=sigma, mode='constant', cval=1)

        gauss = Input.copy()
        gauss[nan_msk] = 0
        gauss = ndimage.gaussian_filter(
                gauss, sigma=sigma, mode='constant', cval=0)
        gauss[nan_msk] = np.nan

        gauss += loss * Input

        return gauss

    def filter_nan_gaussian_conserving2(arr, sigma):
        """Apply a gaussian filter to an array with nans.

        Intensity is only shifted between not-nan pixels and is hence conserved.
        The intensity redistribution with respect to each single point
        is done by the weights of available pixels according
        to a gaussian distribution.
        All nans in arr, stay nans in gauss.
        """
        nan_msk = np.isnan(arr)

        loss = np.zeros(arr.shape)
        loss[nan_msk] = 1
        loss = ndimage.gaussian_filter(
                loss, sigma=sigma, mode='constant', cval=1)

        gauss = arr / (1-loss)
        gauss[nan_msk] = 0
        gauss = ndimage.gaussian_filter(
                gauss, sigma=sigma, mode='constant', cval=0)
        gauss[nan_msk] = np.nan

        return gauss
    

class FibrePlot:
    """Class to plot fibre orientation
    """

    def rotate(self, origin:list, point:list, angle:list):
        """Perform line transformation according to the fibre orientation
            Rotate a point counterclockwise by a given angle around a given origin.
            The angle should be given in radians.

        Args:
            :origin (list): Original/Centre points [x1, y1]
            :point (list): Offset points [x2, y2]
            :angle (list): Angle in radians

        Returns:
            :(x, y): new coordinates
        """

        x = (origin[0] + cos(angle) * (point[0] - origin[0]) - sin(angle) * (point[1] - origin[1]))
        y = origin[1] + sin(angle) * (point[0] - origin[0]) + cos(angle) * (point[1] - origin[1])
        return (x, y)
        
        
    def motorPosition(self,x_pos,y_pos,cat,angle1, angle2, redchi1, redchi2):
        self.x_pos = round(x_pos, 2) * 1
        self.y_pos = round(y_pos, 2) * 1
        
        self.xUnique = np.shape(np.unique(x_pos))[0]
        self.yUnique = np.shape(np.unique(y_pos))[0]
        
        self.xIncrement = round((np.max(self.x_pos) - np.min(self.x_pos))/self.xUnique, 2)
        self.yIncrement = round((np.max(self.y_pos) - np.min(self.y_pos))/self.yUnique, 2)
        
        self.cat = cat
        self.totalPoint = len(self.cat)
        
        self.angle1 = angle1
        self.angle2 = angle2
        
        self.redchi1 = redchi1
        self.redchi2 = redchi2

        
    def Orientation_Plot(self, label:str='OFF'):
        # def forceAspect(ax,aspect=1):
        #     im = ax.get_images()
        #     extent =  im[0].get_extent()
        #     ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
         
        fig = plt.figure()
        plt.axis('equal')
        # plt.axis('off')
        fig.set_size_inches(10, 10)

        k = 0.08
        for i in range(0, self.totalPoint):
            x1 = self.x_pos[i]
            x2 = self.x_pos[i]
            
            y1 = self.y_pos[i] + k
            y2 = self.y_pos[i] - k
            
            if self.cat[i] == 'Fibre':
                if (np.absolute(self.redchi1[i]) <= np.absolute(self.redchi2[i])) or np.isnan(self.redchi2[i]):
                    angle = self.angle1[i]
                    transAngle = math.radians(180-angle)
                elif (np.absolute(self.redchi1[i]) > np.absolute(self.redchi2[i])) or np.isnan(self.redchi1[i]):
                    angle = self.angle2[i]
                    transAngle = math.radians(360-angle)
                    
                xnew1, ynew1 = self.rotate([self.x_pos[i], self.y_pos[i]], [x1, y1], transAngle)
                xnew2, ynew2 = self.rotate([self.x_pos[i], self.y_pos[i]], [x2, y2], transAngle)
                
                plt.plot([xnew1, xnew2],[ynew1, ynew2],'k')

            elif self.cat[i] == 'Resin':
                plt.plot(self.x_pos[i],self.y_pos[i],'.r', markersize=4)
            elif self.cat[i] == 'Off':
                # Uncomment to reveal empty space (off sample)
                #plt.plot(self.x_pos[i],self.y_pos[i],'.g', markersize=4)
                pass

            if label.lower() == 'on':
                plt.text(self.x_pos[i], self.y_pos[i], i, fontsize=8, ha='center', va='center')
            elif label.lower() == 'off':
                continue
            else:
                raise ValueError('Please specify ON or OFF label plot')


        
        
        
        
        






