# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 21:08:35 2023

@author: js2580
"""
# ** Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.signal import find_peaks
import h5py 

# ** Generate Dummy motors positions

# ** Define parent directory path
ParentDir = '.\\'

motor_dir = ParentDir + '\\data\\motor\\'
motorX = 'motorX'
motorY = 'motorY'

header = f'{motorY} {motorX}'


setScan = [413129]
scanlist = setScan

# ** Define Number of rows and columns, starting positions and increment size
RowNo = 10
ColumnNO = 5

ystart = 0
xstart = 0

yincrement = 1
xincrement = 1

# ! Loop through each scanning set
for i, scanNo in enumerate(scanlist):

    print(scanNo)
    pos = []

    for n in range(0, RowNo): #iterate through y-axis
        ypos = ystart + yincrement * n
        for m in range(0, ColumnNO): #iterate through x-axis
            xpos = xstart + xincrement * m
            pos.append([ypos,xpos])

    pos = np.asarray(pos)            
    np.savetxt(f'{motor_dir}\\{scanNo}.dat', pos, fmt="%.2f", header=header, comments="")







