# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 09:39:42 2022

@author: js2580
"""

import os
import sys

filesList = os.listdir('./.')
filesList.remove('Rename_files.py')

# last = 'azimuthal_integration'
initialChar =int(sys.argv[1])
last = str(sys.argv[2])

for count, name in enumerate(filesList):
    initial = name[:initialChar]
    os.rename('./' + name, './' + initial + last)
    
    


