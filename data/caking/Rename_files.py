# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 16:20:01 2023

@author: js2580
"""
import os
import re

suffix = '_caking.nxs'

# Look for all files in the folder
filesList = [f for f in os.listdir() if os.path.isfile(f)]
if 'Rename_files.py' in filesList:
    filesList.remove('Rename_files.py') # Igonre the script itself

pattern = '_eigerScan\d{6}' #\d{6} is 6 digits 
match=(re.search(pattern, filesList[0]))

# Change filenames
initialIndex = match.span(0)[0]
lastIndex = match.span(0)[1]

for count, name in enumerate(filesList):
    initial = name[initialIndex:lastIndex]
    os.rename('./' + name, './' + initial + suffix)
