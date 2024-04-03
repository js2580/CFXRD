# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 20:29:19 2023

@author: js2580
"""

import math
import numpy as np

def two_theta_to_d_spacing(two_theta, wavelength):
    """
    Convert 2θ to d-spacing for X-ray diffraction.

    Parameters:
        two_theta (float): The angle in degrees (2θ).
        wavelength (float): The X-ray wavelength in angstroms.

    Returns:
        float: The corresponding d-spacing.
    """
    two_theta_rad = math.radians(two_theta)
    d_spacing = wavelength / (2 * math.sin(two_theta_rad / 2))
    return d_spacing

def d_spacing_to_two_theta(d_spacing, wavelength):
    """
    Convert d-spacing to 2θ for X-ray diffraction.

    Parameters:
        d_spacing (float): The d-spacing in angstroms.
        wavelength (float): The X-ray wavelength in angstroms.

    Returns:
        float: The corresponding 2θ angle in degrees.
    """
    two_theta_rad = 2 * math.asin(wavelength / (2 * d_spacing))
    two_theta = math.degrees(two_theta_rad)
    return two_theta


def shift(disp=3, R=127.032, twoTheta=np.linspace(7, 43, 43-7+1)):
    # disp = disp #mm
    # R = R #mm
    # twoTheta = np.linspace(7, 43, 43-7+1)
    top = disp*np.sin(np.radians(twoTheta * 2))
    bot = 2*(R - disp * np.sin(np.radians(twoTheta))**2)
    eta = np.arctan(top/bot) * 180 / np.pi
    return eta

wavelength = 0.8156856 #Angstrom

print(d_spacing_to_two_theta(3.520, wavelength))
print(d_spacing_to_two_theta(2.08471, wavelength))

print(shift(disp=-1, R=127.03, twoTheta= np.array([13.306973304218534, 22.563673150600582])))


print(two_theta_to_d_spacing(13.306973304218534-(-0.10098579), wavelength))
print(two_theta_to_d_spacing(22.563673150600582-(-0.1596357), wavelength))