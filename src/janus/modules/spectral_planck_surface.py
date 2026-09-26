#!/usr/bin/env python3
"""
Created on Mon Jan 23 12:16:20 2023

@authors:
Tim Lichtenberg (TL)
Ryan Boukrouche (RB)
"""

import numpy as np


def surf_Planck_nu(atm):
    B = np.zeros(len(atm.band_centres))
    c1 = 1.191042e-5
    c2 = 1.4387752
    for i in range(len(atm.band_centres)):
        nu = atm.band_centres[i]
        B[i] = c1 * nu**3 / (np.exp(c2 * nu / atm.ts) - 1)
    B = (1.0 - atm.albedo_s) * np.pi * B * atm.band_widths / 1000.0
    return B
