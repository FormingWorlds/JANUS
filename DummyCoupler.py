"""
RadInteriorCoupler.py
"""

import numpy as np
import SocRadConv
import SocRadModel
import os, shutil
import coupler_utils
import spider_coupler_utils
import plot_atmosphere


# Define output output output directory
output_dir = os.getcwd()+"/output/"

stellar_toa_heating = 1000.0
p_s = 3.e5

time_current        = 0  # K
surfaceT_current    = 1000  # K
h2o_kg              = 0.0001 # k
co2_kg              = 0.00000 # kg
h2_kg               = 0.0  # kg
ch4_kg              = 0  # kg
co_kg               = 0  # kg
n2_kg               = 0  # kg
o2_kg               = 0  # kg
he_kg               = 0  # kg

h2o_ratio, co2_ratio, h2_ratio, ch4_ratio, co_ratio, n2_ratio, o2_ratio, he_ratio = coupler_utils.CalcMolRatios(h2o_kg, co2_kg, h2_kg, ch4_kg, co_kg, n2_kg, o2_kg, he_kg)

    # Calculate OLR flux for a given surface temperature w/ SOCRATES
heat_flux = str(SocRadConv.RadConvEqm(output_dir, time_current, surfaceT_current, stellar_toa_heating, p_s, h2o_ratio, co2_ratio, h2_ratio, ch4_ratio, co_ratio, n2_ratio, o2_ratio, he_ratio)) # W/m^2


