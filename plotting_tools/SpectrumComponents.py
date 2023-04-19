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
import matplotlib.pyplot as plt

# Define output output output directory
output_dir = os.getcwd()+"/output/"

stellar_toa_heating = 1000.0
p_s                 = 3.e5
mr_list             = []
names_list          = []
all_spectra         = []


mr_list.append([0.0001,0,0,0,0])
names_list.append('H2O')

mr_list.append([0.,0.001,0,0,0])
names_list.append('CO2')

mr_list.append([0,0,0,0.001,0])
names_list.append('CO')

mr_list.append([0,0,0.003,0.0,0])
names_list.append('H2')

mr_list.append([0,0,0,0,0.01])
names_list.append('N2')

mr_list.append([0.0001,0.001,0.003,0.001,0.01])
names_list.append('All')

for i in range(len(mr_list)):
    mr = mr_list[i]
    time_current        = 0     # K
    surfaceT_current    = 300   # K
    h2o_kg              = mr[0] # k
    co2_kg              = mr[1] # kg
    h2_kg               = mr[2] # kg
    ch4_kg              = 0     # kg
    co_kg               = mr[3] # kg
    n2_kg               = mr[4] # kg
    o2_kg               = 1     # kg
    he_kg               = 0     # kg

    h2o_ratio, co2_ratio, h2_ratio, ch4_ratio, co_ratio, n2_ratio, o2_ratio, he_ratio = coupler_CalcMolRatios(h2o_kg, co2_kg, h2_kg, ch4_kg, co_kg, n2_kg, o2_kg, he_kg)

    # Calculate OLR flux for a given surface temperature w/ SOCRATES
    heat_flux, bands, spectrum = SocRadConv.RadConvEqm(output_dir, time_current, surfaceT_current, stellar_toa_heating, p_s, h2o_ratio, co2_ratio, h2_ratio, ch4_ratio, co_ratio, n2_ratio, o2_ratio, he_ratio)
    # heat_flux, bands, spectrum = RadConvEqm(output_dir, time_current, runtime_helpfile, stellar_toa_heating, atm_chemistry, loop_counter, SPIDER_options)  # W/m^2

    all_spectra.append(spectrum)

plt.figure()

for i in range(len(mr_list)-1):
    plt.plot(bands,all_spectra[i],'--',label=names_list[i])

plt.plot(bands,all_spectra[-1],'-',label=names_list[-1],color='k',lw=1.5)

plt.legend()
plt.xlim([0,2500])
plt.savefig('all_spectra.pdf', bbox_inches = 'tight')
