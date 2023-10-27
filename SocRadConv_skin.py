#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:20:27 2023

@authors: 
Mark Hammond (MH)
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
Harrison Nicholls (HN)

AEOLUS radiative-convective model, using SOCRATES for radiative-transfer.
"""

import matplotlib as mpl
mpl.use('Agg')

import time as t
import os, shutil
import numpy as np

from modules.stellar_luminosity import InterpolateStellarLuminosity
from AEOLUS.modules.solve_pt import *
from modules.plot_flux_balance import plot_fluxes
from utils.socrates import CleanOutputDir

import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
from utils.atmosphere_column import atmos
import utils.StellarSpectrum as StellarSpectrum

####################################
##### Stand-alone initial conditions
####################################
if __name__ == "__main__":

    print("Start AEOLUS")

    start = t.time()
    ##### Settings

    # Planet 
    time = { "planet": 0., "star": 4e+9 } # yr,
    star_mass     = 1.0                 # M_sun, mass of star
    mean_distance = 1.0                 # au, orbital distance
    pl_radius     = 7.12e6             # m, planet radius
    pl_mass       = 8.21e24           # kg, planet mass

    # Boundary conditions for pressure & temperature
    T_magma       = 2579.6                # K
    P_top         = 0.1                  # Pa

    # Define volatiles by partial pressures
    P_surf = 0.0
    vol_mixing = {}
    vol_partial = {
        "H2O" : 0.92529e5,
        "NH3" : 0.0,
        "CO2" : 5.98731e5,
        "CH4" : 0.0,
        "CO" : 115.89698e5,
        "O2" : 0.0,
        "N2" : 1.77751e5,
        "H2" : 2.38831e5
        }

    # Rayleigh scattering on/off
    rscatter = False

    # Compute contribution function
    calc_cf = False

    # Pure steam convective adjustment
    pure_steam_adj = False

    # Tropopause calculation
    trppD = False   # Calculate dynamically?
    trppT = 50.0     # Fixed tropopause value if not calculated dynamically

    # Water lookup tables enabled (e.g. for L vs T dependence)
    water_lookup = False
    
    # Surface temperature time-stepping
    surf_dt = False
    cp_dry = False
    # Options activated by surf_dt
    cp_surf = 1e5         # Heat capacity of the ground [J.kg^-1.K^-1]
    mix_coeff_atmos = 1e6 # mixing coefficient of the atmosphere [s]
    mix_coeff_surf  = 1e6 # mixing coefficient at the surface [s]

    # Instellation scaling | 1.0 == no scaling
    Sfrac = 1.0

    ##### Function calls

    # Set up dirs
    dirs = {
            "aeolus": os.getenv('AEOLUS_DIR')+"/",
            "output": os.getenv('AEOLUS_DIR')+"/output/"
            }
    
    # Tidy directory
    if os.path.exists(dirs["output"]):
        shutil.rmtree(dirs["output"])
    os.mkdir(dirs["output"])

    # Create atmosphere object
    atm = atmos(T_magma, P_surf, P_top, pl_radius, pl_mass, 
                vol_partial=vol_partial, trppT=trppT)
    atm.tmp_magma = T_magma

    # Set stellar heating on or off
    atm.toa_heating = 4.349e+04

    # Move/prepare spectral file
    print("Inserting stellar spectrum")

    StellarSpectrum.InsertStellarSpectrum(
        dirs["aeolus"]+"/spectral_files/Mallard/Mallard",
        dirs["aeolus"]+"/spectral_files/stellar_spectra/Sun_t4_4Ga_claire_12.txt",
        dirs["output"]+"runtime_spectral_file"
    )
    
    atm = MCPA_CL(dirs, atm, trppD, rscatter, T_surf_max=T_magma)

    ga.plot_adiabats(atm,filename="output/skin_ga.pdf")
    atm.write_PT(filename="output/skin_pt.tsv")
    atm.write_ncdf("output/skin_atm.nc")
    plot_fluxes(atm,filename="output/skin_fluxes.pdf")

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])

    end = t.time()
    print("Runtime:", round(end - start,2), "s")

