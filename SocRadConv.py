#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:20:27 2023

@authors: 
Mark Hammond (MH)
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
Harrison Nicholls (HN)

JANUS radiative-convective model, using SOCRATES for radiative-transfer.
"""

import matplotlib as mpl
mpl.use('Agg')

import time as t
import os, shutil
import numpy as np

from modules.stellar_luminosity import InterpolateStellarLuminosity
from modules.solve_pt import RadConvEqm
from modules.solve_pt import *
from modules.plot_flux_balance import plot_fluxes
from modules.plot_emission_spectrum import plot_emission
from utils.socrates import CleanOutputDir

import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
from utils.atmosphere_column import atmos
import utils.StellarSpectrum as StellarSpectrum
from utils.ReadSpectralFile import ReadBandEdges

####################################
##### Stand-alone initial conditions
####################################
if __name__ == "__main__":

    print("Start JANUS")

    # Set up dirs
    if os.environ.get('JANUS_DIR') == None:
        raise Exception("Environment variables not set! Have you sourced JANUS.env?")
    dirs = {
            "janus": os.getenv('JANUS_DIR')+"/",
            "output": os.getenv('JANUS_DIR')+"/output/"
            }

    start = t.time()
    ##### Settings

    # Planet 
    time = { "planet": 0., "star": 4e+9 } # yr,
    star_mass     = 1.0                 # M_sun, mass of star
    mean_distance = 1.0                 # au, orbital distance
    pl_radius     = 6.371e6             # m, planet radius
    pl_mass       = 5.972e24            # kg, planet mass

    # Boundary conditions for pressure & temperature
    T_surf        = 2800.0                # K
    P_top         = 1.0                  # Pa

    # Define volatiles by mole fractions
    # P_surf       =  300.0 * 1e5
    # vol_partial = {}
    # vol_mixing = { 
    #                 "CO2"  : 0.0,
    #                 "H2O"  : 1.0,
    #                 "N2"   : 0.0,
                # }
    
    # OR:
    # Define volatiles by partial pressures
    P_surf = 0.0
    vol_mixing = {}
    vol_partial = {
        "H2O" : 1.0e5,
        # "NH3" : 0.,
        "CO2" : 0.0,#2.0e5,
        "CH4" :  3.0e5,
        "CO" : 1.85989e5,
        # "O2" : 0.20e5,
        "N2" : 4.0e5,
        "H2" : 1.0e5,
        }

    # Stellar heating on/off
    stellar_heating = True

    # Rayleigh scattering on/off
    rscatter = False

    # Compute contribution function
    calc_cf = False

    # Pure steam convective adjustment
    pure_steam_adj = False

    # Tropopause calculation
    trppD = False   # Calculate dynamically?
    trppT = 10.0     # Fixed tropopause value if not calculated dynamically

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

    # Tidy directory
    if os.path.exists(dirs["output"]):
        shutil.rmtree(dirs["output"])
    os.mkdir(dirs["output"])

    # Move/prepare spectral file
    print("Inserting stellar spectrum")
    StellarSpectrum.InsertStellarSpectrum(
        # dirs["janus"]+"/spectral_files/shared/old/Dayspring256_old/Dayspring.sf",
        dirs["janus"]+"/spectral_files/Mallard/Mallard.sf",
        dirs["janus"]+"/spectral_files/stellar_spectra/Sun_t4_4Ga_claire_12.txt",
        dirs["output"]
    )

    band_edges = ReadBandEdges(dirs["output"]+"star.sf")

    # Create atmosphere object
    atm = atmos(T_surf, P_surf, P_top, pl_radius, pl_mass, band_edges,
                vol_mixing=vol_mixing, vol_partial=vol_partial, calc_cf=calc_cf, trppT=trppT, req_levels=100, water_lookup=water_lookup)

    # Set stellar heating on or off
    if stellar_heating == False: 
        atm.instellation = 0.
    else:
        atm.instellation = InterpolateStellarLuminosity(star_mass, time, mean_distance)
        print("Instellation:", round(atm.instellation), "W/m^2")

    # Set up atmosphere with general adiabat
    atm_dry, atm = RadConvEqm(dirs, time, atm, standalone=True, cp_dry=cp_dry, trppD=trppD, calc_cf=calc_cf, rscatter=rscatter, pure_steam_adj=pure_steam_adj, surf_dt=surf_dt, cp_surf=cp_surf, mix_coeff_atmos=mix_coeff_atmos, mix_coeff_surf=mix_coeff_surf) 

    # Plot abundances w/ TP structure
    if (cp_dry):
        ga.plot_adiabats(atm_dry,filename=dirs["output"]+"dry_ga.png")
        atm_dry.write_PT(filename=dirs["output"]+"dry_pt.tsv")
        plot_fluxes(atm_dry,filename=dirs["output"]+"dry_fluxes.png")

    ga.plot_adiabats(atm,filename= dirs["output"]+"moist_ga.png")
    atm.write_PT(filename= dirs["output"]+"moist_pt.tsv")
    atm.write_ncdf( dirs["output"]+"moist_atm.nc")
    plot_fluxes(atm,filename= dirs["output"]+"moist_fluxes.png")
    plot_emission(atm, dirs["output"]+"toa_emission.png", planck_surface=True, show_bands=True)

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])

    end = t.time()
    print("Runtime:", round(end - start,2), "s")

