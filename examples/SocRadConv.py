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
import os, shutil, toml
import numpy as np
from importlib.resources import files

from janus.modules.stellar_luminosity import InterpolateStellarLuminosity
from janus.modules.solve_pt import RadConvEqm
from janus.modules.solve_pt import *
from janus.modules.plot_flux_balance import plot_fluxes
from janus.modules.plot_emission_spectrum import plot_emission
from janus.utils.socrates import CleanOutputDir

import janus.utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
from janus.utils.atmosphere_column import atmos
import janus.utils.StellarSpectrum as StellarSpectrum
from janus.utils.ReadSpectralFile import ReadBandEdges

####################################
##### Stand-alone initial conditions
####################################
if __name__ == "__main__":

    print("Start JANUS")

    # Set up dirs
    if os.environ.get('RAD_DIR') == None:
        raise Exception("Socrates environment variables not set! Have you installed Socrates and sourced set_rad_env?")
    dirs = {
            "janus": str(files("janus"))+"/",
            "output": os.path.abspath(os.getcwd())+"/output/"
            }

    start = t.time()
    ##### Settings
    cfg_file =  dirs["janus"]+"data/tests/config_janus.toml"
    with open(cfg_file, 'r'):
          cfg = toml.load(cfg_file)  

    # Planet
    time = { "planet": cfg['planet']['time'], "star": cfg['star']['time']}
    star_mass = cfg['star']['star_mass']
    mean_distance = cfg['star']['mean_distance']
 
    # Define volatiles by partial pressures
    vol_mixing = {}
    vol_partial = {
        "H2O" : 9.0e5,
        "CO2" : 2.0e5,
        "CH4" : 3.0e5,
        "CO" :  5.0e5,
        "N2" :  4.0e5,
        "H2" :  1.0e5,
        }

    # Tidy directory
    if os.path.exists(dirs["output"]):
        shutil.rmtree(dirs["output"])
    os.mkdir(dirs["output"])

    # Move/prepare spectral file
    print("Inserting stellar spectrum")
    StellarSpectrum.InsertStellarSpectrum(
        dirs["janus"]+"data/spectral_files/Dayspring/256/Dayspring.sf",
        dirs["janus"]+"data/spectral_files/stellar_spectra/Sun_t4_4Ga_claire_12.txt",
        dirs["output"]
    )

    band_edges = ReadBandEdges(dirs["output"]+"star.sf")

    # Create atmosphere object
    atm = atmos.from_file(cfg_file, band_edges, vol_mixing=vol_mixing, vol_partial=vol_partial)

    # Set stellar heating on or off
    if cfg['star']['stellar_heating'] == False: 
        atm.instellation = 0.
    else:
        atm.instellation = InterpolateStellarLuminosity(star_mass, time, mean_distance)
        print("Instellation:", round(atm.instellation), "W/m^2")

    # Set up atmosphere with general adiabat
    atm_dry, atm = RadConvEqm(dirs,
                              time,
                              atm,
                              standalone=True,
                              cp_dry=False,
                              trppD=False, # Calculate dynamically?
                              rscatter=False, # Rayleigh scattering on/off
                              pure_steam_adj=False, # Pure steam convective adjustment
                              surf_dt=False, # Surface temperature time-stepping
                              # Options activated by surf_dt
                              cp_surf=1e5, # Heat capacity of the ground [J.kg^-1.K^-1],
                              mix_coeff_atmos=1e6, # mixing coefficient of the atmosphere [s]
                              mix_coeff_surf=1e6 # mixing coefficient at the surface [s]
                              )

    # Plot abundances w/ TP structure
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

