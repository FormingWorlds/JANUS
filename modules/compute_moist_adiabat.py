#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:01:27 2023

@authors: 
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
"""

import copy
import numpy as np
from modules.find_tropopause import find_tropopause
from modules.set_stratosphere import set_stratosphere
from modules.water_cloud import simple_cloud
from modules.relative_humidity import compute_Rh

import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
import utils.socrates as socrates

def compute_moist_adiabat(atm, dirs, standalone, trppD, rscatter=False, do_cloud=False):
    """Compute moist adiabat case

    Parameters
    ----------
        atm : atmos
            Atmosphere object from atmosphere_column.py
        dirs : dict
            Named directories
        standalone : bool
            Running JANUS as standalone code?
        trppD : bool 
            Calculate tropopause dynamically?
        rscatter : bool
            Include Rayleigh scattering?
        do_cloud : bool
            Include water cloud radiation?
            
    """

    atm_moist = copy.deepcopy(atm)

    # Build general adiabat structure
    atm_moist = ga.general_adiabat(atm_moist)
    
    atm_moist.rh = compute_Rh(atm_moist)

    if do_cloud:
        atm_moist = simple_cloud(atm_moist) # Before radiation, set up the cloud for Socrates using the current PT profile

    # Run SOCRATES
    atm_moist = socrates.radCompSoc(atm_moist, dirs, recalc=False, rscatter=rscatter, do_cloud=do_cloud)

    if standalone == True:
        print("w/o stratosphere (net, OLR): " + str(round(atm_moist.net_flux[0], 3)) +" , "+str(round(atm_moist.LW_flux_up[0], 3)) + " W/m^2")

    # Calculate tropopause
    if (trppD == True) or (atm_moist.trppT > atm_moist.minT):
      
        # Find tropopause index
        atm_moist = find_tropopause(atm_moist,trppD, verbose=standalone)

        # Reset stratosphere temperature and abundance levels
        atm_moist = set_stratosphere(atm_moist)

        if do_cloud:
            atm_moist = simple_cloud(atm_moist) # Update cloud location after previous PT changes
        # Recalculate fluxes w/ new atmosphere structure
        atm_moist = socrates.radCompSoc(atm_moist, dirs, recalc=True, rscatter=rscatter, do_cloud=do_cloud)

        if standalone == True:
            print("w/ stratosphere (net, OLR): " + str(round(atm_moist.net_flux[0], 3)) + " , " + str(round(atm_moist.LW_flux_up[0], 3)) + " W/m^2")

    return atm_moist