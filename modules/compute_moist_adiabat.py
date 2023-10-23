#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:01:27 2023

@authors: 
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
"""

import copy
from modules.find_tropopause import find_tropopause
from modules.set_stratosphere import set_stratosphere

import utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
import utils.socrates as socrates

def compute_moist_adiabat(atm, dirs, standalone, trppD, calc_cf=False, rscatter=False):
    """Compute moist adiabat case

    Parameters
    ----------
        atm : atmos
            Atmosphere object from atmosphere_column.py
        dirs : dict
            Named directories
        standalone : bool
            Running AEOLUS as standalone code?
        trppD : bool 
            Calculate tropopause dynamically?
        calc_cf : bool
            Calculate contribution function?
        rscatter : bool
            Include Rayleigh scattering?
            
    """

    atm_moist = copy.deepcopy(atm)

    # Build general adiabat structure
    atm_moist = ga.general_adiabat(atm_moist)

    # Run SOCRATES
    atm_moist = socrates.radCompSoc(atm_moist, dirs, recalc=False, calc_cf=calc_cf, rscatter=rscatter)

    if standalone == True:
        print("w/o stratosphere (net, OLR):", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

    # Calculate tropopause dynamically
    if (trppD == True) or (atm_moist.trppT > 0.0) or (atm_moist.minT > 0.0):
      
        # Find tropopause index
        atm_moist = find_tropopause(atm_moist,trppD, verbose=standalone)

        # Reset stratosphere temperature and abundance levels
        atm_moist = set_stratosphere(atm_moist)

        # Recalculate fluxes w/ new atmosphere structure
        atm_moist = socrates.radCompSoc(atm_moist, dirs, recalc=True, calc_cf=calc_cf, rscatter=rscatter)

        if standalone == True:
            print("w/ stratosphere (net, OLR):", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)), "W/m^2")

    return atm_moist
