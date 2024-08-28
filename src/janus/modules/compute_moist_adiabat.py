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
from janus.modules.find_tropopause import find_tropopause
from janus.modules.set_stratosphere import set_stratosphere
from janus.modules.water_cloud import simple_cloud
from janus.modules.relative_humidity import compute_Rh

import logging
log = logging.getLogger(__name__)

import janus.utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles
import janus.utils.socrates as socrates

def compute_moist_adiabat(atm, dirs, standalone, trppD, rscatter=False):
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
            
    """

    atm_moist = copy.deepcopy(atm)

    # Build general adiabat structure
    atm_moist = ga.general_adiabat(atm_moist)
    
    atm_moist.rh = compute_Rh(atm_moist)

    if atm.do_cloud:
        atm_moist = simple_cloud(atm_moist) # Before radiation, set up the cloud for Socrates using the current PT profile

    # Run SOCRATES
    atm_moist = socrates.radCompSoc(atm_moist, dirs, recalc=False, rscatter=rscatter)

    if standalone == True:
        log.info("w/o stratosphere (net, OLR): %.3f, %.3f W/m^2" % ( atm_moist.net_flux[0] , atm_moist.LW_flux_up[0]))

    # Calculate tropopause
    if (trppD == True) or (atm_moist.trppT > atm_moist.minT):
      
        # Find tropopause index
        atm_moist = find_tropopause(atm_moist,trppD, verbose=standalone)

        # Reset stratosphere temperature and abundance levels
        atm_moist = set_stratosphere(atm_moist)

        if atm.do_cloud:
            atm_moist = simple_cloud(atm_moist) # Update cloud location after previous PT changes
        # Recalculate fluxes w/ new atmosphere structure
        atm_moist = socrates.radCompSoc(atm_moist, dirs, recalc=True, rscatter=rscatter)

        if standalone == True:
            log.info("w/ stratosphere (net, OLR): %.3f, %.3f W/m^2" % (atm_moist.net_flux[0] , atm_moist.LW_flux_up[0] )) 

    return atm_moist
