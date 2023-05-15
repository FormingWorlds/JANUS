#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 13:05:00 2023

@authors: 
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
Harrison Nicholls (HN)
"""

import pickle as pkl

from modules.compute_moist_adiabat import compute_moist_adiabat
from modules.dry_adiabat_timestep import compute_dry_adiabat

def RadConvEqm(dirs, time, atm, standalone, cp_dry, trppD, calc_cf, rscatter, pure_steam_adj=False, surf_dt=False, cp_surf=1e5, mix_coeff_atmos=1e6, mix_coeff_surf=1e6):
    """Solves system for radiative-convective eqm

    Parameters
    ----------
        dirs : dict
            Named directories
        time : dict
            Dictionary of time values, including stellar age and evolution of planet
        atm : atmos
            Atmosphere object from atmosphere_column.py
        standalone : bool
            Running AEOLUS as standalone code?
        cp_dry : bool
            Compute dry adiabat case
        trppD : bool 
            Calculate tropopause dynamically?
        calc_cf : bool
            Calculate contribution function?
        pure_steam_adj : bool
            Use pure steam adjustment?
        surf_dt : float
            Timestep to use for T_surf timestepping cases
        cp_surf : float
            Surface heat capacity in T_surf timestepping cases
        mix_coeff_atmos : float
            Mixing coefficient (atmosphere) for T_surf timestepping cases?
        mix_coeff_surf : float
            Mixing coefficient (surface) for T_surf timestepping cases?
            
    """

    ### Moist/general adiabat
    atm_moist = compute_moist_adiabat(atm, dirs, standalone, trppD, calc_cf, rscatter)

    ### Dry adiabat
    if cp_dry == True:

        # Compute dry adiabat  w/ timestepping
        print("atm.toa_heating in RadConvEqm = ", atm.toa_heating)
        atm_dry   = compute_dry_adiabat(atm, dirs, standalone, calc_cf, rscatter, pure_steam_adj, surf_dt, cp_surf, mix_coeff_atmos, mix_coeff_surf)

        if standalone == True:
            print("Net, OLR => moist:", str(round(atm_moist.net_flux[0], 3)), str(round(atm_moist.LW_flux_up[0], 3)) + " W/m^2", end=" ")
            print("| dry:", str(round(atm_dry.net_flux[0], 3)), str(round(atm_dry.LW_flux_up[0], 3)) + " W/m^2", end=" ")
            print()
    else: 
        atm_dry = {}
    
    # Plot
    if standalone == True:
        #plot_flux_balance(atm_dry, atm_moist, cp_dry, time, dirs)
        # Save to disk
        with open(dirs["output"]+"/"+str(int(time["planet"]))+"_atm.pkl", "wb") as atm_file: 
            pkl.dump(atm_moist, atm_file, protocol=pkl.HIGHEST_PROTOCOL)

    return atm_dry, atm_moist
