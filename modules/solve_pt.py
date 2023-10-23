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
from copy import deepcopy

from modules.compute_moist_adiabat import compute_moist_adiabat
from modules.dry_adiabat_timestep import compute_dry_adiabat

def RadConvEqm(dirs, time, atm, standalone:bool, cp_dry:bool, trppD:bool, calc_cf:bool, rscatter:bool, 
               pure_steam_adj=False, surf_dt=False, cp_surf=1e5, mix_coeff_atmos=1e6, mix_coeff_surf=1e6):
    """Sets the atmosphere to a temperature profile using the general adiabat. 
    
    Optionally does radiative time-stepping, but this is deprecated.

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


def MCPA(dirs, atm, standalone:bool, trppD:bool, rscatter:bool):
    """Calculates the temperature profile using the multiple-condensible pseudoadiabat.

    Prescribes a stratosphere, and also calculates fluxes.

    Parameters
    ----------
        dirs : dict
            Named directories
        atm : atmos
            Atmosphere object from atmosphere_column.py
        standalone : bool
            Running AEOLUS as standalone code?
        trppD : bool 
            Calculate tropopause dynamically?
        rscatter : bool
            Include rayleigh scattering?
            
    """

    ### Moist/general adiabat
    return compute_moist_adiabat(atm, dirs, standalone, trppD, False, rscatter)

def MCPA_CL(dirs, atm_input, standalone:bool, trppD:bool, rscatter:bool, atm_bc:int=0):
    """Calculates the temperature profile using the multiple-condensible pseudoadiabat and steps T_surf to conserve energy.

    Prescribes a stratosphere, and also calculates fluxes. Only works when used with PROTEUS

    Parameters
    ----------
        dirs : dict
            Named directories
        atm : atmos
            Atmosphere object from atmosphere_column.py
        standalone : bool
            Running AEOLUS as standalone code?
        trppD : bool 
            Calculate tropopause dynamically?
        rscatter : bool
            Include rayleigh scattering?
        
        atm_bc : int
            Where to measure boundary condition flux (0: TOA, 1: Surface).
    """

    # Iterate surface temperature
    niter = 10
    conv = False
    for i in range(niter):
        print("Iter surface %02d/%02d" % (i+1,niter))

        # Copy atm object 
        if 'atm' in globals():
            del atm
        atm = deepcopy(atm_input)

        # Handle conduction
        if i > 0:
            T_surf_old = atm.ts
            T_surf_new = atm.tmp_magma - F_atm_new * atm.skin_d / atm.skin_k
            T_surf_new = max(atm.minT, T_surf_new)
            
            if abs(T_surf_new-T_surf_old) < 5.0:
                weights = [0.0, 1.0]  
                conv = True
            else:
                weights = [1.0, 0.0]  
            atm.ts = T_surf_new * weights[0] + T_surf_old * weights[1]

        print("   T_surf = %g K"        % atm.ts)

        # Run computation
        atm.tmpl[-1] = atm.ts
        atm = compute_moist_adiabat(atm, dirs, standalone, trppD, False, rscatter)

        # Parse result
        if atm_bc == 0:
            F_atm_new = atm.net_flux[0]  
        else:
            F_atm_new = atm.net_flux[-1]  

        print("   F_atm  = %1.2e W m-2" % F_atm_new)

        # Save new surface temperature
        print(" ")

        if conv:
            print("Exiting surface loop early")
            break

    return atm

