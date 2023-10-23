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
import scipy.optimize as optimise

from modules.compute_moist_adiabat import compute_moist_adiabat
from modules.dry_adiabat_timestep import compute_dry_adiabat
from utils.atmosphere_column import atmos

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

def MCPA_CL(dirs, atm_input, trppD:bool, rscatter:bool, atm_bc:int=0):
    """Calculates the temperature profile using the multiple-condensible pseudoadiabat and steps T_surf to conserve energy.

    Prescribes a stratosphere, and also calculates fluxes. Only works when used with PROTEUS

    Parameters
    ----------
        dirs : dict
            Named directories
        atm : atmos
            Atmosphere object from atmosphere_column.py
        trppD : bool 
            Calculate tropopause dynamically?
        rscatter : bool
            Include rayleigh scattering?
        
        atm_bc : int
            Where to measure boundary condition flux (0: TOA, 1: Surface).
    """
    
    # Store constants
    alpha =         atm_input.alpha_cloud
    toa_heating =   atm_input.toa_heating
    minT =          atm_input.minT
    nlev_save =     atm_input.nlev_save
    vol_list =      atm_input.vol_list
    pl_m =          atm_input.planet_mass
    pl_r =          atm_input.planet_radius
    ptop =          atm_input.ptop
    psurf =         atm_input.ps
    trppT =         atm_input.trppT
    skin_k =        atm_input.skin_k
    skin_d =        atm_input.skin_d
    tmp_magma =     atm_input.tmp_magma

    # Calculate conductive flux for a given atmos object 'a'
    def skin(a):
        return a.skin_k / a.skin_d * (a.tmp_magma - a.ts)
    
    # Initialise a new atmos object
    def ini_atm(Ts):
        atm = atmos(Ts, psurf, ptop, pl_r, pl_m , vol_mixing=vol_list, trppT=trppT, minT=minT, req_levels=nlev_save)
        atm.toa_heating = toa_heating
        atm.alpha_cloud = alpha
        atm.skin_k = skin_k
        atm.skin_d = skin_d 
        atm.tmp_magma = tmp_magma
        return atm
    
    # We want to optimise this function (returns residual of F_atm and F_skn, given T_surf)
    def func(x):

        print("Evaluated at T_surf = %.1f K" % x)
        atm = compute_moist_adiabat(ini_atm(x), dirs, False, trppD, False, rscatter)

        if atm_bc == 0:
            F_atm = atm.net_flux[0]  
        else:
            F_atm = atm.net_flux[-1]  

        return float(skin(atm) - F_atm)
    
    # Find root (T_surf) satisfying energy balance
    print("Solving for global energy balance with conductive lid")
    x0 = atm_input.ts
    x1 = atm_input.tmp_magma * 0.85
    sol,r = optimise.newton(func, x0, x1=x1, tol=1.0, maxiter=20, full_output = True)

    # Extract solution
    T_surf = float(sol)
    succ  = r.converged

    if not succ:
        print("WARNING: Did not find solution for surface skin balance")
    else:
        print("Found surface solution")
        
    # Get atmosphere state from solution value
    atm = compute_moist_adiabat(ini_atm(T_surf), dirs, False, trppD, False, rscatter)

    if atm_bc == 0:
        F_atm = atm.net_flux[0]  
    else:
        F_atm = atm.net_flux[-1]  

    F_olr = atm.LW_flux_up[0]

    print("    T_surf = %g K"       % T_surf)
    print("    F_atm  = %.4e W m-2" % F_atm)
    print("    F_skn  = %.4e W m-2" % skin(atm))
    print("    F_olr  = %.4e W m-2" % F_olr)

    return atm

