#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 17:18:48 2024

@authors:    
Ryan Boukrouche (RB)
"""

import numpy as np
from janus.utils.find_lcl import find_intersection
import janus.utils.GeneralAdiabat as ga
    
def compute_Rh(atm):
    """Compute the relative humidity.

    Parameters
    ----------
        atm : atmos
            Atmosphere object

    """

    #partial_p = atm.p_vol['H2O']-atm.x_cond['H2O']*atm.p # Remove the condensates so that only the vapor is accounted for.
    partial_p = atm.x_gas['H2O']*atm.p
    p_sat = [ga.p_sat('H2O', atm.tmp[i]) for i in range(len(atm.p))]
    R_h = partial_p/p_sat
    return R_h

def update_for_constant_RH(atm, R_h):
    """Compute the partial pressure and mixing ratio of the condensible that keeps the relative humidity constant below the LCL. Also update the total pressure.

    Parameters
    ----------
        atm : atmos
            Atmosphere object

        R_h : float
            Relative humidity

    """

    p_sat = [ga.p_sat('H2O', atm.tmp[i]) for i in range(len(atm.p))]

    i_lcl, lcl = find_intersection(p_sat, atm.p_vol['H2O'], tolerance=1e0)

    atm.p_vol['H2O'] = R_h * p_sat
    atm.p = sum(atm.p_vol[molecule] for molecule in atm.vol_list.keys()) # Updating the total pressure
    atm.ps = float(sum(atm.p_vol[molecule][-1] for molecule in atm.vol_list.keys())) # Not forgetting the surface pressure

    return atm
