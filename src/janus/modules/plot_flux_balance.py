#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:14:29 2023

@authors: 
Tim Lichtenberg (TL)    
Ryan Boukrouche (RB)
"""

import matplotlib.pyplot as plt
import numpy as np
import math 
#import json

from janus.modules.spectral_planck_surface import surf_Planck_nu
import janus.utils.GeneralAdiabat as ga # Moist adiabat with multiple condensibles

# Plotting
def plot_fluxes(atm,filename='output/fluxes.pdf'):

    fig,ax = plt.subplots()

    ax.axvline(0,color='black',lw=0.8)

    arr_p = atm.pl*1e-5

    ax.plot(atm.flux_up_total,arr_p,color='red',label='UP',lw=1)
    ax.plot(atm.SW_flux_up   ,arr_p,color='red',label='UP SW',linestyle='dotted',lw=2)
    ax.plot(atm.LW_flux_up   ,arr_p,color='red',label='UP LW',linestyle='dashed',lw=1)

    ax.plot(-1.0*atm.flux_down_total,arr_p,color='green',label='DN',lw=2)
    ax.plot(-1.0*atm.SW_flux_down   ,arr_p,color='green',label='DN SW',linestyle='dotted',lw=3)
    ax.plot(-1.0*atm.LW_flux_down   ,arr_p,color='green',label='DN LW',linestyle='dashed',lw=2)

    ax.plot(atm.net_flux ,arr_p,color='black',label='NET')

    ax.set_xlabel("Upward-directed flux [W m-2]")
    ax.set_xscale("symlog")
    ax.set_ylabel("Pressure [Pa]")
    ax.set_yscale("log")
    ax.set_ylim(top=np.amin(arr_p)/1.5, bottom=np.amax(arr_p)*1.5)

    ax.legend(loc='upper left')

    fig.savefig(filename, bbox_inches='tight', dpi=190)
    plt.close()

