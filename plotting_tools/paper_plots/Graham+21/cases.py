# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 16:50:56 2020

@authors: RJ, TL
Plotting script for cases 1–3
"""

import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir('C:/Users/grahamr/Documents/GitHub/soc-rad-conv/')
import GeneralAdiabat as ga
os.chdir('C:/Users/grahamr/Documents/GitHub/soc-rad-conv/plotting_tools/paper_plots/Graham+21')

from matplotlib import cm
import seaborn as sns
from atmosphere_column import atmos

# Line settings
ls_moist    = 2.5
ls_dry      = 2.0
ls_ind      = 1.5

#### LOOP OVER CASES
for set_idx, setting in enumerate([ "case1", "case2", "case3" ]): # "case1", "case2", "case3"
    #for set_idx, setting in enumerate([ "case1"]):#, "case2", "case3" ]): # "case1", "case2", "case3"

    ### Initial conditions in three settings

    # Magma ocean case
    if setting == "case1":
        name = "Magma Ocean"
        time = { "planet": 0., "star": 100e+06 } # yr,
        Mstar           = 1.0 
        mean_distance   = 1.0
        T_surf          = 1500

        # When P_surf is given as str or == 0., then vol_list defined in Pa
        P_surf          = "calc"

        vol_list        = { 
                            "H2O" :  500e+5,
                            "NH3" :  0.,
                            "CO2" :  100e+5,
                            "CH4" :  0.,
                            "CO"  :  0.,
                            "O2"  :  0.,
                            "N2"  :  1e+5,
                            "H2"  :  100e+5,
                          }

        case_color = "CO2"

    # LATE VENEER case
    if setting == "case2":
        name = "Late Veneer"
        time = { "planet": 0., "star": 500e+06 } # yr,
        Mstar           = 1.0 
        mean_distance   = 1.0
        
        # When P_surf is given as str or == 0., then vol_list defined in Pa
        P_surf          = "calc"

        # Plot color
        case_color = "H2O"

        ############## Tsurf = 700 K
        T_surf          = 700
        # # Maximum late veneer, 2 bar CO2
        # vol_list        = { 
        #                   "H2O" :  500e+5,
        #                   "NH3" :  0.07e+5,
        #                   "CO2" :  1.3e+5,
        #                   "CH4" :  0.1e+5,
        #                   "CO"  :  50e+5,
        #                   "O2"  :  0.e+5,
        #                   "N2"  :  0.07e+5,
        #                   "H2"  :  50e+5,
        #                 }
        # # Maximum late veneer, 50 bar CO2
        # vol_list        = { 
        #                   "H2O" :  500e+5,
        #                   "NH3" :  0.08e+5,
        #                   "CO2" :  0e+5,
        #                   "CH4" :  4e+5,
        #                   "CO"  :  50e+5,
        #                   "O2"  :  0e+5,
        #                   "N2"  :  0.09e+5,
        #                   "H2"  :  50e+5,
        #                 }
        # # Vesta, 2 bar CO2
        # vol_list        = { 
        #                   "H2O" :  500e+5,
        #                   "NH3" :  0.0013e+5,
        #                   "CO2" :  0.3e+5,
        #                   "CH4" :  0.05e+5,
        #                   "CO"  :  0.0008e+5,
        #                   "O2"  :  0e+5,
        #                   "N2"  :  0.3e+5,
        #                   "H2"  :  2e+5,
        #                 }
        # Vesta, 50 bar CO2
        vol_list        = { 
                          "H2O" :  500e+5,
                          "NH3" :  0.007e+5,
                          "CO2" :  48e+5,
                          "CH4" :  1e+5,
                          "CO"  :  0.05e+5,
                          "O2"  :  0e+5,
                          "N2"  :  1e+5,
                          "H2"  :  8e+5,
                        }

    # Hadean Earth / OHZ case
    if setting == "case3":

        name = r"Exo-Earth"

        # Planet age and orbit
        time = { "planet": 0., "star": 4567e+06 } # yr,

        # Star mass, M_sun
        Mstar           = 1.0 

        # Planet-star distance, au
        mean_distance   = 1.5

        # Surface pressure range (Pa)
        # When P_surf is given as str or == 0., then vol_list defined in Pa
        P_surf          = "calc"

        # Surface temperature range (K)
        T_surf          = 290

        # Volatiles considered
        vol_list    = { 
                      "H2O" :  0.01e+5,
                      "NH3" :  0.,
                      "CO2" :  29e+5,
                      "CH4" :  0e+5,
                      "CO"  :  0.,
                      "O2"  :  0.,
                      "N2"  :  1e+5,
                      "H2"  :  0e+5,
                    }

        # Plot color
        case_color = "H2"

    alpha_list = [ 0.0, 0.1, 1.0 ] 
    fig, axes = plt.subplots(len(alpha_list), 4, figsize=(13,7),sharex='col',sharey='all')
    plt.subplots_adjust(top = 0.99, bottom=0.01, hspace=0.05, wspace=0.05)

    if np.ndim(axes) == 1:
        axes = np.expand_dims(axes,axis=0)

    for n, alpha in enumerate(alpha_list):

        print("---------", setting, "case:", name, "alpha:", alpha)
        
        atm                = atmos(T_surf, P_surf, vol_list)
        P_surf = atm.ps
        
        # Calculate moist adiabat + condensation
        atm.alpha_cloud         = alpha
        atm                     = ga.general_adiabat(atm)
        
        T_sat_array    = np.linspace(20,3000,1000) 
        p_partial_sum  = np.zeros(len(atm.tmp))
        
        vol_list_sorted = {k: v for k, v in sorted(atm.vol_list.items(), key=lambda item: item[1])}
        
        # Individual species
        for vol in vol_list_sorted.keys():
        
            # Only if volatile is present
            if atm.vol_list[vol] > 1e-10:
        
                # Saturation vapor pressure for given temperature
                #Psat_array = [ ga.p_sat(vol, T)/1e+5 for T in T_sat_array ]
                #axes[n,0].semilogy( T_sat_array, Psat_array, lw=ls_ind, ls=":", color=ga.vol_colors[vol][4]) # , label=r'$p_\mathrm{sat}$'+ga.vol_latex[vol]
        
                # Plot partial pressures
                #axes[n,0].semilogy(atm.tmp, atm.p_vol[vol]/1e+5, color=ga.vol_colors[vol][4], lw=ls_ind, ls="-",alpha=0.99) # , label=r'$p$'+ga.vol_latex[vol]
                
                # Plot dew-point temperatures as functions of pressure, given the partial pressure for a given species at a given pressure level
                axes[n,0].semilogy(np.vectorize(ga.Tdew)(vol,atm.pl_vol[vol]),atm.pl/1e+5, color = ga.vol_colors[vol][4], lw = ls_ind, ls='-',alpha=0.99)
                # Sum up partial pressures
                p_partial_sum += atm.p_vol[vol]
        
                # Plot individual molar concentrations
                axes[n,1].semilogy(atm.x_cond[vol],atm.p/1e+5, color=ga.vol_colors[vol][4], lw=ls_ind, ls="--") # , label=ga.vol_latex[vol]+" cond."
                axes[n,1].semilogy(atm.x_gas[vol],atm.p/1e+5, color=ga.vol_colors[vol][4], lw=ls_ind, ls="-") # , label=ga.vol_latex[vol]+" gas" # , label=ga.vol_latex[vol]
                
        # General moist adiabat
        axes[n,0].semilogy(atm.tmp, atm.p/1e+5, color=ga.vol_colors["black_1"], lw=ls_moist,label="Adiabat",alpha=0.99)
        
        # Phase molar concentrations
        axes[n,1].semilogy(atm.xd+atm.xv,atm.p/1e+5, color=ga.vol_colors["black_2"], lw=ls_ind, ls=":") # , label=r"Gas phase"
        
        # Molar masses
        axes[n,2].semilogy(atm.mu*1000, atm.p/1e+5, color=ga.vol_colors["black_2"],lw=ls_ind)
        
        # Specific heats
        axes[n,3].semilogy(atm.cp, atm.p/1e+5, color=ga.vol_colors["black_2"],lw=ls_ind)
        
        fs_l = 16
        fs_m = 14
        fs_s = 12
        
        axes[n,0].invert_yaxis()
       
        if n == 0:
            pass
        axes[n,0].set_xlim([20,np.max(atm.ts)])
        

        axes[n,0].set_ylim(top=atm.ptop/1e+5,bottom=atm.ps/1e+5)
        
        axes[n,2].set_yscale('log')
        
        axes[n,1].set_xscale("log")
        axes[n,1].set_xlim([1e-4, 1.05])
        axes[n,1].set_xticks([1e-3, 1e-2, 1e-1, 1e0])
        
        axes[n,0].tick_params(axis='both', which='major', labelsize=fs_s)
        axes[n,0].tick_params(axis='both', which='minor', labelsize=fs_s)
        axes[n,1].tick_params(axis='both', which='major', labelsize=fs_s)
        axes[n,1].tick_params(axis='both', which='minor', labelsize=fs_s)
        axes[n,2].tick_params(axis='both', which='major', labelsize=fs_s)
        axes[n,2].tick_params(axis='both', which='minor', labelsize=fs_s)
        axes[n,3].tick_params(axis='both', which='major', labelsize=fs_s)
        axes[n,3].tick_params(axis='both', which='minor', labelsize=fs_s)
        
        axes[n,3].set_xscale('log')
        axes[n,3].set_xlim([1e1, 3e3])
        
        xp = 0.02
        yp = 0.01
        
        if n == 0:
            char = "A"
        if n == 1:
            char = "B"
        if n == 2:
            char = "C"
        if n == 3:
            char = "D"
        axes[n,0].text(xp, yp, char+str(1), color="k", rotation=0, ha="left", va="bottom", fontsize=fs_l, transform=axes[n,0].transAxes, bbox=dict(facecolor='white', edgecolor='none', alpha=0.5, boxstyle='round,pad=0.05'))
        axes[n,1].text(xp, yp, char+str(2), color="k", rotation=0, ha="left", va="bottom", fontsize=fs_l, transform=axes[n,1].transAxes, bbox=dict(facecolor='white', edgecolor='none', alpha=0.5, boxstyle='round,pad=0.05'))
        axes[n,2].text(xp, yp, char+str(3), color="k", rotation=0, ha="left", va="bottom", fontsize=fs_l, transform=axes[n,2].transAxes, bbox=dict(facecolor='white', edgecolor='none', alpha=0.5, boxstyle='round,pad=0.05'))
        axes[n,3].text(xp, yp, char+str(4), color="k", rotation=0, ha="left", va="bottom", fontsize=fs_l, transform=axes[n,3].transAxes, bbox=dict(facecolor='white', edgecolor='none', alpha=0.5, boxstyle='round,pad=0.05'))

        if n == 1:
            axes[n,0].set_ylabel(r'Pressure $P$ (bar)', fontsize=fs_l)
        
        # Alpha annotation
        xp = 1.01
        yp = 0.5
        axes[n,3].text(xp, yp, r"$\alpha_\mathrm{c}$ = "+str(alpha), color="k", rotation=90, ha="left", va="center", fontsize=fs_m, transform=axes[n,3].transAxes)

    # Titles
    axes[0,2].set_title(name, color=ga.vol_colors[case_color][5], fontsize=fs_l, x=0.0, y=1.03, pad=1)
    axes[n,0].set_xlabel('Temperature\n'+r'$T$ (K)', fontsize=fs_m) 
    axes[n,1].set_xlabel('Molar concentration\n$X^{\mathrm{i}}_{\mathrm{phase}}$ (non-dim.)', fontsize=fs_m)
    axes[n,2].set_xlabel('Molar mass\n'+r'$\mu$ (g mol$^{-1}$)',fontsize=fs_m)
    axes[n,3].set_xlabel('Specific heat\n'+r'$\widehat{\rm c_p}$ (J K$^{-1}$ mol$^{-1}$)',fontsize=fs_m)

    # Line styles
    axes[0,0].semilogy([0,0],[0,0], color=ga.vol_colors["black_3"], lw=ls_ind, ls="-", label=r'$T_{\rm dew}^i$')
    axes[0,0].semilogy([0,0],[0,0], color=ga.vol_colors["black_3"], lw=ls_ind, ls=":", label=r'$p_\mathrm{sat}$')
    axes[0,0].legend(fontsize=fs_s, ncol=1, loc=1)#, title=r"$\alpha_\mathrm{c}$")
    axes[0,1].semilogy([0,0],[0,0], color=ga.vol_colors["black_3"], lw=ls_ind, ls="-", label=r'$x^i_\mathrm{v}$')
    axes[0,1].semilogy([0,0],[0,0], color=ga.vol_colors["black_3"], lw=ls_ind, ls="--", label=r'$x^i_\mathrm{c}$')
    axes[0,1].semilogy([0,0],[0,0], color=ga.vol_colors["black_2"], lw=ls_ind, ls=":", label=r'$X_\mathrm{v}$')
    axes[0,1].legend(fontsize=fs_s, ncol=1, loc=2)

    # Volatile legend
    xp = 1.0
    yp = 1.0
    for idx, vol in enumerate(vol_list.keys()):
        if vol_list[vol] > 0.:
            axes[1,0].text(xp, yp, ga.vol_latex[vol], color=ga.vol_colors[vol][4], rotation=0, ha="right", va="top", fontsize=fs_l, transform=axes[1,0].transAxes, bbox=dict(facecolor='white', edgecolor='none', alpha=0.5, boxstyle='round,pad=0.05'))
            yp -= 0.1

    sns.despine()
        
    #plt.show()
        
    plt.savefig(setting+'.pdf', bbox_inches='tight')
    #plt.close(fig)  

