import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from scipy import interpolate
import copy
import pathlib
import pickle as pkl
import json
import glob, re, os
import seaborn as sns
import phys
import GeneralAdiabat as ga
import SocRadModel
import SocRadConv
from atmosphere_column import atmos


#### PLOT SETTINGS

ls_moist    = 2.5
ls_dry      = 2.0
ls_ind      = 1.5

fig, ((Aax1, Aax2, Aax3), (Bax1, Bax2, Bax3)) = plt.subplots(2, 3, figsize=(13,8))
# sns.set_style("ticks")
sns.despine()

# Line settings
col_idx  = 4
col_vol = "H2O"

ls_list = [ "-", "--", ":", "-." ]
lw      = 1.5
fs_l = 16
fs_m = 14
fs_s = 12

dirs =  {
               "output":   os.getcwd()+"/", 
               "rad_conv": "/Users/tim/bitbucket/pcd_couple-interior-atmosphere/atm_rad_conv"
            }

#### LOOP OVER PARAMETERS

for set_idx, setting in enumerate([ "set1", "set2"]):

    # Retained condensate fraction
    for alpha_idx, alpha_cloud in enumerate([ 0.0, 0.5, 1.0 ]):

        ls = ls_list[alpha_idx]

        ### Initial conditions

        if setting == "set1":

            # Planet age and orbit
            time = { "planet": 0., "star": 1e+09 } # yr,

            # Star mass, M_sun
            Mstar           = 1.0 

            # Planet-star distance, au
            mean_distance   = 1.0

            # Surface pressure range (Pa)
            P_surf          = 30.1e+5

            # Surface temperature range (K)
            T_surf          = 273.15

            # Volatiles considered
            vol_dict    = { 
                          "H2O" :  0.01e+5/P_surf,
                          "NH3" :  0.,
                          "CO2" :  29e+5/P_surf,
                          "CH4" :  0.,
                          "CO"  :  0.,
                          "O2"  :  0.,
                          "N2"  :  1e+5/P_surf,
                          "H2"  :  0e+5/P_surf,
                        }

            # Plot axes
            ax1 = Aax1
            ax2 = Aax2
            ax3 = Aax3

            # Plot color
            col_vol = "H2"

        if setting == "set2":
            time = { "planet": 0., "star": 1e+09 } # yr,
            Mstar           = 1.0 
            mean_distance   = 1.0
            P_surf          = 156e+5
            T_surf          = 700
            vol_dict    = { 
                          "H2O" :  100e+5/P_surf,
                          "NH3" :  0.,
                          "CO2" :  100e+5/P_surf,
                          "CH4" :  0.,
                          "CO"  :  0.,
                          "O2"  :  0.,
                          "N2"  :  10e+5/P_surf,
                          "H2"  :  0e+5/P_surf,
                        }

            ax1 = Bax1
            ax2 = Bax2
            ax3 = Bax3
            col_vol = "H2O"

        # Create atmosphere object
        atm                = atmos(T_surf, P_surf, vol_dict)

        # Set fraction of condensate retained in column
        atm.alpha_cloud       = alpha_cloud

        # Stellar heating
        atm.toa_heating = SocRadConv.InterpolateStellarLuminosity(Mstar, time, mean_distance, atm.albedo_pl)

        print("TOA heating:", round(atm.toa_heating), "W/m^2")

        # Compute heat flux
        atm_dry, atm_moist = SocRadConv.RadConvEqm(dirs, time, atm, [], [], standalone=True, cp_dry=False, trpp=True) 


        ## FLUXES
        ax1.semilogy(atm_moist.net_flux,atm_moist.pl, color=ga.vol_colors[col_vol][col_idx], ls=ls, lw=2)
        # ax1.legend(ncol=6, fontsize=10, loc=3)
        ax1.invert_yaxis()
        ax1.set_xscale("symlog") # https://stackoverflow.com/questions/3305865/what-is-the-difference-between-log-and-symlog
        ax1.set_xlabel(r'Net outgoing flux $F_\mathrm{net}^{\uparrow}$ (W m$^{-2}$)')
        ax1.set_ylabel(r'Pressure $P$ (Pa)')
        ax1.set_ylim(top=atm_moist.ptop, bottom=atm_moist.ps)

        

        # HEATING versus pressure
        ax2.plot(atm_moist.net_heating, atm_moist.p, lw=2, color=ga.vol_colors[col_vol][col_idx], label=str(atm.alpha_cloud), ls=ls)
        trpp_idx = int(atm_moist.trpp[0])
        if trpp_idx > 0:
            ax2.axhline(atm_moist.pl[trpp_idx], color=ga.vol_colors[col_vol2][col_idx], lw=1.0, ls=ls, label=r'Tropopause')
        ax2.invert_yaxis()
        # ax2.set_ylabel(r'Pressure $P$ (Pa)')
        ax2.set_xlabel(r'Heating rate (K/day)')
        ax2.set_yscale("log") 
        x_minmax = np.max([abs(np.min(atm_moist.net_heating[10:])), abs(np.max(atm_moist.net_heating[10:]))])
        x_minmax = np.max([ 20, x_minmax ])
        if not math.isnan(x_minmax):
            ax2.set_xlim(left=-x_minmax, right=x_minmax)
        ax2.set_ylim(top=atm_moist.ptop, bottom=atm_moist.ps)


        

        # SPECTRUM
        ax3.plot(atm_moist.band_centres, atm_moist.net_spectral_flux[:,0]/atm_moist.band_widths, color=ga.vol_colors[col_vol][col_idx], ls=ls, label=str(atm.alpha_cloud))
        ax3.set_ylabel(r'Spectral flux density (W m$^{-2}$ cm$^{-1}$)')
        ax3.set_xlabel(r'Wavenumber (cm$^{-1}$)')
        ymax_plot = 1.2*np.max(atm_moist.net_spectral_flux[:,0]/atm_moist.band_widths)
        ax3.set_ylim(bottom=0, top=ymax_plot)
        ax3.set_xlim(left=0, right=np.max(np.where(atm_moist.net_spectral_flux[:,0]/atm_moist.band_widths > ymax_plot/1000., atm_moist.band_centres, 0.)))
        ax3.legend(fontsize=10, loc=1, title=r"$\alpha_\mathrm{c}$")


# Zero line
Aax1.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)
Bax1.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)
Aax2.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)
Bax2.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)

# Subplot labels
lx = 0.98
ly = 0.015
Aax1.text(lx, ly, r'A1', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Aax1.transAxes)
Aax2.text(lx, ly, r'A2', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Aax2.transAxes)
Aax3.text(lx, ly, r'A3', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Aax3.transAxes)
Bax1.text(lx, ly, r'B1', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Bax1.transAxes)
Bax2.text(lx, ly, r'B2', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Bax2.transAxes)
Bax3.text(lx, ly, r'B3', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Bax3.transAxes)

sns.despine()

plt.savefig(dirs["output"]+"RT.pdf", bbox_inches="tight")
plt.close(fig)
