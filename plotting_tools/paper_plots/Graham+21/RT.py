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

# fig, ((Aax1, Aax2, Aax3), (Bax1, Bax2, Bax3), (Cax1, Cax2, Cax3)) = plt.subplots(3, 3, figsize=(13,9))
fig, ((Aax1, Bax1, Cax1), (Aax3, Bax3, Cax3)) = plt.subplots(2, 3, figsize=(13,6))
plt.subplots_adjust(top = 0.99, bottom=0.01, hspace=0.2, wspace=0.2)
sns.despine()

# Line settings
col_idx  = 5
col_vol = "H2O"

ls_list = [ "-", "--", ":", "-." ]
lw      = 1.5
fs_l = 15
fs_m = 13
fs_s = 11

dirs =  {
               "output":   os.getcwd()+"/", 
               "rad_conv": "/Users/tim/bitbucket/pcd_couple-interior-atmosphere/atm_rad_conv"
            }

#### LOOP OVER PARAMETERS
for set_idx, setting in enumerate([ "set1", "set2", "set3" ]): # "set1", "set2", "set3"

    # Retained condensate fraction
    for alpha_idx, alpha_cloud in enumerate([ 0.0, 0.1, 1.0 ]): # 0.0, 0.1, 1.0

        ls = ls_list[alpha_idx]

        ### Initial conditions

        # Magma ocean case
        if setting == "set1":
            name = "Magma ocean"
            time = { "planet": 0., "star": 100e+06 } # yr,
            Mstar           = 1.0 
            mean_distance   = 1.0
            P_surf          = "calc"
            T_surf          = 1500
            vol_dict    = { 
                          "H2O" :  500e+5,
                          "NH3" :  0.,
                          "CO2" :  100e+5,
                          "CH4" :  0.,
                          "CO"  :  0.,
                          "O2"  :  0.,
                          "N2"  :  1e+5,
                          "H2"  :  100e+5,
                        }

            ax1 = Aax1
            # ax2 = Cax2
            ax3 = Aax3
            col_vol = "CO2"

        # LATE VENEER case
        if setting == "set2":
            name = "Late veneer"
            time = { "planet": 0., "star": 500e+06 } # yr,
            Mstar           = 1.0 
            mean_distance   = 1.0
            
            P_surf          = "calc"

            ############## Tsurf = 700 K
            T_surf          = 700
            # Maximum late veneer, 2 bar CO2
            vol_dict        = { 
                              "H2O" :  500e+5,
                              "NH3" :  0.07e+5,
                              "CO2" :  1.3e+5,
                              "CH4" :  0.1e+5,
                              "CO"  :  50e+5,
                              "O2"  :  0.e+5,
                              "N2"  :  0.07e+5,
                              "H2"  :  50e+5,
                            }
            # # Maximum late veneer, 50 bar CO2
            # vol_dict        = { 
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
            # vol_dict        = { 
            #                   "H2O" :  500e+5,
            #                   "NH3" :  0.0013e+5,
            #                   "CO2" :  0.3e+5,
            #                   "CH4" :  0.05e+5,
            #                   "CO"  :  0.0008e+5,
            #                   "O2"  :  0e+5,
            #                   "N2"  :  0.3e+5,
            #                   "H2"  :  2e+5,
            #                 }
            # # Vesta, 50 bar CO2
            # vol_dict        = { 
            #                   "H2O" :  500e+5,
            #                   "NH3" :  0.007e+5,
            #                   "CO2" :  48e+5,
            #                   "CH4" :  1e+5,
            #                   "CO"  :  0.05e+5,
            #                   "O2"  :  0e+5,
            #                   "N2"  :  1e+5,
            #                   "H2"  :  8e+5,
            #                 }

            # ############## Tsurf = 500 K
            # T_surf          = 500
            # # Maximum late veneer, 2 bar CO2
            # vol_dict        = { 
            #                   "H2O" :  60e+5,
            #                   "NH3" :  0.07e+5,
            #                   "CO2" :  1.3e+5,
            #                   "CH4" :  0.1e+5,
            #                   "CO"  :  50e+5,
            #                   "O2"  :  0.e+5,
            #                   "N2"  :  0.07e+5,
            #                   "H2"  :  50e+5,
            #                 }
            # # Maximum late veneer, 50 bar CO2
            # vol_dict        = { 
            #                   "H2O" :  60e+5,
            #                   "NH3" :  0.08e+5,
            #                   "CO2" :  0e+5,
            #                   "CH4" :  4e+5,
            #                   "CO"  :  50e+5,
            #                   "O2"  :  0e+5,
            #                   "N2"  :  0.09e+5,
            #                   "H2"  :  50e+5,
            #                 }
            # # Vesta, 2 bar CO2
            # vol_dict        = { 
            #                   "H2O" :  60e+5,
            #                   "NH3" :  0.0013e+5,
            #                   "CO2" :  0.3e+5,
            #                   "CH4" :  0.05e+5,
            #                   "CO"  :  0.0008e+5,
            #                   "O2"  :  0e+5,
            #                   "N2"  :  0.3e+5,
            #                   "H2"  :  2e+5,
            #                 }
            # # Vesta, 50 bar CO2
            # vol_dict        = { 
            #                   "H2O" :  60e+5,
            #                   "NH3" :  0.007e+5,
            #                   "CO2" :  48e+5,
            #                   "CH4" :  1e+5,
            #                   "CO"  :  0.05e+5,
            #                   "O2"  :  0e+5,
            #                   "N2"  :  1e+5,
            #                   "H2"  :  8e+5,
            #                 }


            # # ############## Tsurf = 320 K
            # T_surf          = 320
            # # Maximum late veneer, 2 bar CO2
            # vol_dict        = { 
            #                   "H2O" :  0.1e+5,
            #                   "NH3" :  0.07e+5,
            #                   "CO2" :  1.3e+5,
            #                   "CH4" :  0.1e+5,
            #                   "CO"  :  50e+5,
            #                   "O2"  :  0.e+5,
            #                   "N2"  :  0.07e+5,
            #                   "H2"  :  50e+5,
            #                 }
            # # Maximum late veneer, 50 bar CO2
            # vol_dict        = { 
            #                   "H2O" :  0.1e+5,
            #                   "NH3" :  0.08e+5,
            #                   "CO2" :  0e+5,
            #                   "CH4" :  4e+5,
            #                   "CO"  :  50e+5,
            #                   "O2"  :  0e+5,
            #                   "N2"  :  0.09e+5,
            #                   "H2"  :  50e+5,
            #                 }
            # # Vesta, 2 bar CO2
            # vol_dict        = { 
            #                   "H2O" :  0.1e+5,
            #                   "NH3" :  0.0013e+5,
            #                   "CO2" :  0.3e+5,
            #                   "CH4" :  0.05e+5,
            #                   "CO"  :  0.0008e+5,
            #                   "O2"  :  0e+5,
            #                   "N2"  :  0.3e+5,
            #                   "H2"  :  2e+5,
            #                 }
            # # Vesta, 50 bar CO2
            # vol_dict        = { 
            #                   "H2O" :  0.1e+5,
            #                   "NH3" :  0.007e+5,
            #                   "CO2" :  48e+5,
            #                   "CH4" :  1e+5,
            #                   "CO"  :  0.05e+5,
            #                   "O2"  :  0e+5,
            #                   "N2"  :  1e+5,
            #                   "H2"  :  8e+5,
            #                 }

            ax1 = Bax1
            # ax2 = Bax2
            ax3 = Bax3
            col_vol = "H2O"

        # Hadean Earth / OHZ case
        if setting == "set3":

            name = r"Exo-Earth"

            # Planet age and orbit
            time = { "planet": 0., "star": 4567e+06 } # yr,

            # Star mass, M_sun
            Mstar           = 1.0 

            # Planet-star distance, au
            mean_distance   = 1.5

            # Surface pressure range (Pa)
            P_surf          = "calc"

            # Surface temperature range (K)
            T_surf          = 290

            # Volatiles considered
            vol_dict    = { 
                          "H2O" :  0.01e+5,
                          "NH3" :  0.,
                          "CO2" :  29e+5,
                          "CH4" :  0e+5,
                          "CO"  :  0.,
                          "O2"  :  0.,
                          "N2"  :  1e+5,
                          "H2"  :  0e+5,
                        }

            # Plot axes
            ax1 = Cax1
            # ax2 = Cax2
            ax3 = Cax3

            # Plot color
            col_vol = "H2"


        print("---------", setting, "case:", name, "alpha:", alpha_cloud)

        # Create atmosphere object
        atm                = atmos(T_surf, P_surf, vol_dict)

        # print(P_surf, atm.vol_list)

        # Set fraction of condensate retained in column
        atm.alpha_cloud       = alpha_cloud

        # Stellar heating
        atm.toa_heating = SocRadConv.InterpolateStellarLuminosity(Mstar, time, mean_distance, atm.albedo_pl)

        print("TOA heating:", round(atm.toa_heating), "W/m^2")

        # Compute heat flux
        atm_dry, atm_moist = SocRadConv.RadConvEqm(dirs, time, atm, [], [], standalone=True, cp_dry=False, trpp=True) 


        ### FLUXES
        # ax1.semilogy(atm_moist.net_flux,atm_moist.pl, color=ga.vol_colors[col_vol][col_idx], ls=ls, lw=2)
        ax1.semilogy(atm_moist.flux_up_total,atm_moist.pl/1e+5, color=ga.vol_colors[col_vol][col_idx], ls=ls, lw=1.0)
        # ax1.semilogy(atm_moist.SW_flux_down*(-1),atm_moist.pl, color="k", ls=":", lw=1)
        # ax1.legend(ncol=6, fontsize=10, loc=3)
        ax1.invert_yaxis()
        ax1.set_xscale("symlog") # https://stackoverflow.com/questions/3305865/what-is-the-difference-between-log-and-symlog
        ax1.set_xlabel(r'Outgoing flux $F_\mathrm{t}^{\uparrow}$ (W m$^{-2}$)')
        # ax1.set_ylabel(r'Pressure $P$ (bar)')
        ax1.set_ylim(top=atm_moist.ptop/1e+5, bottom=atm_moist.ps/1e+5)
        
        # Annotate settig name
        ax1.text(0.98, 0.98, name, color=ga.vol_colors[col_vol][col_idx], rotation=0, ha="right", va="top", fontsize=fs_l, transform=ax1.transAxes)


        # ### HEATING versus pressure
        # ax2.plot(atm_moist.net_heating, atm_moist.p, lw=2, color=ga.vol_colors[col_vol][col_idx], label=str(atm.alpha_cloud), ls=ls)
        # trpp_idx = int(atm_moist.trpp[0])
        # if trpp_idx > 0:
        #     ax2.axhline(atm_moist.pl[trpp_idx], color=ga.vol_colors[col_vol][col_idx], lw=1.0, ls=ls, label=r'Tropopause')
        # ax2.invert_yaxis()
        # ax2.set_ylabel(r'Pressure $P$ (Pa)')
        # ax2.set_xlabel(r'Heating rate $\mathcal{H}$ (K/day)')
        # ax2.set_yscale("log") 
        # x_minmax = np.max([abs(np.min(atm_moist.net_heating[10:])), abs(np.max(atm_moist.net_heating[10:]))])
        # # x_minmax = np.max([ 20, x_minmax ])
        # if not math.isnan(x_minmax):
        #     ax2.set_xlim(left=-x_minmax, right=x_minmax)
        # ax2.set_ylim(top=atm_moist.ptop, bottom=atm_moist.ps)


        ### SPECTRUM

        # # Wavenumber
        # ax3.plot(atm_moist.band_centres, atm_moist.flux_up_total[:,0]/atm_moist.band_widths, color=ga.vol_colors[col_vol][col_idx], ls=ls, label=str(atm.alpha_cloud))
        # ax3.set_ylabel(r'Spectral flux density (W m$^{-2}$ cm$^{-1}$)')
        # ax3.set_xlabel(r'Wavenumber (cm$^{-1}$)')
        # ymax_plot = 1.2*np.max(atm_moist.flux_up_total[:,0]/atm_moist.band_widths)
        # ax3.set_ylim(bottom=0, top=ymax_plot)
        # ax3.set_xlim(left=0, right=np.max(np.where(atm_moist.flux_up_total[:,0]/atm_moist.band_widths > ymax_plot/1000., atm_moist.band_centres, 0.)))
        # ax3.legend(fontsize=10, loc=1, title=r"$\alpha_\mathrm{c}$")

        # Wavelength
        spectral_flux_up_total = atm_moist.LW_spectral_flux_up + atm_moist.SW_spectral_flux_up
        OLR_cm_moist = spectral_flux_up_total[:,0]/atm_moist.band_widths
        wavelength_moist  = [ 1e+4/i for i in atm_moist.band_centres ]          # microns
        OLR_micron_moist  = [ 1e+4*i for i in OLR_cm_moist ]                    # microns
        ax3.plot(wavelength_moist, OLR_micron_moist, color=ga.vol_colors[col_vol][col_idx], lw=1.0, label=str(atm.alpha_cloud), ls=ls)
        # ax3.set_ylabel(r'Spectral flux density (W m$^{-2}$ $\mu$m$^{-1}$)')
        ax3.set_xlabel(r'Wavelength $\lambda$ ($\mu$m)')
        ax3.set_xscale("log")
        ax3.set_yscale("log") 
        ax3.set_xlim(left=1, right=30)
        # ax3.set_ylim(bottom=1e-20, top=1e5)
        # ax4.set_yticks([1e-10, 1e-5, 1e0, 1e5])
        ax3.set_xticks([1, 3, 6, 10, 20, 30])
        ax3.set_xticklabels(["1", "3", "6", "10", "20", "30"])
        ax3.legend(fontsize=10, loc=5, title=r"$\alpha_\mathrm{c}$")

# Zero line
# Aax1.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5, zorder=1)
# Bax1.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5, zorder=1)
# Cax1.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5, zorder=1)
# Aax2.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)
# Bax2.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)

# Specific axes ranges
Aax3.set_ylim(bottom=1e-12, top=1e4)
Bax3.set_ylim(bottom=1e-19, top=1e4)
Cax3.set_ylim(bottom=1e-20, top=1e4)

# Subplot labels
lx = 0.98
ly = 0.015
Aax1.text(lx, ly, r'A1', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Aax1.transAxes)
# Aax2.text(lx, ly, r'A2', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Aax2.transAxes)
Aax3.text(lx, ly, r'A2', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Aax3.transAxes)
Bax1.text(lx, ly, r'B1', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Bax1.transAxes)
# Bax2.text(lx, ly, r'B2', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Bax2.transAxes)
Bax3.text(lx, ly, r'B2', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Bax3.transAxes)
Cax1.text(lx, ly, r'C1', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Cax1.transAxes)
# Cax2.text(lx, ly, r'C2', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Cax2.transAxes)
Cax3.text(lx, ly, r'C2', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_l, transform=Cax3.transAxes)

# Axes labels
Aax1.set_ylabel(r'Pressure $P$ (bar)')
Aax3.set_ylabel(r'Spectral flux density (W m$^{-2}$ $\mu$m$^{-1}$)')

sns.despine()

plt.savefig(dirs["output"]+"RT.pdf", bbox_inches="tight")
plt.close(fig)
