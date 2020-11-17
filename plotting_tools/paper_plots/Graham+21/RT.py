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


### Initial conditions

# Planet age and orbit
time = { "planet": 0., "star": 1e+09 } # yr,

# Star mass, M_sun
Mstar       = 1.0 

# Planet-star distance, au
mean_distance    = 1.0

# Surface pressure range (Pa)
P_surf    = 306e+5

# Surface temperature range (K)
T_surf    = 700

# Volatiles considered
vol_dict    = { 
              "H2O" :  200e+5/P_surf,
              "NH3" :  0.,
              "CO2" :  100e+5/P_surf,
              "CH4" :  0.,
              "CO"  :  0.,
              "O2"  :  0.,
              "N2"  :  1e+5/P_surf,
              "H2"  :  5e+5/P_surf,
            }

ls_list = [ "-", "--", ":", "-." ]
lw      = 1.5
col_idx = 4

# Font sizes 
fs_l = 16
fs_m = 14
fs_s = 12

legend1_handles = []
legend2_handles = []

fig, ax1 = plt.subplots(1, 1, figsize=(9,6))
sns.set_style("ticks")
sns.despine()

dirs =  {
           "output":   os.getcwd()+"/", 
           "rad_conv": "/Users/tim/bitbucket/pcd_couple-interior-atmosphere/atm_rad_conv"
        }

# Create atmosphere object
atm                = atmos(T_surf, P_surf, vol_dict)

# Set fraction of condensate retained in column
atm.alpha_cloud         = 0.

# Stellar heating
atm.toa_heating = SocRadConv.InterpolateStellarLuminosity(Mstar, time, mean_distance, atm.albedo_pl)

print("TOA heating:", round(atm.toa_heating), "W/m^2")

# Compute heat flux
atm_dry, atm_moist = SocRadConv.RadConvEqm(dirs, time, atm, [], [], standalone=True, cp_dry=False, trpp=True) 

ls_moist    = 2.5
ls_dry      = 2.0
ls_ind      = 1.5

fig, (ax2, ax4, ax3) = plt.subplots(1, 3, figsize=(13,4))
# sns.set_style("ticks")
# sns.despine()

# Line settings
col_idx  = 3
col_vol1 = "H2O"
col_vol2 = "N2"
col_vol3 = "H2"
col_vol4 = "O2"

# Fluxes vs. pressure

# Zero line
ax2.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)
# ax2.axvline(-1e+3, color=ga.vol_colors["qgray_light"], lw=0.5)
# ax2.axvline(1e+3, color=ga.vol_colors["qgray_light"], lw=0.5)

# Net thermal flux
ax2.semilogy(atm_moist.LW_flux_up-atm_moist.LW_flux_down,atm_moist.pl, color=ga.vol_colors[col_vol1][col_idx], ls="-", label=r'$F_\mathrm{t}^{\uparrow}$')

# Net flux
ax2.semilogy(atm_moist.net_flux,atm_moist.pl, color=ga.vol_colors[col_vol1][col_idx+2], ls="-", lw=2, label=r'$F_\mathrm{net}^{\uparrow}$')


ax2.legend(ncol=6, fontsize=10, loc=3)
ax2.invert_yaxis()
ax2.set_xscale("symlog") # https://stackoverflow.com/questions/3305865/what-is-the-difference-between-log-and-symlog
ax2.set_xlabel(r'Net outgoing flux $F_\mathrm{net}^{\uparrow}$ (W m$^{-2}$)')
ax2.set_ylabel(r'Pressure $P$ (Pa)')
ax2.set_ylim(top=atm_moist.ptop, bottom=atm_moist.ps)

# Wavenumber vs. OLR
ax3.plot(atm_moist.band_centres, SocRadConv.surf_Planck_nu(atm_moist)/atm_moist.band_widths, color="gray",ls='--',label=str(round(atm_moist.ts))+' K blackbody')
ax3.plot(atm_moist.band_centres, atm_moist.net_spectral_flux[:,0]/atm_moist.band_widths, color=ga.vol_colors[col_vol1][col_idx+1])
ax3.set_ylabel(r'Spectral flux density (W m$^{-2}$ cm$^{-1}$)')
ax3.set_xlabel(r'Wavenumber (cm$^{-1}$)')
ax3.legend(loc=1)
ymax_plot = 1.2*np.max(atm_moist.net_spectral_flux[:,0]/atm_moist.band_widths)
ax3.set_ylim(bottom=0, top=ymax_plot)
ax3.set_xlim(left=0, right=np.max(np.where(atm_moist.net_spectral_flux[:,0]/atm_moist.band_widths > ymax_plot/1000., atm_moist.band_centres, 0.)))
# ax3.set_xlim(left=0, right=30000)


# Heating versus pressure
ax4.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)
ax4.plot(atm_moist.LW_heating, atm_moist.p, ls="--", color=ga.vol_colors[col_vol1][col_idx+1], label=r'LW')
ax4.plot(atm_moist.net_heating, atm_moist.p, lw=2, color=ga.vol_colors[col_vol1][col_idx+1], label=r'Net')
ax4.plot(atm_moist.SW_heating, atm_moist.p, ls=":", color=ga.vol_colors[col_vol1][col_idx+1], label=r'SW')

# Plot tropopause
trpp_idx = int(atm_moist.trpp[0])
if trpp_idx > 0:
    ax4.axhline(atm_moist.pl[trpp_idx], color=ga.vol_colors[col_vol2][col_idx], lw=1.0, ls="-.", label=r'Tropopause')

ax4.invert_yaxis()
ax4.legend(ncol=4, fontsize=10, loc=3)
ax4.set_ylabel(r'Pressure $P$ (Pa)')
ax4.set_xlabel(r'Heating (K/day)')
# ax4.set_xscale("log")
ax4.set_yscale("log") 
x_minmax = np.max([abs(np.min(atm_moist.net_heating[10:])), abs(np.max(atm_moist.net_heating[10:]))])
x_minmax = np.max([ 20, x_minmax ])
if not math.isnan(x_minmax):
    ax4.set_xlim(left=-x_minmax, right=x_minmax)
ax4.set_ylim(top=atm_moist.ptop, bottom=atm_moist.ps)
# ax4.set_yticks([1e-10, 1e-5, 1e0, 1e5])
# ax4.set_xticks([0.1, 0.3, 1, 3, 10, 30, 100])
# ax4.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30", "100"])

# # Wavelength versus OLR log plot
# OLR_cm_moist = atm_moist.LW_spectral_flux_up[:,0]/atm_moist.band_widths
# wavelength_moist  = [ 1e+4/i for i in atm_moist.band_centres ]          # microns
# OLR_micron_moist  = [ 1e+4*i for i in OLR_cm_moist ]                    # microns

# ax4.plot(wavelength_moist, OLR_micron_moist, color=ga.vol_colors[col_vol1][col_idx+1])
# ax4.set_ylabel(r'Spectral flux density (W m$^{-2}$ $\mu$m$^{-1}$)')
# ax4.set_xlabel(r'Wavelength $\lambda$ ($\mu$m)')
# ax4.set_xscale("log")
# ax4.set_yscale("log") 
# ax4.set_xlim(left=0.1, right=100)
# ax4.set_ylim(bottom=1e-20, top=1e5)
# # ax4.set_yticks([1e-10, 1e-5, 1e0, 1e5])
# ax4.set_xticks([0.1, 0.3, 1, 3, 10, 30, 100])
# ax4.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30", "100"])

sns.despine()

plt.savefig(dirs["output"]+"RT.pdf", bbox_inches="tight")
plt.close(fig)
