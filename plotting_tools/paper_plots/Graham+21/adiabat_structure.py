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
import os
#os.chdir('C:/Users/grahamr/Documents/GitHub/soc-rad-conv')
import phys
import GeneralAdiabat as ga
import SocRadModel
#import SocRadConv
#os.chdir('C:/Users/grahamr/Documents/GitHub/soc-rad-conv/plotting_tools/paper_plots/Graham+21')
from atmosphere_column import atmos


### Initial conditions

# Planet age and orbit
time = { "planet": 0., "star": 100e+6 } # yr,

# Star age, yr
Tstar_range = [ 0.100e+9 ]          # yr , 4.567e+9

# Star mass, M_sun
Mstar       = 1.0 

# Planet-star distance, au
distance    = 1.0

# Surface pressure range (Pa)
P_surf    = 1e7

# Surface temperature range (K)
T_surf    = 700

# Volatiles considered
vol_dict    = { 
              "H2O" :  300e+5/P_surf,
              "CO2" :  100e+5/P_surf,
              "NH3" :  0.,
              "CH4" :  0.,
              "CO"  :  0.,
              "O2"  :  0.,
              "N2"  :  1e+5/P_surf,
              "H2"  :  0.,
            }

# # Set up plot
# fig, ax1 = plt.subplots(1, 1, figsize=(7,6))
# fig2, ax2 = plt.subplots(1, 1, figsize=(7,6))
# sns.set_style("ticks")
# sns.despine()

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

# # Loop through volatiles, options: "H2O", "CO2", "H2", "N2", "CH4", "CO", "O2"
# vol_list = [ "H2O", "CO2", "CH4", "O2", "CO", "N2", "H2" ]
# # vol_list = [ "H2O" ]
# for vol_idx, vol in enumerate(reversed(vol_list)): 
# # for vol_idx, vol in enumerate(reversed([ "H2" ])): 

#     # Set current volatile to 1, others to zero
#     for vol1 in vol_dict.keys():
#         if vol1 == vol:
#             vol_dict[vol1] = 1.0
#         else:
#             vol_dict[vol1] = 0.0

# Create atmosphere object
atm                = atmos(T_surf, P_surf, vol_dict)

# Set fraction of condensate retained in column
atm.alpha_cloud         = 0.

# Calculate moist adiabat + condensation
atm                     = ga.general_adiabat(atm)

ls_moist    = 6
ls_dry      = 2.0
ls_ind      = 3

fig, ax1 = plt.subplots(1, 1, figsize=(7,7))#, sharey=True)
# sns.set_style("ticks")
# sns.despine()


# For reference p_sat lines
T_sat_array    = np.linspace(20,3000,1000) 
p_partial_sum  = np.zeros(len(atm.tmp))

vol_list_sorted = {k: v for k, v in sorted(atm.vol_list.items(), key=lambda item: item[1])}

# Individual species
for vol in vol_list_sorted.keys():
# for vol in [ "N2", "CO2", "H2O" ]:

    # Only if volatile is present
    if atm.vol_list[vol] > 1e-10:

      # Plot partial pressures
        #ax1.semilogy(atm.tmpl, atm.pl_vol[vol]/1e+5, color=ga.vol_colors[vol][4], lw=1.5, ls="--",alpha=0.99 , label=r'$p$'+ga.vol_latex[vol]) # 

        # Saturation vapor pressure for given temperature
        #Psat_array = [ ga.p_sat(vol, T)/1e+5 for T in T_sat_array ]
        #ax1.semilogy( T_sat_array, Psat_array, lw=1, ls=":", color=ga.vol_colors[vol][4], label=r'$p_\mathrm{sat}$'+ga.vol_latex[vol])

        if vol == "CO2":
          zo = 5
        else:
          zo = 3
        
        # Plot dew-point temperatures as functions of pressure, given the partial pressure for a given species at a given pressure level
        ax1.semilogy(np.vectorize(ga.Tdew)(vol,atm.pl_vol[vol]),atm.pl/1e+5, color = ga.vol_colors[vol][4], lw = ls_ind, ls='--',alpha=0.99, label=r'$T_{\rm dew}$('+ga.vol_latex[vol]+')', zorder=zo)
                
        # # Sum up partial pressures
        # p_partial_sum += atm.pl_vol[vol]

        # # Plot individual molar concentrations
        # ax2.semilogy(atm.x_cond[vol],atm.p, color=ga.vol_colors[vol][4], lw=ls_ind, ls="--", label=ga.vol_latex[vol]+" cond.")
        # ax2.semilogy(atm.x_gas[vol],atm.p, color=ga.vol_colors[vol][4], lw=ls_ind, ls="-", label=ga.vol_latex[vol]+" gas")
        
# # # Plot sum of partial pressures as check
# ax1.semilogy(atm.tmp, p_partial_sum, color="green", lw=ls_dry, ls="--", label=r'$\sum p^\mathrm{i}$',alpha=0.99)

# General moist adiabat
ax1.semilogy(atm.tmpl, atm.pl/1e+5, color=ga.vol_colors["black_1"], lw=ls_moist, alpha=0.99, zorder=1) # ,label="General adiabat"

# # Phase molar concentrations
# ax2.semilogy(atm.xd+atm.xv,atm.p, color=ga.vol_colors["black_2"], lw=ls_ind, ls=":", label=r"Gas phase")

sns.despine()

fs_l = 16
fs_m = 14
fs_s = 12

ax1.invert_yaxis()
ax1.set_xlabel(r'Temperature $T$ (K)', fontsize=fs_l)
ax1.set_ylabel(r'Pressure $P$ (bar)', fontsize=fs_l)
# ax1.set_title('Adiabats & individual Clausius-Clapeyron slopes', fontsize=fs_l)
l1 = ax1.legend(loc=1, ncol=3, fontsize=fs_s, title="Dew point temperatures")
plt.setp(l1.get_title(),fontsize=fs_s)
ax1.set_xlim([0,np.max(atm.ts)])

ax1.set_ylim(top=atm.ptop/1e+5)
ax1.set_ylim(bottom=atm.ps/1e+5)

ax1.tick_params(axis='both', which='major', labelsize=fs_m)
ax1.tick_params(axis='both', which='minor', labelsize=fs_m)

ax1.text(0.525, 0.270, r'Multi-component adiabat', color="k", rotation=-35.5, ha="center", va="center", fontsize=fs_l, transform=ax1.transAxes)

ax1.text(1.0, 0.04, '1', color="k", rotation=0, ha="center", va="center", fontsize=fs_s, transform=ax1.transAxes, bbox=dict(boxstyle="Round, pad=0.15", fc="white", ec="k", lw=1, alpha=0.99))
ax1.text(0.775, 0.11, '2', color="k", rotation=0, ha="center", va="center", fontsize=fs_s, transform=ax1.transAxes, bbox=dict(boxstyle="Round, pad=0.15", fc="white", ec="k", lw=1, alpha=0.99))
ax1.text(0.21, 0.43, '3', color="k", rotation=0, ha="center", va="center", fontsize=fs_s, transform=ax1.transAxes, bbox=dict(boxstyle="Round, pad=0.15", fc="white", ec="k", lw=1, alpha=0.99))
ax1.text(0.12, 0.830, '4', color="k", rotation=0, ha="center", va="center", fontsize=fs_s, transform=ax1.transAxes, bbox=dict(boxstyle="Round, pad=0.15", fc="white", ec="k", lw=1, alpha=0.99))

# ax1.text(0.95, 0.7, r'$P_\mathrm{surf}$ = '+str(int(atm.ps/1e+5))+" bar\n"+r'$\alpha_\mathrm{c}$ = %.1f'%atm.alpha_cloud, color="k", rotation=0, ha="right", va="center", fontsize=fs_m, transform=ax1.transAxes)
ax1.text(0.835, 0.75, r'$T_\mathrm{s}$ = '+str(int(np.max(atm.ts)))+" K", color="k", rotation=0, ha="center", va="center", fontsize=fs_m, transform=ax1.transAxes)
ax1.text(0.85, 0.70, r'$P_\mathrm{s}$ = '+str(int(atm.ps/1e+5))+" bar", color="k", rotation=0, ha="center", va="center", fontsize=fs_m, transform=ax1.transAxes)

ax1.text(0.808, 0.65, r'$p_\mathrm{s}$N$_2$ = '+str(round(atm.pl_vol["N2"][-1]/1e+5))+" bar", color=ga.vol_colors["N2"][4], rotation=0, ha="center", va="center", fontsize=fs_m, transform=ax1.transAxes)
ax1.text(0.815, 0.60, r'$p_\mathrm{s}$CO$_2$ = '+str(round(atm.pl_vol["CO2"][-1]/1e+5))+" bar", color=ga.vol_colors["CO2"][4], rotation=0, ha="center", va="center", fontsize=fs_m, transform=ax1.transAxes)
ax1.text(0.8145, 0.55, r'$p_\mathrm{s}$H$_2$O = '+str(round(atm.pl_vol["H2O"][-1]/1e+5))+" bar", color=ga.vol_colors["H2O"][4], rotation=0, ha="center", va="center", fontsize=fs_m, transform=ax1.transAxes)
# ax1.text(0.814, 0.50, r'$\alpha_\mathrm{c}$ = %.1f'%atm.alpha_cloud, color="k", rotation=0, ha="center", va="center", fontsize=fs_m, transform=ax1.transAxes)
ax1.text(0.815, 0.50, r'$\alpha_\mathrm{c}$ = %.1f'%atm.alpha_cloud, color="k", rotation=0, ha="center", va="center", fontsize=fs_m, transform=ax1.transAxes)
# fig.suptitle(r'$\alpha$=%.1f'%atm.alpha_cloud)


sns.despine()
plt.show()
plt.savefig(dirs["output"]+"adiabat_structure.pdf", bbox_inches="tight")
#plt.close(fig)
