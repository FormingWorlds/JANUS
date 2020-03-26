import numpy as np
import math,phys
import GeneralAdiabat as ga # Moist adiabat with multiple condensibles
import matplotlib.pyplot as plt
import matplotlib
import SocRadModel
from atmosphere_column import atmos
import pandas as pd
from scipy import interpolate
import seaborn as sns
import copy
import SocRadConv

### Initial conditions

# Planet age and orbit
time_current  = 1e+7                # yr
time_offset   = 1e+9                # yr
mean_distance = 1.0                 # au
toa_heating   = 0.                  # W m^-2

# Star age range, yr
Tstar_range = [ 0.100e+9 ]          # yr , 4.567e+9

# Star mass range, M_sun
Mstar_range = [ 1.0, 0.1 ]

# Planet-star distance range, au
distance_range = [ 1.0, 0.01 ]

# Surface pressure range (Pa)
prs_range   = [ 1e+5 ]

# Surface temperature range (K)
tmp_range   = [ 800 ]

# Volatile molar concentrations: ! must sum to one !
vol_list    = { 
              "H2O" : .0, 
              "CO2" : .0,
              "H2"  : .0, 
              "N2"  : .0,  
              "CH4" : .0, 
              "O2"  : .0, 
              "CO"  : .0 
            }

# Set up plot
fig, ax1 = plt.subplots(1, 1, figsize=(7,6))
sns.set_style("ticks")
sns.despine()

ls_list = [ "-", "--", ":", "-." ]
lw      = 2.0
col_idx = 4

# Font sizes 
fs_l = 16
fs_m = 14
fs_s = 12

legendA1_handles = []
legendA2_handles = []
# legendB1_handles = []
# legendB2_handles = []

##### PLOT A

# Loop through volatiles, options: "H2O", "CO2", "H2", "N2", "CH4", "CO", "O2"
for vol_idx, vol in enumerate(reversed([ "CO2", "N2" ])): 

    # Set current volatile to 1, others to zero
    for vol1 in vol_list.keys():
        if vol1 == vol:
            vol_list[vol1] = 1.0
        else:
            vol_list[vol1] = 0.0

    # Loop through surface pressures
    for prs_idx, P_surf in enumerate(prs_range):

        # Loop through surface temperatures
        for T_surf in tmp_range:

            # Create atmosphere object
            atm = atmos(T_surf, P_surf, vol_list)

            # Compute heat flux
            atm = SocRadConv.RadConvEqm("./output", time_current, atm, toa_heating, [], [], standalone=False, cp_dry=False)

            print(vol, "@", round(P_surf)/1e+5, "bar,", T_surf, "K")

            if prs_idx == 0:
                l1, = ax1.plot(atm.band_centres,atm.LW_spectral_flux_up[:,0]/atm.band_widths, color=ga.vol_colors[vol][col_idx], ls=ls_list[prs_idx], lw=lw, label=ga.vol_latex[vol]) 
                legendA1_handles.append(l1)
            else:
                l1, = ax1.plot(atm.band_centres,atm.LW_spectral_flux_up[:,0]/atm.band_widths, color=ga.vol_colors[vol][col_idx], ls=ls_list[prs_idx], lw=lw) 
            
            # # Add settings legend
            # if vol_idx == 0: 
            #     l2, = ax1.plot(atm.band_centres,SocRadConv.surf_Planck_nu(atm)/atm.band_widths, color=ga.vol_colors["qgray"], ls=ls_list[prs_idx], label=r'T$_\mathrm{s}$ = '+str(round(atm.ts))+' K\nP$_\mathrm{s}$ = '+str(round(P_surf/1e+5))+' bar')
            #     legendA2_handles.append(l2)

# Overplot blackbody
l1a, = ax1.plot(atm.band_centres,SocRadConv.surf_Planck_nu(atm)/atm.band_widths, color=ga.vol_colors["qgray"], lw=lw, ls="--", label=r'Blackbody')
legendA1_handles.append(l1a)

#### PLOT A settings
# Legend for the main volatiles
legendA1 = ax1.legend(handles=legendA1_handles[::-1], loc=1, ncol=1, fontsize=fs_m)
ax1.add_artist(legendA1)
# # Legend for the line styles
# legendA1 = ax1.legend(handles=legendA2_handles, loc=1, ncol=1, fontsize=fs_m)

ax1.text(0.02, 0.98, r'$T_\mathrm{s}$ = '+str(round(atm.ts))+' K\n$P_\mathrm{s}$ = '+str(round(P_surf/1e+5))+' bar', color="k", rotation=0, ha="left", va="top", fontsize=fs_m, transform=ax1.transAxes)

ax1.set_xlim([np.min(atm.band_centres),np.max(atm.band_centres)])
ax1.set_ylabel(r'Spectral exitance (W m$^{-2}$ cm$^{-1}$)', fontsize=fs_l)
ax1.set_xlabel(r'Wavenumber $\tilde{\nu}$ (cm$^{-1}$)', fontsize=fs_l)

ax1.set_xlim(left=0, right=5000)
ax1.set_ylim(bottom=0)

ax1.tick_params(axis='both', which='major', labelsize=fs_s)
ax1.tick_params(axis='both', which='minor', labelsize=fs_s)

plt.savefig("./output"+'/window_showcase.pdf', bbox_inches="tight")
plt.close(fig)

