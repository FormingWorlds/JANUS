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
# Stellar heating
toa_heating, star_luminosity = SocRadConv.InterpolateStellarLuminosity(1.0, time_current, time_offset, mean_distance)
toa_heating   = 0.

# Surface pressure range (Pa)
prs_range = [ 260e+5, 10e5  ]

# Surface temperature range (K)
tmp_range = np.linspace(100, 3000, 30)

# Volatile molar concentrations: ! must sum to one !
vol_list = { 
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

legend1_handles = []
legend2_handles = []

# Loop through volatiles
for vol in [ "H2O", "CO2", "H2", "N2", "CO", "O2"  ]: # "H2O", "CO2", "H2", "N2", "CO", "O2"

    # Set current volatile to 1, others to zero
    for vol1 in vol_list.keys():
        if vol1 == vol:
            vol_list[vol1] = 1.0
        else:
            vol_list[vol1] = 0.0

    # Loop through surface pressures
    for prs_idx, P_surf in enumerate(prs_range):

        OLR_array = []

        # Loop through surface temperatures
        for T_surf in tmp_range:

            # Create atmosphere object
            atm = atmos(T_surf, P_surf, vol_list)

            # Compute heat flux
            LW_flux_up, band_centres, LW_spectral_flux_up_per_band_widths = SocRadConv.RadConvEqm("./output", time_current, atm, toa_heating, [], [], standalone=False, cp_dry=False)

            OLR_array.append(LW_flux_up)

        # Print + plot OLR
        print(vol, "@", round(P_surf)/1e+5, "bar, OLRs:", OLR_array, "W/m^2")
        if prs_idx == 0:
            l1, = ax1.plot(tmp_range,OLR_array, color=ga.vol_colors[vol+"_1"], ls=ls_list[prs_idx], lw=lw, label=ga.vol_latex[vol])
        else:
            l1, = ax1.plot(tmp_range,OLR_array, color=ga.vol_colors[vol+"_1"], ls=ls_list[prs_idx], lw=lw)
        
        l2, = ax1.plot([0],[0], color="gray", ls=ls_list[prs_idx], lw=lw, label=r"$P_\mathrm{surf}$ = "+str(round(P_surf/1e+5))+" bar")
        legend1_handles.append(l1)
        legend2_handles.append(l2)

# ax1.legend()

# Create and add legend for the main volatiles
legend1 = ax1.legend(handles=legend1_handles, loc=2)
ax = plt.gca().add_artist(legend1)

# Second legend for the line styles
ax1.legend(handles=legend2_handles, loc=4)

# Plot settings
ax1.set_xlabel(r'Surface temperature $T_\mathrm{surf}$ (K)')
ax1.set_ylabel(r'OLR (W m$^{-2}$)')
ax1.set_yscale("log") 

# ax1.set_xscale("log")
ax1.set_xlim(left=np.min(tmp_range), right=np.max(tmp_range))
ax1.set_xticks([np.min(tmp_range), 500, 1000, 1500, 2000, 2500, np.max(tmp_range)])
# ax1.set_ylim(bottom=1e-20, top=1e5)
# ax1.set_yticks([1e-10, 1e-5, 1e0, 1e5])

plt.savefig("./output"+'/radiation_limits.pdf', bbox_inches="tight")
plt.close(fig)

