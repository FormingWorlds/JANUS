import numpy as np
import math, phys, os
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

dirs = {"output": os.getcwd()+"/output", "rad_conv": os.getcwd()}

# Planet age and orbit
time = { "planet": 0., "star": 4567e+6 } # yr,
mean_distance = 1.0                 # au

# Star age range, yr
star_age_range = [ 0.100e+9 ]          # yr , 4.567e+9

# Star mass range, M_sun
Mstar_range = [ 1.0 ]

# Planet-star distance range, au
distance_range = [ 1.0 ]

# Surface pressure range (Pa)
prs_range   = [ 260e+5  ]

# Surface temperature range (K)
tmp_range   = np.linspace(200, 3000, 3)

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
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,6))
sns.set_style("ticks")
sns.despine()

ls_list = [ "-", "--", ":", "-." ]
lw      = 2.0
col_idx = 5

legendA1_handles = []
legendA2_handles = []
legendB1_handles = []
legendB2_handles = []

##### PLOT A

# Loop through volatiles, options: "H2O", "CO2", "H2", "N2", "CH4", "CO", "O2"
for vol_idx, vol in enumerate([ "CO2" ]): 

    # Set current volatile to 1, others to zero
    for vol1 in vol_list.keys():
        if vol1 == vol:
            vol_list[vol1] = 1.0
        else:
            vol_list[vol1] = 1e-10

    # Loop through surface pressures
    for prs_idx, P_surf in enumerate(prs_range):

        OLR_array = []

        # Loop through surface temperatures
        for T_surf in tmp_range:

            print("###", vol, P_surf, T_surf)

            # Create atmosphere object
            atm = atmos(T_surf, P_surf, vol_list)

            # Compute heat flux
            atm = SocRadConv.RadConvEqm(dirs, time, atm, [], [], standalone=False, cp_dry=False)

            OLR_array.append(atm.net_flux[0])

        # OLR
        print(vol, "@", round(P_surf)/1e+5, "bar, OLRs:", OLR_array, "W/m^2")
        if prs_idx == 0:
            l1, = ax1.plot(tmp_range,OLR_array, color=ga.vol_colors[vol][col_idx], ls=ls_list[prs_idx], lw=lw, label=ga.vol_latex[vol])
        else:
            l1, = ax1.plot(tmp_range,OLR_array, color=ga.vol_colors[vol][col_idx], ls=ls_list[prs_idx], lw=lw)
        legendA1_handles.append(l1)
        
        # Add P_surf legend
        if vol_idx == 0: 
            l2, = ax1.plot([0],[0], color="gray", ls=ls_list[prs_idx], lw=lw, label=r"$P_\mathrm{s}$ = "+str(round(P_surf/1e+5))+" bar")
            legendA2_handles.append(l2)


##### PLOT B

## Vary star parameters
P_surf  = 260e+5     # Pa

# Loop through volatiles
for vol_idx, vol in enumerate([ "H2O", "CO2" ]): # "H2O", "CO2", "H2", "N2", "CH4", "CO", "O2"

    # Set current volatile to 1, others to zero
    for vol1 in vol_list.keys():
        if vol1 == vol:
            vol_list[vol1] = 1.0
        else:
            vol_list[vol1] = 0.0

    # Loop through distance range
    for distance_idx, distance in enumerate(distance_range): 

        # Loop through star masses
        for Mstar_idx, Mstar in enumerate(Mstar_range): 

            # Set current volatile to 1, others to zero
            for vol1 in vol_list.keys():
                if vol1 == vol:
                    vol_list[vol1] = 1.0
                else:
                    vol_list[vol1] = 0.0

            # Loop through star ages
            for star_age_idx, star_age in enumerate(star_age_range):

                print("###", vol, distance, star_age)

                time["star"] = star_age

                LW_flux_up_array    = []
                net_flux_array      = []

                atm.toa_heating = SocRadConv.InterpolateStellarLuminosity(Mstar, time, distance, atm.albedo_pl)

                # Loop through surface temperatures
                for T_surf in tmp_range:

                    # Create atmosphere object
                    atm = atmos(T_surf, P_surf, vol_list)

                    # Compute heat flux
                    atm = SocRadConv.RadConvEqm(dirs, time, atm, [], [], standalone=False, cp_dry=False)

                    LW_flux_up_array.append(atm.LW_flux_up[0])
                    net_flux_array.append(atm.net_flux[0])

                print(vol, "@", distance, Mstar, star_age/1e+6, atm.toa_heating, LW_flux_up_array, net_flux_array)

                l1, = ax2.plot(tmp_range, net_flux_array, color=ga.vol_colors[vol][col_idx-Mstar_idx*2], ls=ls_list[distance_idx], lw=lw, label=ga.vol_latex[vol]+", "+str(Mstar)+" $M_{\odot}$")
                l2, = ax2.plot([0],[0], color=ga.vol_colors["qgray"], ls=ls_list[distance_idx], lw=lw, label=r"$a$ = "+str(distance)+" au")

                if distance_idx == 0: legendB1_handles.append(l1)
                if Mstar_idx == 0 and vol_idx == 0: legendB2_handles.append(l2)


##### PLOT A settings
# Legend for the main volatiles
legendA1 = ax1.legend(handles=legendA1_handles, loc=2, ncol=2, fontsize=10)
ax1.add_artist(legendA1)
# Legend for the line styles
legendA1 = ax1.legend(handles=legendA2_handles, loc=4, ncol=1, fontsize=10)

ax1.set_xlabel(r'Surface temperature, $T_\mathrm{s}$ (K)')
ax1.set_ylabel(r'Net outgoing radiation, $F^{\uparrow}_\mathrm{net}$ (W m$^{-2}$)')
ax1.set_yscale("symlog")
ax1.set_xlim(left=np.min(tmp_range), right=np.max(tmp_range))
ax1.set_xticks([np.min(tmp_range), 500, 1000, 1500, 2000, 2500, np.max(tmp_range)])
# ax1.set_ylim(bottom=1e-20, top=1e5)
# ax1.set_yticks([1e-10, 1e-5, 1e0, 1e5])


##### PLOT B settings
# Legend for the main volatiles
legendB1 = ax2.legend(handles=legendB1_handles, loc=2, ncol=2, fontsize=10)
ax2.add_artist(legendB1)
# Legend for the line styles
legendB1 = ax2.legend(handles=legendB2_handles, loc=4, ncol=1, fontsize=10)


ax2.set_xlabel(r'Surface temperature, $T_\mathrm{s}$ (K)')
ax2.set_ylabel(r'Net outgoing radiation, $F^{\uparrow}_\mathrm{net}$ (W m$^{-2}$)')
ax2.set_yscale("symlog")
ax2.set_xlim(left=np.min(tmp_range), right=np.max(tmp_range))
ax2.set_xticks([np.min(tmp_range), 500, 1000, 1500, 2000, 2500, np.max(tmp_range)])

plt.savefig(dirs["output"]+'/radiation_limits.pdf', bbox_inches="tight")
plt.close(fig)

