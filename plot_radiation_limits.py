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

# Color definitions: https://www.codecademy.com/articles/seaborn-design-ii
no_colors   = 7
vol_colors = {
    "H2O"     : sns.color_palette("Blues", no_colors),
    "CO2"     : sns.color_palette("Reds", no_colors),
    "H2"      : sns.color_palette("Greens", no_colors),
    "N2"      : sns.cubehelix_palette(7),
    "O2"      : sns.light_palette("darkturquoise", no_colors),
    "CH4"     : sns.color_palette("Oranges", no_colors),
    "CO"      : sns.light_palette("#731d1d", no_colors),
    "S"       : sns.light_palette("#EBB434", no_colors),
    "black_1" : "#000000",
    "black_2" : "#323232",
    "black_3" : "#7f7f7f",
    "qgray"          : "#768E95",
    "qgray2"         : "#888888",
    "qblue"          : "#4283A9", # http://www.color-hex.com/color/4283a9
    "qgreen"         : "#62B4A9", # http://www.color-hex.com/color/62b4a9
    "qred"           : "#E6767A",
    "qturq"          : "#2EC0D1",
    "qorange"        : "#ff7f0e",
    "qmagenta"       : "#9A607F",
    "qyellow"        : "#EBB434",
    "qgray_dark"     : "#465559",
    "qblue_dark"     : "#274e65",
    "qgreen_dark"    : "#3a6c65",
    "qred_dark"      : "#b85e61",
    "qturq_dark"     : "#2499a7",
    "qmagenta_dark"  : "#4d303f",
    "qyellow_dark"   : "#a47d24",
    "qgray_light"    : "#acbbbf",
    "qblue_light"    : "#8db4cb",
    "qgreen_light"   : "#a0d2cb",
    "qred_light"     : "#eb9194",
    "qturq_light"    : "#57ccda",
    "qmagenta_light" : "#c29fb2",
    "qyellow_light"  : "#f1ca70",
}

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
prs_range   = [ 260e+5, 10e5  ]

# Surface temperature range (K)
tmp_range   = np.linspace(100, 3000, 20)

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
for vol_idx, vol in enumerate([ "H2O", "CO2", "H2", "N2", "CH4", "CO", "O2" ]): 

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
            atm = SocRadConv.RadConvEqm("./output", time_current, atm, toa_heating, [], [], standalone=False, cp_dry=False)

            OLR_array.append(atm.net_flux[0])

        # OLR
        print(vol, "@", round(P_surf)/1e+5, "bar, OLRs:", OLR_array, "W/m^2")
        if prs_idx == 0:
            l1, = ax1.plot(tmp_range,OLR_array, color=vol_colors[vol][col_idx], ls=ls_list[prs_idx], lw=lw, label=ga.vol_latex[vol])
        else:
            l1, = ax1.plot(tmp_range,OLR_array, color=vol_colors[vol][col_idx], ls=ls_list[prs_idx], lw=lw)
        legendA1_handles.append(l1)
        
        # Add P_surf legend
        if vol_idx == 0: 
            l2, = ax1.plot([0],[0], color="gray", ls=ls_list[prs_idx], lw=lw, label=r"$P_\mathrm{s}$ = "+str(round(P_surf/1e+5))+" bar")
            legendA2_handles.append(l2)


##### PLOT B

## Vary star parameters
P_surf  = 260e+5     # Pa

# Loop through volatiles
for vol_idx, vol in enumerate([ "H2O", "H2" ]): # "H2O", "CO2", "H2", "N2", "CH4", "CO", "O2"

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
            for Tstar_idx, Tstar in enumerate(Tstar_range):

                LW_flux_up_array    = []
                net_flux_array      = []

                toa_heating, star_luminosity = SocRadConv.InterpolateStellarLuminosity(Mstar, time_current, Tstar, distance)

                # Loop through surface temperatures
                for T_surf in tmp_range:

                    # Create atmosphere object
                    atm = atmos(T_surf, P_surf, vol_list)

                    # Compute heat flux
                    atm = SocRadConv.RadConvEqm("./output", time_current, atm, toa_heating, [], [], standalone=False, cp_dry=False)

                    LW_flux_up_array.append(atm.LW_flux_up[0])
                    net_flux_array.append(atm.net_flux[0])

                print(vol, "@", distance, Mstar, Tstar, toa_heating, LW_flux_up_array, net_flux_array)

                l1, = ax2.plot(tmp_range, net_flux_array, color=vol_colors[vol][col_idx-Mstar_idx*2], ls=ls_list[distance_idx], lw=lw, label=ga.vol_latex[vol]+", "+str(Mstar)+" $M_{\odot}$")
                l2, = ax2.plot([0],[0], color=vol_colors["qgray"], ls=ls_list[distance_idx], lw=lw, label=r"$a$ = "+str(distance)+" au")

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

plt.savefig("./output"+'/radiation_limits.pdf', bbox_inches="tight")
plt.close(fig)

