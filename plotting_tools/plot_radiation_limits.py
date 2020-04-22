import numpy as np
import math, phys, os, glob
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
from natsort import natsorted # https://pypi.python.org/pypi/natsort
import pickle as pkl
import matplotlib.transforms as mtransforms # https://matplotlib.org/examples/pylab_examples/fancybox_demo.html
from matplotlib.patches import FancyBboxPatch

def CleanOutputDir( output_dir ):

    types = ("*.json", "*.log", "*.csv", "*.pkl", "current??.????", "profile.*") 
    files_to_delete = []
    for files in types:
        files_to_delete.extend(glob.glob(output_dir+"/"+files))

    # print("Remove old output files:")
    for file in natsorted(files_to_delete):
        os.remove(file)
    #     print(os.path.basename(file), end =" ")
    # print("\n==> Done.")

### Initial conditions

dirs = {"output": os.getcwd()+"/../output", "data_dir": os.getcwd()+"/../output/radiation_limits_data", "rad_conv": os.getcwd()+"/.."}

# Check if data dirs exists, otherwise create
if not os.path.exists(dirs["data_dir"]):
    os.makedirs(dirs["data_dir"])
    print("--> Create data directory:", dirs["data_dir"])

# Planet age and orbit
time = { "planet": 0., "star": 4567e+6 } # yr,

# Star age range, yr
star_age_range = [ 4.567e+9 ]          # yr: 0.100e+9, 4.567e+9

# Star mass range, M_sun
Mstar_range = [ 1.0 ]

# Planet-star distance range, au
distance_range = [ 1.0, 0.5 ]

# Surface pressure range (Pa)
prs_range   = [ 260e+5, 1e+5 ]

# Surface temperature range (K)
tmp_range   = np.linspace(200, 3000, 100)

# With / without stratosphere?
trpp        = True

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
lw      = 1.5
col_idx = 5

legendA1_handles = []
legendA2_handles = []
legendB1_handles = []
legendB2_handles = []

##### PLOT A
print("############# PLOT A #############")

# Loop through volatiles, options: "H2O", "CO2", "H2", "N2", "CH4", "CO", "O2"
for vol_idx, vol in enumerate([ "H2O", "CO2", "H2", "CH4" ]): 

    # Set current volatile to 1, others to zero
    for vol1 in vol_list.keys():
        
        # Pure cases
        if vol1 == vol:
            vol_list[vol1] = 1.0
        else:
            vol_list[vol1] = 0.0
        if vol == "H2O-CO2":
            if vol1 == "H2O" or vol1 == "CO2":
                vol_list[vol1] = 0.5
            else:
                vol_list[vol1] = 0.0
        if vol == "H2-CO":
            if vol1 == "H2" or vol1 == "CO":
                vol_list[vol1] = 0.5
            else:
                vol_list[vol1] = 0.0
        if vol == "H2-CH4":
            if vol1 == "H2" or vol1 == "CH4":
                vol_list[vol1] = 0.5
            else:
                vol_list[vol1] = 0.0
        if vol == "H2O-H2":
            if vol1 == "H2O" or vol1 == "H2":
                vol_list[vol1] = 0.5
            else:
                vol_list[vol1] = 0.0
        if vol == "H2-N2":
            if vol1 == "H2" or vol1 == "N2":
                vol_list[vol1] = 0.5
            else:
                vol_list[vol1] = 0.0

    # Loop through surface pressures
    for prs_idx, P_surf in enumerate(prs_range):

        # Check if data present, otherwise create
        file_name = "a_"+vol+"_Ps"+str(round(P_surf))+"_p"+str(np.size(prs_range))+"_T"+str(np.size(tmp_range))+".pkl"
        file_path = dirs["data_dir"]+"/"+file_name

        # If data exists, read it from file
        if os.path.isfile(file_path):

            # Read pickle file
            a_dict_stream = open(file_path, 'rb')
            a_dict = pkl.load(a_dict_stream)
            a_dict_stream.close()

            print("--> Read in file:", file_name)

        # Else: calculate it
        else:

            OLR_array = []

            # Loop through surface temperatures
            for T_surf in tmp_range:

                print("###", vol, P_surf, T_surf)

                # Create atmosphere object
                atm = atmos(T_surf, P_surf, vol_list)

                # Compute heat flux
                atm_dry, atm = SocRadConv.RadConvEqm(dirs, time, atm, [], [], standalone=False, cp_dry=False, trpp=trpp)

                # OLR !
                OLR_array.append(atm.LW_flux_up[0])

            # Save data to disk
            a_dict = { vol: OLR_array, "tmp_range": tmp_range }
            with open(file_path, "wb") as dict_file: 
                pkl.dump(a_dict, dict_file)

        # OLR
        print(vol, "@", round(P_surf)/1e+5, "bar, OLRs:", a_dict[vol], "W/m^2")
        if prs_idx == 0:
            l1, = ax1.plot(a_dict["tmp_range"], a_dict[vol], color=ga.vol_colors[vol][col_idx], ls=ls_list[prs_idx], lw=lw, label=ga.vol_latex[vol])
        else:
            l1, = ax1.plot(a_dict["tmp_range"], a_dict[vol], color=ga.vol_colors[vol][col_idx], ls=ls_list[prs_idx], lw=lw)
        legendA1_handles.append(l1)
        
        # Add P_surf legend
        if vol_idx == 0: 
            l2, = ax1.plot([0],[0], color="gray", ls=ls_list[prs_idx], lw=lw, label=r"$P_\mathrm{s}$ = "+str(round(P_surf/1e+5))+" bar")
            legendA2_handles.append(l2)

        # Clean SOCRATES dir
        CleanOutputDir( dirs["rad_conv"] )

##### PLOT B
print("############# PLOT B #############")

## Vary star parameters
P_surf  = 10e+5    # Pa

b_ymax = 0
b_ymin = 0

# Loop through volatiles, options: "H2O", "CO2", "H2", "N2", "CH4", "CO", "O2"
for vol_idx, vol in enumerate([ "H2O", "CO2", "H2", "CH4" ]):

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

                # Check if data present, otherwise create
                file_name = "b_"+vol+"_Ps"+str(round(P_surf))+"_d"+str(distance)+"_Mstar"+str(Mstar)+"_tstar"+str(np.size(star_age))+"_nT"+str(np.size(tmp_range))+".pkl"
                file_path = dirs["data_dir"]+"/"+file_name

                # If data exists, read it from file
                if os.path.isfile(file_path):

                    # Read pickle file
                    b_dict_stream = open(file_path, 'rb')
                    b_dict = pkl.load(b_dict_stream)
                    b_dict_stream.close()

                    print("--> Read in file:", file_name)

                else:

                    time["star"] = star_age

                    LW_flux_up_array    = []
                    net_flux_array      = []

                    # Loop through surface temperatures
                    for T_surf in tmp_range:

                        print("###", vol, distance, star_age, T_surf)

                        # Create atmosphere object
                        atm = atmos(T_surf, P_surf, vol_list)

                        atm.toa_heating = SocRadConv.InterpolateStellarLuminosity(Mstar, time, distance, atm.albedo_pl)

                        # Compute heat flux
                        atm_dry, atm = SocRadConv.RadConvEqm(dirs, time, atm, [], [], standalone=False, cp_dry=False, trpp=trpp)

                        LW_flux_up_array.append(atm.LW_flux_up[0])
                        net_flux_array.append(atm.net_flux[0])

                    # Save data to disk
                    b_dict = { vol: net_flux_array, "tmp_range": tmp_range }
                    with open(file_path, "wb") as dict_file: 
                        pkl.dump(b_dict, dict_file)

                print(vol, "@", distance, Mstar, star_age/1e+6, b_dict[vol])

                l1, = ax2.plot(b_dict["tmp_range"], b_dict[vol], color=ga.vol_colors[vol][col_idx-Mstar_idx*2], ls=ls_list[distance_idx], lw=lw, label=ga.vol_latex[vol])
                l2, = ax2.plot([0],[0], color=ga.vol_colors["qgray"], ls=ls_list[distance_idx], lw=lw, label=r"$a$ = "+str(distance)+" au")

                if distance_idx == 0: legendB1_handles.append(l1)
                if Mstar_idx == 0 and vol_idx == 0: legendB2_handles.append(l2)

                # Set ylim range
                b_ymin = np.min([ b_ymin, np.min(b_dict[vol]) ])
                b_ymax = np.max([ b_ymax, np.max(b_dict[vol]) ])

        # Clean SOCRATES dir
        CleanOutputDir( dirs["rad_conv"] )

########## Read in literature data

Goldblatt13_Ts  = []
Goldblatt13_OLR = []
with open(dirs["rad_conv"]+"/plotting_tools/comparison_data/Goldblatt13_data.txt", 'r') as data_file:
    for line in data_file:
        if not line.startswith('#'):
            line = line.rstrip('\n')
            line = line.split(",")
            Goldblatt13_Ts.append(float(line[0]))
            Goldblatt13_OLR.append(float(line[1]))
Kopparapu13_Ts  = []
Kopparapu13_OLR = []
with open(dirs["rad_conv"]+"/plotting_tools/comparison_data/Kopparapu13_data.txt", 'r') as data_file:
    for line in data_file:
        if not line.startswith('#'):
            line = line.rstrip('\n')
            line = line.split(",")
            Kopparapu13_Ts.append(float(line[0]))
            Kopparapu13_OLR.append(float(line[1]))
Hamano15_Ts  = []
Hamano15_OLR = []
with open(dirs["rad_conv"]+"/plotting_tools/comparison_data/Hamano15_data.txt", 'r') as data_file:
    for line in data_file:
        if not line.startswith('#'):
            line = line.rstrip('\n')
            line = line.split(",")
            Hamano15_Ts.append(float(line[0]))
            Hamano15_OLR.append(float(line[1]))

########## GENERAL PLOT SETTINGS

label_fs    = 12
legend_fs   = 11
ticks_fs    = 10
annotate_fs = 14

##### PLOT A settings
# Legend for the main volatiles
legendA1 = ax1.legend(handles=legendA1_handles, loc=2, ncol=2, fontsize=legend_fs)
ax1.add_artist(legendA1)
# Legend for the line styles
legendA2 = ax1.legend(handles=legendA2_handles, loc=4, ncol=2, fontsize=legend_fs)

ax1.set_xlabel(r'Surface temperature, $T_\mathrm{s}$ (K)', fontsize=label_fs)
ax1.set_ylabel(r'Outgoing longwave radiation, $F^{\uparrow}_\mathrm{OLR}$ (W m$^{-2}$)', fontsize=label_fs)
ax1.set_yscale("log")
ax1.set_xlim(left=np.min(tmp_range), right=np.max(tmp_range))
ax1.set_xticks([np.min(tmp_range), 500, 1000, 1500, 2000, 2500, np.max(tmp_range)])
# ax1.set_ylim(bottom=1e-20, top=1e5)
# ax1.set_yticks([1e-10, 1e-5, 1e0, 1e5])

##### PLOT B settings
# Legend for the main volatiles
legendB1 = ax2.legend(handles=legendB1_handles, loc=2, ncol=2, fontsize=legend_fs)
ax2.add_artist(legendB1)
# Legend for the line styles
legendB2 = ax2.legend(handles=legendB2_handles, loc=4, ncol=2, fontsize=legend_fs)

ax2.set_xlabel(r'Surface temperature, $T_\mathrm{s}$ (K)', fontsize=label_fs)
ax2.set_ylabel(r'Net outgoing radiation, $F^{\uparrow}_\mathrm{net}$ (W m$^{-2}$)', fontsize=label_fs)
ax2.set_yscale("symlog")
ax2.set_xlim(left=np.min(tmp_range), right=np.max(tmp_range))
ax2.set_ylim(bottom=b_ymin*2., top=b_ymax*10)
ax2.set_xticks([np.min(tmp_range), 500, 1000, 1500, 2000, 2500, np.max(tmp_range)])

# Annotate surface pressure
ax2.text(0.98, 0.98, r'$P_{\mathrm{s}} = $'+str(round(P_surf/1e+5))+" bar", va="top", ha="right", fontsize=annotate_fs, transform=ax2.transAxes, bbox=dict(fc='white', ec="white", alpha=0.9, boxstyle='round', pad=0.1), color=ga.vol_colors["black_1"] )

ax1.text(0.02, 0.015, 'A', color="k", rotation=0, ha="left", va="bottom", fontsize=20, transform=ax1.transAxes, bbox=dict(fc='white', ec="white", alpha=0.7, pad=0.1, boxstyle='round'))
ax2.text(0.02, 0.015, 'B', color="k", rotation=0, ha="left", va="bottom", fontsize=20, transform=ax2.transAxes, bbox=dict(fc='white', ec="white", alpha=0.7, pad=0.1, boxstyle='round'))

# Indicate cooling/heating regimes
ax2.fill_between(tmp_range, 0, +1e+10, alpha=0.05, color="blue")
ax2.text(np.max(tmp_range)*0.99, 0.3, "Net cooling", va="bottom", ha="right", fontsize=legend_fs, color=ga.vol_colors["qblue_dark"])
ax2.text(np.max(tmp_range)*0.99, -0.3, "Net heating", va="top", ha="right", fontsize=legend_fs, color=ga.vol_colors["qred_dark"])
ax2.fill_between(tmp_range, 0, -1e+10, alpha=0.05, color="red")

### Plot and annotate literature comparison
ax1.plot(Goldblatt13_Ts, Goldblatt13_OLR, color=ga.vol_colors["qgray"], ls=":", lw=1.0, zorder=0.1)
ax1.text(1900, 320, "Goldblatt+ 13", va="bottom", ha="right", fontsize=legend_fs-3, color=ga.vol_colors["qgray"])
# ax1.plot(Kopparapu13_Ts, Kopparapu13_OLR, color=ga.vol_colors["qgray"], ls="-.", lw=1.0, zorder=0.1)
ax1.plot(Hamano15_Ts, Hamano15_OLR, color=ga.vol_colors["qgray"], ls="-.", lw=1.0, zorder=0.1)
ax1.text(2180, 330, "Hamano+ 15", va="top", ha="left", fontsize=legend_fs-3, color=ga.vol_colors["qgray"])


plt.savefig(dirs["output"]+'/radiation_limits.pdf', bbox_inches="tight")
plt.close(fig)

