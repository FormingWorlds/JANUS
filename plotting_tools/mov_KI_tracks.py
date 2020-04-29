import numpy as np
import math, phys, os, glob, re
import GeneralAdiabat as ga # Moist adiabat with multiple condensibles
import matplotlib.pyplot as plt
import matplotlib
import SocRadModel
from atmosphere_column import atmos
import pandas as pd
from scipy import interpolate
# import seaborn as sns
import copy
import SocRadConv
# from natsort import natsorted # https://pypi.python.org/pypi/natsort
import pickle as pkl
import matplotlib.transforms as mtransforms # https://matplotlib.org/examples/pylab_examples/fancybox_demo.html
from matplotlib.patches import FancyBboxPatch
from random import seed
from random import randint

# Sting sorting not based on natsorted package
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

# Disable and enable print: https://stackoverflow.com/questions/8391411/suppress-calls-to-print-python
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__

def CleanOutputDir( output_dir ):

    types = ("*.json", "*.log", "*.csv", "*.pkl", "current??.????", "profile.*") 
    files_to_delete = []
    for files in types:
        files_to_delete.extend(glob.glob(output_dir+"/"+files))

    # print("Remove old output files:")
    for file in natural_sort(files_to_delete):
        os.remove(file)
    #     print(os.path.basename(file), end =" ")
    # print("\n==> Done.")

########## Read in literature data

def literature_comparison():

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

    ### Plot and annotate literature comparison
    ax1.plot(Goldblatt13_Ts, Goldblatt13_OLR, color=ga.vol_colors["qgray"], ls=":", lw=1.0, zorder=0.1)
    ax1.text(1900, 320, "Goldblatt+ 13", va="bottom", ha="right", fontsize=7, color=ga.vol_colors["qgray"], bbox=dict(fc='white', ec="white", alpha=0.5, pad=0.05, boxstyle='round'))
    # ax1.plot(Kopparapu13_Ts, Kopparapu13_OLR, color=ga.vol_colors["qgray"], ls="-.", lw=1.0, zorder=0.1)
    ax1.plot(Hamano15_Ts, Hamano15_OLR, color=ga.vol_colors["qgray"], ls="-.", lw=1.0, zorder=0.1)
    # ax1.text(2180, 330, "Hamano+ 15", va="top", ha="left", fontsize=7, color=ga.vol_colors["qgray"], bbox=dict(fc='white', ec="white", alpha=0.5, pad=0.05, boxstyle='round'))
    ax1.text(2180, 350, "Kopparapu+ 13 / Hamano+ 15", va="top", ha="left", fontsize=7, color=ga.vol_colors["qgray"], bbox=dict(fc='white', ec="white", alpha=0.5, pad=0.05, boxstyle='round'))
    # ax1.text(1500, 330, "Kopparapu+ 13", va="top", ha="left", fontsize=7, color=ga.vol_colors["qgray"], bbox=dict(fc='white', ec="white", alpha=0.5, pad=0.05, boxstyle='round'))

def define_mixing_ratios(vol_tmp, vol_list):

    # Set current volatile to 1, others to zero
    for vol_curr in vol_list.keys():

        ### Pure cases
        if vol_curr == vol_tmp:
            vol_list[vol_curr] = 1.0
        else:
            vol_list[vol_curr] = 0.0

        ### Mixed  cases
        if "-" in vol_tmp:
            vol_split = vol_tmp.split("-")
            vol1 = vol_split[0]
            vol2 = vol_split[-1]
            if vol_curr == vol1 or vol_curr == vol2:
                vol_list[vol_curr] = 0.5
            else:
                vol_list[vol_curr] = 0.0

    # Color ranges for pure species
    if vol_tmp in vol_list.keys(): 
        vol_color = ga.vol_colors[vol_tmp][5]
    # Specific ones for mixtures
    else:
        try:
            vol_color = ga.vol_colors[vol_tmp]
        except:
            int_value = randint(1, 10) # range of colorscale
            vol_color = ga.vol_colors["mixtures"][int_value]

    return vol_list, vol_color

### Initial conditions

dirs = { "output": os.getcwd()+"/../output/mov_figs", "data_dir": os.getcwd()+"/../output/radiation_limits_data", "rad_conv": os.getcwd()+"/.."}

# Check if data dirs exists, otherwise create
for dat_dir in dirs.values():
    if not os.path.exists(dat_dir):
        os.makedirs(dat_dir)
        print("--> Create directory:", dat_dir)

# General plotting settings

ls_list = [ "-", "--", ":", "-.", (0, (2, 1)), (0, (5, 1)) ]
lw      = 1.5
col_idx = 5

label_fs     = 12
legend_fs    = 10
ticks_fs     = 12
annotate_fs  = 12
numbering_fs = 18

# Define resolution of movie images (low quality: 100; high: 300)
img_dpi = 250

# Set up volatile mixing ratio dict
vol_list    = { 
              "H2O" : .0, 
              "CO2" : .0,
              "H2"  : .0, 
              "N2"  : .0,  
              "CH4" : .0, 
              "O2"  : .0, 
              "CO"  : .0 
            }

##### Parameter ranges considered

# Planet age and orbit
time = { "planet": 0., "star": 4567e+6 } # yr,

# Star age range, yr
star_age_range = [ 0.100e+9 ]          # yr: 0.100e+9, 4.567e+9

# Star mass range, M_sun
Mstar_range = [ 1.0 ]

# Planet-star distance range, au
distance_range = [ 1.0 ]

# Surface pressure range (Pa) for plot A
prs_range    = [ 1e+5 ]

# Surface temperature range (K)
tmp_range   = np.arange(200, 3001, 50)
# tmp_range   = [ 200, 1500, 3000 ] # rapid testing option

##### Define volatile combinations plotted, options: 
#   Single species: "H2O", "CO2", "H2", "CH4", "N2", "CO", "O2"
#   Mixtures: "H2O-CO2", "H2-CO", "H2-CH4", "H2O-H2", "H2-N2", "CO2-N2"
vol_array = [ "H2O", "CO2" ]

##### Define atmosphere settings to be considered
#   Options: 
#       "trpp"  : With tropopause/stratosphere included
#       "moist" : Pure moist adiabat structure 
#       "tstep" : With timestepping
atm_options = [ "trpp" ] # "trpp", "moist", "tstep"

##### Define batch name list
batch_name_list = [ "comp_H2O_1bar", "comp_CH4_1bar", "comp_H2_1bar", "comp_CO2_1bar" ]

# Loop over settings in chain execution
for batch_name in batch_name_list:

    # Special batches for chain execution
    if batch_name ==  "Earth_H2O_1bar":
        vol_array   = [ "H2O", "H2O-H2", "H2O-CO2", "H2O-N2" ]
        atm_options = [ "trpp" ]
        prs_range    = [ 1e+5 ]
    if batch_name ==  "Earth_H2_1bar":
        vol_array   = [ "H2", "H2-H2O", "H2-CH4", "H2-CO", "H2-N2" ]
        atm_options = [ "trpp" ]
        prs_range    = [ 1e+5 ]
    if batch_name ==  "Earth_CO2_1bar":
        vol_array   = [ "CO2", "CO2-H2O", "CO2-H2", "CO2-N2"  ]
        atm_options = [ "trpp" ]
        prs_range    = [ 1e+5 ]
    if batch_name ==  "Earth_CH4_1bar":
        vol_array   = [ "CH4", "CH4-CO", "CH4-H2", "CH4-N2"  ]
        atm_options = [ "trpp" ]
        prs_range    = [ 1e+5 ]
    if batch_name ==  "Earth_H2O_260bar":
        vol_array   = [ "H2O", "H2O-H2", "H2O-CO2", "H2O-N2" ]
        atm_options = [ "trpp" ]
        prs_range    = [ 260e+5 ]
    if batch_name ==  "Earth_H2_260bar":
        vol_array   = [ "H2", "H2-H2O", "H2-CH4", "H2-CO", "H2-N2" ]
        atm_options = [ "trpp" ]
        prs_range    = [ 260e+5 ]
    if batch_name ==  "Earth_CO2_260bar":
        vol_array   = [ "CO2", "CO2-H2O", "CO2-H2", "CO2-N2"  ]
        atm_options = [ "trpp" ]
        prs_range    = [ 260e+5 ]
    if batch_name ==  "Earth_CH4_260bar":
        vol_array   = [ "CH4", "CH4-CO", "CH4-H2", "CH4-N2"  ]
        atm_options = [ "trpp" ]
        prs_range    = [ 260e+5 ]
    if batch_name ==  "comp_H2O_1bar":
        vol_array   = [ "H2O" ]
        atm_options = [ "trpp", "moist", "tstep" ]
        prs_range    = [ 1e+5 ]
        tmp_range   = np.arange(200, 3001, 100)
    if batch_name ==  "comp_H2_1bar":
        vol_array   = [ "H2" ]
        atm_options = [ "trpp", "moist", "tstep" ]
        prs_range    = [ 1e+5 ]
        tmp_range   = np.arange(200, 3001, 100)
    if batch_name ==  "comp_CO2_1bar":
        vol_array   = [ "CO2" ]
        atm_options = [ "trpp", "moist", "tstep" ]
        prs_range    = [ 1e+5 ]
        tmp_range   = np.arange(200, 3001, 100)
    if batch_name ==  "comp_CH4_1bar":
        vol_array   = [ "CH4" ]
        atm_options = [ "trpp", "moist", "tstep" ]
        prs_range    = [ 1e+5 ]
        tmp_range   = np.arange(200, 3001, 100)

    print(">>>>>>>>> Setting: ", batch_name)

    # Add current vol-setting to dirs and create save folder
    dirs["batch_dir"] = dirs["output"] + "/" + batch_name

    # Check if data dirs exists, otherwise create
    for dat_dir in dirs.values():
        if not os.path.exists(dat_dir):
            os.makedirs(dat_dir)
            print("--> Create directory:", dat_dir)

    # Mov fig numbering
    fig_counter = 0

    # Linestyle counter for ax1, ax2
    ls_idx = 0

    # Reset previously plotted data for next batch
    data_dict = {}

    # Format tmp range
    tmp_range   = [ int(round(Ts)) for Ts in tmp_range ]
    tmp_range   = tmp_range[::-1] # movies from high T to low T

    # Loop through different atmosphere settings
    for atm_option in atm_options: 

        print(">>>>> Atmosphere option: ", atm_option)

        # Implement setting
        if atm_option == "trpp":
            trpp    = True
            cp_dry  = False
            set_name = "stratosphere"
        if atm_option == "moist":
            trpp    = False
            cp_dry  = False
            set_name = "moist"
        if atm_option == "tstep":
            trpp    = False
            cp_dry  = True
            set_name = "time-stepped"

        # Loop through volatiles
        #   Options: [ "H2O", "CO2", "H2", "N2", "CH4", "CO", "O2" ]
        for vol_idx, vol in enumerate(vol_array): 

            print(">>>>>", vol)

            # Random seed for color generation
            seed(int(fig_counter+ls_idx))

            # Define mixing ratios and color based on function
            vol_list, vol_color = define_mixing_ratios(vol, vol_list)

            # Loop through surface pressures
            for prs_idx, P_surf in enumerate(prs_range):

                # Loop through distance range
                for dist_idx, dist in enumerate(distance_range): 

                    # Loop through star masses
                    for Mstar_idx, Mstar in enumerate(Mstar_range): 

                        # Loop through star ages
                        for tstar_idx, tstar in enumerate(star_age_range):

                            time["star"] = tstar

                            OLR_array    = []
                            NET_array    = []
                            tmp_array    = []

                            # Define handle to store data
                            handle = str(atm_option) \
                                     + "_"+vol  \
                                     + "_Ps-" + str(round(P_surf/1e+5))  \
                                     + "_d-"+str(dist)  \
                                     + "_Mstar-" + str(Mstar)  \
                                     + "_tstar-" + str(round(tstar/1e+6))
                            handle_legend = ga.vol_latex[vol] + ", " + str(round(P_surf/1e+5)) + " bar, " + str(dist) + " au, " + set_name

                            # Loop through surface temperatures
                            for T_surf in tmp_range:

                                # Purge previous plot settings
                                legendA1_handles = []
                                legendA2_handles = []
                                legendB1_handles = []
                                legendB2_handles = []
                                plt.close("all")

                                # Set up new plot
                                # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,6))
                                # print(plt.subplots(2, 3, figsize=(18,0)))
                                fig, ((ax1, ax3, ax5, ax7), (ax2, ax4, ax6, ax8)) = plt.subplots(2, 4, figsize=(22,10))
                                # sns.set_style("ticks")
                                # sns.despine()

                                # Print old stuff
                                for handle_old in data_dict.keys():
                                    l1o, = ax1.plot(data_dict[handle_old]["TMP"], data_dict[handle_old]["OLR"], color=data_dict[handle_old]["color"], ls=data_dict[handle_old]["ls"], lw=lw, label=data_dict[handle_old]["legend"], zorder=1, alpha=0.5)
                                    legendA1_handles.append(l1o)
                                    l3o, = ax2.plot(data_dict[handle_old]["TMP"], data_dict[handle_old]["NET"], color=data_dict[handle_old]["color"], ls=data_dict[handle_old]["ls"], lw=lw, label=data_dict[handle_old]["legend"], zorder=1, alpha=0.5)
                                    legendB1_handles.append(l3o)

                                print("###", atm_option, ":", vol, "@", round(P_surf/1e+5), "bar,", dist, "au,", Mstar, "M_sun,", round(tstar/1e+6), "Myr,", int(T_surf), "K", end=" ")

                                # Define file name and path
                                file_name = handle  \
                                            + "_Ts-" + str(int(T_surf))  \
                                            + ".pkl"
                                file_path = dirs["data_dir"]+"/"+file_name

                                # If data exists, read it from file
                                if os.path.isfile(file_path):

                                    # Read pickle file
                                    atm_stream = open(file_path, 'rb')
                                    atm = pkl.load(atm_stream)
                                    atm_stream.close()

                                    print("--> Read file:", file_name, end=" ")
                                # Else: compute anew
                                else:

                                    print("--> No file, create:", file_name, end=" ")

                                    # Create atmosphere object
                                    atm = atmos(T_surf, P_surf, vol_list)

                                    # Compute stellar heating
                                    atm.toa_heating = SocRadConv.InterpolateStellarLuminosity(Mstar, time, dist, atm.albedo_pl)

                                    # Compute atmosphere structure and fluxes
                                    atm_dry, atm = SocRadConv.RadConvEqm(dirs, time, atm, [], [], standalone=False, cp_dry=cp_dry, trpp=trpp)

                                    # Use timestepped atm
                                    if atm_option == "tstep":
                                        atm = atm_dry

                                    # Save data to disk
                                    with open(file_path, "wb") as atm_file: 
                                        pkl.dump(atm, atm_file)

                                # Add current T_surf data
                                OLR_array.append(atm.LW_flux_up[0])
                                NET_array.append(atm.net_flux[0])
                                tmp_array.append(T_surf)

                                ############################################################
                                #################### PLOT A ################################
                                ############################################################

                                # labelA = ga.vol_latex[vol] + ", " + str(dist) + " au, " + str(round(P_surf/1e+5)) + " bar, " + set_name

                                l1, = ax1.plot(tmp_array, OLR_array, color=vol_color, ls=ls_list[ls_idx], lw=lw, label=handle_legend)

                                # Fill color and P_surf legends
                                legendA1_handles.append(l1)

                                # Annotate marker at tip
                                ax1.plot(tmp_array[-1], OLR_array[-1], marker='o', ms=3, color=vol_color)
                                
                                # Literature comparison for corresponding correct settings
                                if P_surf == 260e+5 and vol == "H2O":
                                    literature_comparison()

                                legendA1 = ax1.legend(handles=legendA1_handles, loc=2, ncol=1, fontsize=legend_fs, framealpha=0.3)
                                ax1.add_artist(legendA1)

                                ax1.text(0.98, 0.985, 'A', color="k", rotation=0, ha="right", va="top", fontsize=numbering_fs, transform=ax1.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))

                                ax1.set_xlabel(r'Surface temperature, $T_\mathrm{s}$ (K)', fontsize=label_fs)
                                ax1.set_ylabel(r'Outgoing longwave radiation, $F^{\uparrow}_\mathrm{OLR}$ (W m$^{-2}$)', fontsize=label_fs)
                                ax1.set_yscale("log")
                                ax1.set_xlim(left=np.min(tmp_range), right=np.max(tmp_range))
                                ax1.set_ylim(bottom=1e-1, top=1e+7)
                                ax1.set_xticks([np.min(tmp_range), 500, 1000, 1500, 2000, 2500, np.max(tmp_range)])
                                # ax1.set_yticks([1e-10, 1e-5, 1e0, 1e5])
                                ax1.tick_params(axis='both', which='major', labelsize=ticks_fs)
                                ax1.tick_params(axis='both', which='minor', labelsize=ticks_fs)

                                ############################################################
                                #################### PLOT B ################################
                                ############################################################

                                # labelB = ga.vol_latex[vol] + ", " + str(dist) + " au, " + str(round(P_surf/1e+5)) + " bar, " + set_name
                                l3, = ax2.plot(tmp_array, NET_array, color=vol_color, ls=ls_list[ls_idx], lw=lw, label=handle_legend)

                                # Fill legend
                                legendB1_handles.append(l3)

                                # Annotate marker at tip
                                ax2.plot(tmp_array[-1], NET_array[-1], marker='o', ms=3, color=vol_color)
                                
                                # Legend for the main volatiles
                                legendB1 = ax2.legend(handles=legendB1_handles, loc=2, ncol=1, fontsize=legend_fs, framealpha=0.3)
                                ax2.add_artist(legendB1)

                                ax2.set_xlabel(r'Surface temperature, $T_\mathrm{s}$ (K)', fontsize=label_fs)
                                ax2.set_ylabel(r'Net outgoing radiation, $F^{\uparrow}_\mathrm{net}$ (W m$^{-2}$)', fontsize=label_fs)
                                ax2.set_yscale("symlog")
                                ax2.set_xlim(left=np.min(tmp_range), right=np.max(tmp_range))
                                ax2.set_ylim(bottom=-1e+4, top=1e+7)
                                ax2.set_xticks([np.min(tmp_range), 500, 1000, 1500, 2000, 2500, np.max(tmp_range)])
                                ax2.tick_params(axis='both', which='major', labelsize=ticks_fs)
                                ax2.tick_params(axis='both', which='minor', labelsize=ticks_fs)

                                # # Annotate subplot numbering
                                ax2.text(0.98, 0.985, 'B', color="k", rotation=0, ha="right", va="top", fontsize=numbering_fs, transform=ax2.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))

                                # Indicate cooling/heating regimes
                                # ax2.fill_between(tmp_range, 0, +1e+10, alpha=0.05, color="blue")
                                ax2.text(np.max(tmp_range)*0.99, 0.3, "Net cooling", va="bottom", ha="right", fontsize=legend_fs-1, color=ga.vol_colors["qblue_dark"])
                                ax2.axhline([0], color=ga.vol_colors["qgray"], lw=1.0, ls=":")
                                ax2.text(np.max(tmp_range)*0.99, -0.3, "Net heating", va="top", ha="right", fontsize=legend_fs-1, color=ga.vol_colors["qred_dark"])
                                # ax2.fill_between(tmp_range, 0, -1e+10, alpha=0.05, color="red")

                                ##### Atm structure plots

                                # Line settings
                                col_idx  = 3
                                col_vol1 = "H2O"
                                col_vol2 = "N2"
                                col_vol3 = "H2"
                                col_vol4 = "O2"

                                atm_moist = atm

                                ls_moist    = 2.0
                                ls_dry      = 2.0
                                ls_ind      = 1.0

                                ############################################################
                                #################### PLOT C ################################
                                ############################################################

                                # Temperature vs. pressure
                                if cp_dry == False:
                                    adiabat_label = r'General adiabat'
                                else:
                                    adiabat_label = r'Timestepped'
                                ax3.semilogy(atm_moist.tmp,atm_moist.p, color=ga.vol_colors["qgray_dark"], lw=ls_moist, ls="-", label=adiabat_label)

                                # Print active species
                                active_species = r""
                                for vol_tmp in atm_moist.vol_list:
                                    if atm_moist.vol_list[vol_tmp] > 1e-5:
                                        active_species = active_species + ga.vol_latex[vol_tmp] + ", "
                                active_species = active_species[:-2]
                                ax3.text(0.02, 0.02, r"Active species: "+active_species, va="bottom", ha="left", fontsize=10, transform=ax3.transAxes, bbox=dict(fc='white', ec="white", alpha=0.5, boxstyle='round', pad=0.1), color=ga.vol_colors["black_1"] )
                                
                                # For reference p_sat lines
                                T_sat_array    = np.linspace(20,3000,1000) 
                                p_partial_sum  = np.zeros(len(atm.tmp))

                                vol_list_sorted = {k: v for k, v in sorted(atm.vol_list.items(), key=lambda item: item[1])}

                                # Individual species
                                for vol_tmp in vol_list_sorted.keys():

                                    # Only if volatile is present
                                    if atm.vol_list[vol_tmp] > 1e-10:
                                
                                        # Saturation vapor pressure
                                        Psat_array = [ ga.p_sat(vol_tmp, T) for T in T_sat_array ]
                                        ax3.semilogy( T_sat_array, Psat_array, label=r'$p_\mathrm{sat}$'+ga.vol_latex[vol_tmp], lw=ls_ind, ls=":", color=ga.vol_colors[vol_tmp][4])

                                        # Plot partial pressure
                                        ax3.semilogy(atm.tmp, atm.p_vol[vol_tmp], color=ga.vol_colors[vol_tmp][4], lw=ls_ind, ls="-", label=r'$p$'+ga.vol_latex[vol_tmp],alpha=0.99)

                                ax3.text(0.98, 0.985, 'C', color="k", rotation=0, ha="right", va="top", fontsize=numbering_fs, transform=ax3.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))
                                ax3.legend(loc="best", bbox_to_anchor=(0.42, 0.5, 0.5, 0.5)) # (x, y, width, height), https://matplotlib.org/3.2.1/api/_as_gen/matplotlib.pyplot.legend.html
                                ax3.invert_yaxis()
                                ax3.set_xlabel(r'Temperature $T$ (K)', fontsize=label_fs)
                                ax3.set_ylabel(r'Pressure $P$ (Pa)', fontsize=label_fs)
                               
                                # ax3.set_ylim(bottom=atm_moist.ps*1.01, top=1e0)
                                ax3.set_ylim(bottom=atm.ps, top=atm.ptop)
                                ax3.set_ylim(bottom=atm.ps, top=atm.ptop)

                                ############################################################
                                #################### PLOT D ################################
                                ############################################################

                                # Zero line
                                ax4.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)
                                # ax2.axvline(-1e+3, color=ga.vol_colors["qgray_light"], lw=0.5)
                                # ax2.axvline(1e+3, color=ga.vol_colors["qgray_light"], lw=0.5)

                                # LW down
                                ax4.semilogy(atm_moist.LW_flux_down*(-1),atm_moist.pl, color=vol_color, ls=(0, (2, 1)), label=r'$F_\mathrm{LW}^{\downarrow}$')
                                # ls=(0, (3, 1, 1, 1))
                                
                                # SW down
                                ax4.semilogy(atm_moist.SW_flux_down*(-1),atm_moist.pl, color=vol_color, ls=":", label=r'$F_\mathrm{SW}^{\downarrow}$')
                                
                                # Net flux
                                ax4.semilogy(atm_moist.net_flux,atm_moist.pl, color=vol_color, ls="-", lw=2, label=r'$F_\mathrm{net}^{\uparrow}$')

                                # SW up
                                ax4.semilogy(atm_moist.SW_flux_up,atm_moist.pl, color=vol_color, ls="-.", label=r'$F_\mathrm{SW}^{\uparrow}$')
                                
                                # LW up
                                ax4.semilogy(atm_moist.LW_flux_up,atm_moist.pl, color=vol_color, ls=(0, (5, 1)), label=r'$F_\mathrm{LW}^{\uparrow}$')
                                
                                ax4.legend(ncol=5, fontsize=8.5, loc=3)
                                ax4.invert_yaxis()
                                ax4.set_xscale("symlog") # https://stackoverflow.com/questions/3305865/what-is-the-difference-between-log-and-symlog
                                ax4.set_xlabel(r'Outgoing flux $F^{\uparrow}$ (W m$^{-2}$)', fontsize=label_fs)
                                ax4.set_ylabel(r'Pressure $P$ (Pa)', fontsize=label_fs)
                                
                                ax4.set_ylim(bottom=atm.ps, top=atm.ptop)
                                ax4.set_ylim(bottom=atm.ps, top=atm.ptop)

                                ax4.text(0.98, 0.985, 'D', color="k", rotation=0, ha="right", va="top", fontsize=numbering_fs, transform=ax4.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))


                                ############################################################
                                #################### PLOT E: Abundances ####################
                                ############################################################

                                # Abundances
                                for vol_tmp in vol_list_sorted.keys():

                                    # Only if volatile is present
                                    if atm.vol_list[vol_tmp] > 1e-10:

                                        # Plot individual molar concentrations
                                        ax5.semilogy(atm.x_cond[vol_tmp], atm.p, color=ga.vol_colors[vol_tmp][4], lw=ls_ind, ls="--", label=ga.vol_latex[vol_tmp]+" cond.")
                                        ax5.semilogy(atm.x_gas[vol_tmp], atm.p, color=ga.vol_colors[vol_tmp][4], lw=ls_ind, ls="-", label=ga.vol_latex[vol_tmp]+" gas")

                                # Phase molar concentrations
                                ax5.semilogy(atm.xd+atm.xv,atm.p, color=ga.vol_colors["black_2"], lw=ls_moist, ls=":", label=r"Gas phase")

                                ax5.set_ylim(bottom=atm.ps, top=atm.ptop)
                                ax5.set_ylim(bottom=atm.ps, top=atm.ptop)

                                ax5.legend(loc=2, ncol=2)
                                ax5.set_xscale("log")
                                ax5.set_xlim([1e-4, 1.05])
                                ax5.set_xticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0])
                                ax5.set_xticklabels(["$10^{-4}$", "0.001", "0.01", "0.1", "1"])

                                ax5.set_ylabel(r'Pressure $P$ (Pa)', fontsize=label_fs)
                                ax5.set_xlabel(r'Molar concentration, $X^{\mathrm{i}}_{\mathrm{phase}}$', fontsize=label_fs)

                                ax5.text(0.98, 0.985, 'E', color="k", rotation=0, ha="right", va="top", fontsize=numbering_fs, transform=ax5.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))

                                
                                ############################################################
                                #################### PLOT F: Heating #######################
                                ############################################################

                                # Heating versus pressure
                                ax6.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)

                                ax6.plot(atm_moist.LW_heating, atm_moist.p, ls="--", color=vol_color, label=r'LW')
                                ax6.plot(atm_moist.net_heating, atm_moist.p, lw=ls_moist, color=vol_color, label=r'Net')
                                ax6.plot(atm_moist.SW_heating, atm_moist.p, ls=":", color=vol_color, label=r'SW')

                                # Plot tropopause
                                trpp_idx = int(atm_moist.trpp[0])
                                if trpp_idx > 0:
                                    ax6.axhline(atm_moist.pl[trpp_idx], color=vol_color, lw=1.0, ls="-.", label=r'Tropopause')
                                
                                ax6.invert_yaxis()
                                ax6.legend(ncol=4, fontsize=9.0, loc=3)
                                ax6.set_ylabel(r'Pressure $P$ (Pa)', fontsize=label_fs)
                                ax6.set_xlabel(r'Heating (K/day)', fontsize=label_fs)
                                # ax6.set_xscale("log")
                                ax6.set_yscale("log") 
                                x_minmax = np.max([abs(np.min(atm_moist.net_heating[10:])), abs(np.max(atm_moist.net_heating[10:]))])
                                x_minmax = np.max([ 20, x_minmax ])
                                if not math.isnan(x_minmax):
                                    ax6.set_xlim(left=-x_minmax, right=x_minmax)
                                if x_minmax > 1e5:
                                     ax6.set_xscale("symlog") 
                                
                                ax6.set_ylim(bottom=atm.ps, top=atm.ptop)
                                ax6.set_ylim(bottom=atm.ps, top=atm.ptop)

                                ax6.text(0.98, 0.985, 'F', color="k", rotation=0, ha="right", va="top", fontsize=numbering_fs, transform=ax6.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))

                                ############################################################
                                #################### PLOT G: Wavenumber ####################
                                ############################################################

                                # Wavenumber vs. OLR
                                ax7.plot(atm_moist.band_centres, SocRadConv.surf_Planck_nu(atm_moist)/atm_moist.band_widths, color="gray", ls='--', label=str(int(round(atm_moist.ts)))+' K blackbody')
                                ax7.plot(atm_moist.band_centres, atm_moist.net_spectral_flux[:,0]/atm_moist.band_widths, color=vol_color, lw=ls_ind)
                                ax7.set_ylabel(r'Spectral flux density (W m$^{-2}$ cm$^{-1}$)', fontsize=label_fs)
                                ax7.set_xlabel(r'Wavenumber (cm$^{-1}$)', fontsize=label_fs)
                                ax7.legend(loc=2)
                                ymax_plot = 1.2*np.max(atm_moist.net_spectral_flux[:,0]/atm_moist.band_widths)
                                ax7.set_ylim(bottom=0, top=ymax_plot)
                                ax7.set_xlim(left=0, right=np.max(np.where(atm_moist.net_spectral_flux[:,0]/atm_moist.band_widths > ymax_plot/1000., atm_moist.band_centres, 0.)))

                                ax7.text(0.98, 0.985, 'G', color="k", rotation=0, ha="right", va="top", fontsize=numbering_fs, transform=ax7.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))

                                ############################################################
                                #################### PLOT H: Wavelength ####################
                                ############################################################

                                # Wavelength versus OLR log plot
                                OLR_cm_moist = atm_moist.net_spectral_flux[:,0]/atm_moist.band_widths
                                wavelength_moist  = [ 1e+4/i for i in atm_moist.band_centres ] # microns
                                OLR_micron_moist  = [ 1e+4*i for i in OLR_cm_moist ] # microns
                                
                                ax8.plot(wavelength_moist, OLR_micron_moist, color=vol_color, lw=ls_ind)
                                ax8.set_ylabel(r'Spectral flux density (W m$^{-2}$ $\mu$m$^{-1}$)', fontsize=label_fs)
                                ax8.set_xlabel(r'Wavelength $\lambda$ ($\mu$m)', fontsize=label_fs)
                                ax8.set_xscale("log")
                                ax8.set_yscale("log") 
                                ax8.set_xlim(left=0.1, right=100)
                                ax8.set_ylim(bottom=np.max(OLR_micron_moist)*1e-3, top=np.max(OLR_micron_moist))
                                # ax8.set_yticks([1e-10, 1e-5, 1e0, 1e5])
                                ax8.set_xticks([0.1, 0.3, 1, 3, 10, 30, 100])
                                ax8.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30", "100"])

                                ax8.text(0.98, 0.985, 'H', color="k", rotation=0, ha="right", va="top", fontsize=numbering_fs, transform=ax8.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))

                                ############################################################
                                ############################ FINAL FIGURE PLOT SETTINGS ####
                                ############################################################

                                fig_name = batch_name + "_" + "{:04d}".format(fig_counter) + ".png"
                                print("-->", fig_name)
                                plt.savefig( dirs["batch_dir"] + "/" + fig_name, bbox_inches="tight", dpi=img_dpi)
                                # plt.close(fig)

                                fig_counter += 1

                            ##### Clean SOCRATES dir after each run
                            CleanOutputDir( dirs["rad_conv"] )

                            ### Add to "former" lists
                            data_dict[handle] = {}
                            data_dict[handle]["OLR"]    = OLR_array
                            data_dict[handle]["NET"]    = NET_array
                            data_dict[handle]["TMP"]    = tmp_array
                            data_dict[handle]["ls"]     = ls_list[ls_idx]
                            data_dict[handle]["legend"] = handle_legend
                            data_dict[handle]["color"]  = vol_color
                            

                            # Increase linestyle index for subplots A+B
                            ls_idx += 1

        
