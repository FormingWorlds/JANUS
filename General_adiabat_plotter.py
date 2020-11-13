# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 16:50:56 2020

@author: RJ
Plotting script for the general adiabat
"""
import numpy as np
import matplotlib.pyplot as plt
import GeneralAdiabat as ga
from matplotlib import cm

#%%
no_colors   = 7

vol_colors = {
    "H2O"            : cm.get_cmap('PuBu', no_colors)(range(no_colors)),
    "CO2"            : cm.get_cmap("Reds", no_colors)(range(no_colors)),
    "H2"             : cm.get_cmap("Greens", no_colors)(range(no_colors)),
    "N2"             : cm.get_cmap("Purples", no_colors)(range(no_colors)),
    "O2"             : cm.get_cmap("Wistia", no_colors+2)(range(no_colors+2)),
    "CH4"            : cm.get_cmap("RdPu", no_colors)(range(no_colors)),
    "CO"             : cm.get_cmap("pink_r", no_colors)(range(no_colors)),
    "S"              : cm.get_cmap("YlOrBr", no_colors)(range(no_colors)),
    "He"             : cm.get_cmap("Greys", no_colors)(range(no_colors)),
    "NH3"            : cm.get_cmap("cool", no_colors)(range(no_colors)),
    "mixtures"       : cm.get_cmap("Set3", 9)(range(no_colors)),
    "H2O-CO2"        : cm.get_cmap("Set3", 9)(range(no_colors))[1],
    "CO2-H2O"        : cm.get_cmap("Set3", 9)(range(no_colors))[1],
    "H2O-H2"         : cm.get_cmap("Set3", 9)(range(no_colors))[2],
    "H2-H2O"         : cm.get_cmap("Set3", 9)(range(no_colors))[2],
    "H2-CO"          : cm.get_cmap("Set3", 9)(range(no_colors))[3],
    "CO-H2"          : cm.get_cmap("Set3", 9)(range(no_colors))[3],
    "H2-CO2"         : cm.get_cmap("Set3", 9)(range(no_colors))[4],
    "CO2-H2"         : cm.get_cmap("Set3", 9)(range(no_colors))[4],
    "H2-CH4"         : cm.get_cmap("Set3", 9)(range(no_colors))[5],
    "CH4-H2"         : cm.get_cmap("Set3", 9)(range(no_colors))[5],
    "H2-N2"          : cm.get_cmap("Set2", 9)(range(no_colors))[0],
    "N2-H2"          : cm.get_cmap("Set2", 9)(range(no_colors))[0],
    "CO2-N2"         : cm.get_cmap("Set2", 9)(range(no_colors))[1],
    "N2-CO2"         : cm.get_cmap("Set2", 9)(range(no_colors))[1],
    # "H2O"            : sns.color_palette("PuBu", no_colors),
    # "CO2"            : sns.color_palette("Reds", no_colors),
    # "H2"             : sns.color_palette("Greens", no_colors),
    # "N2"             : sns.color_palette("Purples", no_colors), # sns.cubehelix_palette(7)
    # "O2"             : sns.light_palette("darkturquoise", no_colors),
    # "CH4"            : sns.color_palette("RdPu", no_colors),
    # "CO"             : sns.light_palette("#731d1d", no_colors),
    # "S"              : sns.light_palette("#EBB434", no_colors),
    # "He"             : sns.color_palette("Greys", no_colors),
    # "NH3"            : sns.light_palette("teal", no_colors),
    # "mixtures"       : sns.color_palette("Set2", 11),
    # "H2O-CO2"        : sns.color_palette("Set2", 11)[9],
    # "CO2-H2O"        : sns.color_palette("Set2", 11)[9],
    # "H2O-H2"         : sns.color_palette("Set2", 11)[3],
    # "H2-H2O"         : sns.color_palette("Set2", 11)[3],
    # "H2-CO"          : sns.color_palette("Set2", 11)[7],
    # "CO-H2"          : sns.color_palette("Set2", 11)[7],
    # "H2-CO2"         : sns.color_palette("Set2", 11)[8],
    # "CO2-H2"         : sns.color_palette("Set2", 11)[8],
    # "H2-CH4"         : sns.color_palette("Set2", 11)[2],
    # "CH4-H2"         : sns.color_palette("Set2", 11)[2],
    # "H2-N2"          : sns.color_palette("Set2", 11)[4],
    # "N2-H2"          : sns.color_palette("Set2", 11)[4],
    # "CO2-N2"         : sns.color_palette("Set2", 11)[5],
    # "N2-CO2"         : sns.color_palette("Set2", 11)[5],
    "black_1"        : "#000000",
    "black_2"        : "#323232",
    "black_3"        : "#7f7f7f",
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

# Volatile Latex names
vol_latex = {
    "H2O"     : r"H$_2$O",
    "CO2"     : r"CO$_2$",
    "H2"      : r"H$_2$" ,
    "CH4"     : r"CH$_4$",
    "CO"      : r"CO",
    "N2"      : r"N$_2$",
    "S"       : r"S",
    "O2"      : r"O$_2$",
    "He"      : r"He",
    "NH3"     : r"NH$_3$",
    "H2O-CO2" : r"H$_2$O–CO$_2$",
    "H2O-H2"  : r"H$_2$O–H$_2$",
    "H2O-CO"  : r"H$_2$O–CO",
    "H2O-CH4" : r"H$_2$O–CH$_4$",
    "H2O-N2"  : r"H$_2$O–N$_2$",
    "H2O-O2"  : r"H$_2$O–O$_2$",
    "H2-H2O"  : r"H$_2$–H$_2$O",
    "H2-CO"   : r"H$_2$–CO",
    "H2-CH4"  : r"H$_2$–CH$_4$",
    "H2-CO2"  : r"H$_2$–CO$_2$",
    "H2-N2"   : r"H$_2$–N$_2$",
    "H2-O2"   : r"H$_2$-O$_2$",
    "CO2-N2"  : r"CO$_2$–N$_2$",
    "CO2-H2O" : r"CO$_2$–H$_2$O",
    "CO2-CO"  : r"CO$_2$–CO",
    "CO2-CH4"  : r"CO$_2$–CH$_4$",
    "CO2-O2"  : r"CO$_2$–O$_2$",
    "CO2-H2"  : r"CO$_2$–H$_2$",
    "CO-H2O" : r"CO–H$_2$O",
    "CO-CO2" : r"CO–CO$_2$",
    "CO-H2"  : r"CO–H$_2$",
    "CO-CH4" : r"CO–CH$_4$",
    "CO-N2"  : r"CO–N$_2$",
    "CO-O2"  : r"CO–O$_2$",
    "CH4-H2O" : r"CH$_4$–H$_2$O",
    "CH4-CO2" : r"CH$_4$–CO$_2$",
    "CH4-H2"  : r"CH$_4$–H$_2$",
    "CH4-CO"  : r"CH$_4$–CO",
    "CH4-CH4" : r"CH$_4$–CH$_4$",
    "CH4-N2"  : r"CH$_4$–N$_2$",
    "CH4-O2"  : r"CH$_4$–O$_2$",
    "N2-H2O" : r"N$_2$–H$_2$O",
    "N2-CO2" : r"N$_2$–CO$_2$",
    "N2-H2"  : r"N$_2$–H$_2$",
    "N2-CO"  : r"N$_2$–CO",
    "N2-CH4" : r"N$_2$–CH$_4$",
    "N2-N2"  : r"N$_2$–N$_2$",
    "N2-O2"  : r"N$_2$–O$_2$",
    "O2-H2O" : r"O$_2$–H$_2$O",
    "O2-CO2" : r"O$_2$–CO$_2$",
    "O2-H2"  : r"O$_2$–H$_2$",
    "O2-CO"  : r"O$_2$–CO",
    "O2-CH4" : r"O$_2$–CH$_4$",
    "O2-N2"  : r"O$_2$–N$_2$",
    "O2-O2"  : r"O$_2$–O$_2$",
}
#%%
P_surf                  = 10e+5      # Pa
T_surf                  = 280         # K

# Volatile molar concentrations: ! must sum to one !
vol_list = { 
              "H2O" :  0.01,    # 300e+5/P_surf --> specific p_surf
              "CO2" :  0.89,    # 100e+5/P_surf
              "H2"  : .0, 
              "N2"  : 0.1,     # 1e+5/P_surf
              "CH4" : .0, 
              "O2"  : .0, 
              "CO"  : .0, 
              "He"  : .0,
              "NH3" : .0, 
            }
# Create atmosphere object 
#atm                = ga.atmos(T_surf, P_surf, vol_list)

# Calculate moist adiabat + condensation
#atm                     = ga.general_adiabat(atm, rainout=True)

#%%
ls_moist    = 2.5
ls_dry      = 2.0
ls_ind      = 1.5


#fig, (ax1) = plt.subplots(1, 1, figsize=(13,6))

# sns.set_style("ticks")
# sns.despine()
   
# For reference p_sat lines
alpha_list = [0.0,1.0] 
fig, axes = plt.subplots(len(alpha_list), 3, figsize=(13,6))
if np.ndim(axes) == 1:
    axes = np.expand_dims(axes,axis=0)
# Create atmosphere object 
for n,alpha in enumerate(alpha_list):
    
    atm                = ga.atmos(T_surf, P_surf, vol_list)
    
    # Calculate moist adiabat + condensation
    atm.alpha_cloud         = alpha
    atm                     = ga.general_adiabat(atm)
    
    
    T_sat_array    = np.linspace(20,3000,1000) 
    p_partial_sum  = np.zeros(len(atm.tmp))
    
    vol_list_sorted = {k: v for k, v in sorted(atm.vol_list.items(), key=lambda item: item[1])}
    
    # Individual species
    for vol in vol_list_sorted.keys():
    
        # Only if volatile is present
        if atm.vol_list[vol] > 1e-10:
    
            # Saturation vapor pressure for given temperature
            Psat_array = [ ga.p_sat(vol, T) for T in T_sat_array ]
            axes[n,0].semilogy( T_sat_array, Psat_array, label=r'$p_\mathrm{sat}$'+vol_latex[vol], lw=ls_ind, ls=":", color=vol_colors[vol][4])
    
            # Plot partial pressures
            axes[n,0].semilogy(atm.tmp, atm.p_vol[vol], color=vol_colors[vol][4], lw=ls_ind, ls="-", label=r'$p$'+vol_latex[vol],alpha=0.99)
    
            # Sum up partial pressures
            p_partial_sum += atm.p_vol[vol]
    
            # Plot individual molar concentrations
            axes[n,1].semilogy(atm.x_cond[vol],atm.p, color=vol_colors[vol][4], lw=ls_ind, ls="--", label=vol_latex[vol]+" cond.")
            axes[n,1].semilogy(atm.x_gas[vol],atm.p, color=vol_colors[vol][4], lw=ls_ind, ls="-", label=vol_latex[vol]+" gas")
            
    # # Plot sum of partial pressures as check
    # ax1.semilogy(atm.tmp, p_partial_sum, color="green", lw=ls_dry, ls="-", label=r'$\sum p^\mathrm{i}$',alpha=0.99)
    
    # # Dry adiabat function from RTB book
    # ax1.semilogy( dry_adiabat( atm.ts, atm.p, atm.cp ), atm.p , color=vol_colors["black_3"], ls="-.", lw=ls_dry, label=r'Dry adiabat function') # Functional form
    
    # General moist adiabat
    axes[n,0].semilogy(atm.tmp, atm.p, color=vol_colors["black_1"], lw=ls_moist,label="Adiabat",alpha=0.99)
    
    # Phase molar concentrations
    axes[n,1].semilogy(atm.xd+atm.xv,atm.p, color=vol_colors["black_2"], lw=ls_ind, ls=":", label=r"Gas phase")
    
    fs_l = 16
    fs_m = 14
    fs_s = 12
    
    axes[n,0].invert_yaxis()
    axes[n,0].set_ylabel(r'Pressure, $P$ (Pa)', fontsize=fs_l)
    # ax1.set_title('Adiabats & individual Clausius-Clapeyron slopes', fontsize=fs_l)
    if n == 0:
        axes[n,0].legend(loc=1, ncol=np.min([len(atm.vol_list)+1,2]), fontsize=fs_s)
        axes[n,1].legend(loc=2, ncol=2, fontsize=fs_s)

    axes[n,0].set_xlim([0,np.max(atm.ts)])
    
    axes[n,1].invert_yaxis()
    #axes[n,1].set_title('Phase & species abundances', fontsize=fs_l)
    
    #axes[n,1].set_ylabel(r'Pressure, $P$ (Pa)', fontsize=fs_l)
    
    axes[n,0].set_ylim(top=atm.ptop)
    axes[n,0].set_ylim(bottom=atm.ps)
    axes[n,1].set_ylim(top=atm.ptop)
    axes[n,1].set_ylim(bottom=atm.ps)
    
    axes[n,1].set_xscale("log")
    axes[n,1].set_xlim([1e-4, 1.05])
    axes[n,1].set_xticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0])
    if n+1 == len(alpha_list):
        
        axes[n,0].set_xlabel(r'Temperature, $T$ (K)', fontsize=fs_l)

        
        axes[n,1].set_xlabel(r'Molar concentration, $X^{\mathrm{i}}_{\mathrm{phase}}$', fontsize=fs_l)
    axes[n,1].set_xticklabels(["$10^{-4}$", "0.001", "0.01", "0.1", "1"])    
    # ax2.set_xlim(right=1.1)
    
    axes[n,0].tick_params(axis='both', which='major', labelsize=fs_m)
    axes[n,0].tick_params(axis='both', which='minor', labelsize=fs_m)
    axes[n,1].tick_params(axis='both', which='major', labelsize=fs_m)
    axes[n,1].tick_params(axis='both', which='minor', labelsize=fs_m)
    
    axes[n,0].text(0.02, 0.015, 'A', color="k", rotation=0, ha="left", va="bottom", fontsize=fs_l+3, transform=axes[n,0].transAxes)
    axes[n,1].text(0.02, 0.015, 'B', color="k", rotation=0, ha="left", va="bottom", fontsize=fs_l+3, transform=axes[n,1].transAxes)
    axes[n,0].grid()
    axes[n,1].grid()
    plt.show()
    
    #plt.savefig('./output/general_adiabat.pdf', bbox_inches='tight')
    #plt.close(fig)  

