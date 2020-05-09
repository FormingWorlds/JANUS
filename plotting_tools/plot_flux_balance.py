
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from scipy import interpolate
import seaborn as sns
import copy
import pathlib
import pickle as pkl
import json

try:
    import phys
    import GeneralAdiabat as ga # Moist adiabat with multiple condensibles
    import SocRadModel
    from atmosphere_column import atmos
except:
    import atm_rad_conv.phys as phys
    import atm_rad_conv.GeneralAdiabat as ga
    import atm_rad_conv.SocRadModel as SocRadModel
    from atm_rad_conv.atmosphere_column import atmos

def surf_Planck_nu(atm):
    h   = 6.63e-34
    c   = 3.0e8
    kb  = 1.38e-23
    B   = np.zeros(len(atm.band_centres))
    c1  = 1.191042e-5
    c2  = 1.4387752
    for i in range(len(atm.band_centres)):
        nu      = atm.band_centres[i]
        B[i]    = (c1*nu**3 / (np.exp(c2*nu/atm.ts)-1))
    B   = (1.-atm.albedo_s) * np.pi * B * atm.band_widths/1000.0
    return B

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(13,10))
sns.set_style("ticks")
sns.despine()

# Line settings
col_idx  = 3
col_vol1 = "H2O"
col_vol2 = "N2"
col_vol3 = "H2"
col_vol4 = "O2"


time = 0
atm_file = "/Users/tim/bitbucket/pcd_couple-interior-atmosphere/atm_rad_conv/output/"+str(int(time))+"_atm.pkl"

# Read pickle file
atm_file_stream = open(atm_file,'rb')
atm = pkl.load(atm_file_stream)
atm_file_stream.close()

# print(atm.p)

# Temperature vs. pressure
ax1.semilogy(atm.tmp,atm.p, color=ga.vol_colors[col_vol1][col_idx+1], ls="-", label=r'Moist adiabat')
ax1.legend()
ax1.invert_yaxis()
ax1.set_xlabel(r'Temperature $T$ (K)')
ax1.set_ylabel(r'Pressure $P$ (Pa)')
# ax1.set_ylim(bottom=atm.ps*1.01)
ax1.set_ylim(top=atm.ptop, bottom=atm.ps)

# Print active species
active_species = r""
for vol in atm.vol_list:
    if atm.vol_list[vol] > 1e-5:
        active_species = active_species + ga.vol_latex[vol] + ", "
active_species = active_species[:-2]
ax1.text(0.02, 0.02, r"Active species: "+active_species, va="bottom", ha="left", fontsize=10, transform=ax1.transAxes, bbox=dict(fc='white', ec="white", alpha=0.5, boxstyle='round', pad=0.1), color=ga.vol_colors["black_1"] )

# Fluxes vs. pressure

# Zero line
ax2.axvline(0, color=ga.vol_colors["qgray_light"], lw=0.5)
# ax2.axvline(-1e+3, color=ga.vol_colors["qgray_light"], lw=0.5)
# ax2.axvline(1e+3, color=ga.vol_colors["qgray_light"], lw=0.5)

# LW down
ax2.semilogy(atm.LW_flux_down*(-1),atm.pl, color=ga.vol_colors[col_vol2][col_idx+1], ls="--", label=r'$F_\mathrm{LW}^{\downarrow}$')
# ls=(0, (3, 1, 1, 1))

# SW down
ax2.semilogy(atm.SW_flux_down*(-1),atm.pl, color=ga.vol_colors[col_vol2][col_idx], ls=":", label=r'$F_\mathrm{SW}^{\downarrow}$')

# Net flux
ax2.semilogy(atm.net_flux,atm.pl, color=ga.vol_colors[col_vol1][6], ls="-", lw=2, label=r'$F_\mathrm{net}$')

# SW up
ax2.semilogy(atm.SW_flux_up,atm.pl, color=ga.vol_colors[col_vol1][col_idx], ls=":", label=r'$F_\mathrm{SW}^{\uparrow}$')

# LW up
ax2.semilogy(atm.LW_flux_up,atm.pl, color=ga.vol_colors[col_vol1][col_idx], ls=(0, (5, 1)), label=r'$F_\mathrm{LW}^{\uparrow}$')

ax2.legend(ncol=6, fontsize=10, loc=3)
ax2.invert_yaxis()
ax2.set_xscale("symlog") # https://stackoverflow.com/questions/3305865/what-is-the-difference-between-log-and-symlog
ax2.set_xlabel(r'Outgoing flux $F^{\uparrow}$ (W m$^{-2}$)')
ax2.set_ylabel(r'Pressure $P$ (Pa)')
ax2.set_ylim(top=atm.ptop, bottom=atm.ps)

# Wavenumber vs. OLR
ax3.plot(atm.band_centres, surf_Planck_nu(atm)/atm.band_widths, color="gray",ls='--',label=str(round(atm.ts))+' K blackbody')
ax3.plot(atm.band_centres, atm.net_spectral_flux[:,0]/atm.band_widths, color=ga.vol_colors[col_vol1][col_idx+1])
ax3.set_ylabel(r'Spectral flux density (W m$^{-2}$ cm$^{-1}$)')
ax3.set_xlabel(r'Wavenumber (cm$^{-1}$)')
ax3.legend(loc=1)
ymax_plot = 1.2*np.max(atm.net_spectral_flux[:,0]/atm.band_widths)
ax3.set_ylim(bottom=0, top=ymax_plot)
ax3.set_xlim(left=0, right=np.max(np.where(atm.net_spectral_flux[:,0]/atm.band_widths > ymax_plot/1000., atm.band_centres, 0.)))
# ax3.set_xlim(left=0, right=30000)

# print(np.shape(atm.LW_flux_up))
cff_sum = np.sum(atm.cff[:,:] * atm.LW_flux_up[0], axis=0)
cff_sum = cff_sum / np.trapz(cff_sum,atm.p/np.max(atm.p), axis=0)
print(np.shape(cff_sum))
# Contribution function
print(np.shape(atm.cff))
for i in range(0, atm.nbands, 50):
    ax4.semilogy(atm.cff[i,:], atm.p/np.max(atm.p))
# ax4.semilogy(cff_sum, atm.p/np.max(atm.p))
ax4.set_ylim(1, np.min(atm.p)/np.max(atm.p))


plt.savefig("/Users/tim/bitbucket/pcd_couple-interior-atmosphere/atm_rad_conv/output/"+"cff_comparison.pdf", bbox_inches="tight")
plt.close(fig)