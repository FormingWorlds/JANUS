#!/usr/bin/env python3

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt

import time as t
import os, shutil
import numpy as np

from modules.stellar_luminosity import InterpolateStellarLuminosity
from modules.radcoupler import RadConvEqm

from utils.atmosphere_column import atmos
import utils.StellarSpectrum as StellarSpectrum


def run_once(T_surf, dirs):

    # Planet 
    time = { "planet": 0., "star": 4567e+6 } # yr,
    star_mass     = 1.0                 # M_sun, mass of star
    mean_distance = 1.0                 # au, orbital distance
    pl_radius     = 6.371e6             # m, planet radius
    pl_mass       = 5.972e24            # kg, planet mass

    # Boundary conditions for pressure & temperature
    P_top         = 1.0                  # Pa

    # Define volatiles by mole fractions
    P_surf       = 300 * 1e5
    vol_mixing = {
                    "H2O" : 1.0,
                    "CO2" : 0.0,
                    "N2"  : 0.0
                }
    
    # Rayleigh scattering on/off
    rscatter = False

    # Compute contribution function
    calc_cf = False

    # Tropopause calculation
    trppD = False   # Calculate dynamically?
    trppT = 0.0     # Fixed tropopause value if not calculated dynamically
    
    # Instellation scaling | 1.0 == no scaling
    Sfrac = 1.0

    ##### Function calls

    # Create atmosphere object
    atm            = atmos(T_surf, P_surf, P_top, pl_radius, pl_mass, vol_mixing=vol_mixing, trppT=trppT)

    # Compute stellar heating
    _, atm.toa_heating = InterpolateStellarLuminosity(star_mass, time, mean_distance, atm.albedo_pl, Sfrac)

    # Do rad trans
    _, atm_moist = RadConvEqm(dirs, time, atm, standalone=True, cp_dry=False, trppD=trppD, calc_cf=calc_cf, rscatter=rscatter) 

    return [T_surf, atm_moist.LW_flux_up[0]]


if __name__=='__main__':

    print("Start")
    print(" ")

    # Set up dirs
    dirs = {
            "rad_conv": os.getenv('AEOLUS_DIR')+"/",
            "output": os.getenv('AEOLUS_DIR')+"/output/"
            }
    
    # Tidy directory
    if os.path.exists(dirs["output"]):
        shutil.rmtree(dirs["output"])
    os.mkdir(dirs["output"])

    # Setup spectral file
    print("Inserting stellar spectrum")
    StellarSpectrum.InsertStellarSpectrum(
        dirs["rad_conv"]+"/spectral_files/Oak/Oak",
        dirs["rad_conv"]+"/spectral_files/stellar_spectra/Sun_t4_4Ga_claire_12.txt",
        dirs["output"]+"runtime_spectral_file"
    )
    print(" ")

    
    # Run AEOLUS in a loop to generate runaway curve
    print("Running AEOLUS...")
    Ts_arr = []
    OLR_arr = []
    for Ts in np.linspace(200, 2200, 30):
        print("T_surf = %d" % Ts)
        out = run_once(Ts, dirs)
        Ts_arr.append(out[0])
        OLR_arr.append(out[1])
        print(" ")
    
    # Get literature data
    g2013 = np.loadtxt(dirs["rad_conv"]+"plotting_tools/comparison_data/Goldblatt13_data.txt",
                          dtype=float, skiprows=2, delimiter=',').T 
    k2013 = np.loadtxt(dirs["rad_conv"]+"plotting_tools/comparison_data/Kopparapu13_data.txt",
                          dtype=float, skiprows=2, delimiter=',').T 
    h2015 = np.loadtxt(dirs["rad_conv"]+"plotting_tools/comparison_data/Hamano15_data.txt",
                          dtype=float, skiprows=2, delimiter=',').T 

    # Setup plot
    print("Making plot")
    fig,ax = plt.subplots(1,1)

    # Plot data
    lw = 2
    ax.plot(k2013[0], k2013[1], color='tab:red',   lw=lw, label='Kopparapu+2013')
    ax.plot(g2013[0], g2013[1], color='tab:green', lw=lw, label='Goldblatt+2013')
    ax.plot(h2015[0], h2015[1], color='tab:blue',  lw=lw, label='Hamano+2015')
    ax.plot(Ts_arr, OLR_arr,    color='black',     lw=lw, label='AEOLUS')

    # Setup figure and save
    fig.legend(loc='upper center')
    ax.set_xlabel("Surface temperature [K]")
    ax.set_ylabel("OLR [W m-2]")
    fig.savefig(dirs["output"]+"runaway_demo.pdf")
    print(" ")

    # Done
    print("Done!")

