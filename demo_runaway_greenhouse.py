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
                    "CO2"  : 0.0,
                    "H2O"  : 1.0,
                    "N2"   : 0.0,
                    "H2"   : 0.0, 
                    "NH3"  : 0.0,
                    "CH4"  : 0.0, 
                    "O2"   : 0.0, 
                    "CO"   : 0.0
                }


    
    # Rayleigh scattering on/off
    rscatter = False

    # Compute contribution function
    calc_cf = False

    # Pure steam convective adjustment
    pure_steam_adj = False

    # Tropopause calculation
    trppD = False   # Calculate dynamically?
    trppT = 0.0     # Fixed tropopause value if not calculated dynamically
    
    # Surface temperature time-stepping
    surf_dt = False
    cp_dry = False
    # Options activated by surf_dt
    cp_surf = 1e5         # Heat capacity of the ground [J.kg^-1.K^-1]
    mix_coeff_atmos = 1e6 # mixing coefficient of the atmosphere [s]
    mix_coeff_surf  = 1e6 # mixing coefficient at the surface [s]

    # Instellation scaling | 1.0 == no scaling
    Sfrac = 1.0

    ##### Function calls

    
    # Create atmosphere object
    atm            = atmos(T_surf, P_surf, P_top, pl_radius, pl_mass, vol_mixing=vol_mixing, trppT=trppT)

    # Compute stellar heating
    _, atm.toa_heating = InterpolateStellarLuminosity(star_mass, time, mean_distance, atm.albedo_pl, Sfrac)

    # Do rad trans
    _, atm_moist = RadConvEqm(dirs, time, atm, standalone=True, cp_dry=cp_dry, trppD=trppD, calc_cf=calc_cf, rscatter=rscatter, pure_steam_adj=pure_steam_adj, surf_dt=surf_dt, cp_surf=cp_surf, mix_coeff_atmos=mix_coeff_atmos, mix_coeff_surf=mix_coeff_surf) 


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
    for Ts in np.linspace(200, 2200, 40):
        print("T_surf = %d" % Ts)
        out = run_once(Ts, dirs)
        Ts_arr.append(out[0])
        OLR_arr.append(out[1])
        print(" ")
    
    # Get literature data
    k2013 = np.loadtxt(dirs["rad_conv"]+"tools/kopparapu+2013_runaway_curve.csv",
                          dtype=float, skiprows=4, delimiter=',').T 
    g2013 = np.loadtxt(dirs["rad_conv"]+"tools/goldblatt+2013_runaway_curve.csv",
                          dtype=float, skiprows=4, delimiter=',').T 

    # Setup plot
    print("Making plot")
    fig,ax = plt.subplots(1,1)

    # Plot literature data
    ax.plot(k2013[0], k2013[1], color='tab:red',   label='Kopparapu+2013')
    ax.plot(g2013[0], g2013[1], color='tab:green', label='Goldblatt+2013')

    # Plot our data
    ax.plot(Ts_arr, OLR_arr, color='black', label='AEOLUS')

    # Setup figure and save
    fig.legend()
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("OLR [W m-2]")
    fig.savefig(dirs["output"]+"runaway_demo.pdf")
    print(" ")

    # Done
    print("Done!")

