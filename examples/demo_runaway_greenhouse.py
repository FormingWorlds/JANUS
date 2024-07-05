#!/usr/bin/env python3

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams.update({'font.size': 12})

import matplotlib.pyplot as plt
import os, shutil, toml
import numpy as np
from matplotlib.ticker import MultipleLocator
from importlib.resources import files

from janus.modules.stellar_luminosity import InterpolateStellarLuminosity
from janus.modules.solve_pt import RadConvEqm
from janus.utils.socrates import CleanOutputDir

from janus.utils.atmosphere_column import atmos
import janus.utils.StellarSpectrum as StellarSpectrum
from janus.utils.ReadSpectralFile import ReadBandEdges
from janus.utils.data import DownloadSpectralFiles

if __name__=='__main__':

    print("Start")
    print(" ")

    # Set up dirs
    if os.environ.get('RAD_DIR') == None:
        raise Exception("Socrates environment variables not set! Have you installed Socrates and sourced set_rad_env?")
    if os.environ.get('FWL_DATA') == None:
        raise Exception("The FWL_DATA environment variable where spectral data will be downloaded needs to be set up!")
    dirs = {
            "janus": str(files("janus"))+"/",
            "output": os.path.abspath(os.getcwd())+"/output/"
            }

    # Tidy directory
    if os.path.exists(dirs["output"]):
        shutil.rmtree(dirs["output"])
    os.mkdir(dirs["output"])

    #Download required spectral files
    DownloadSpectralFiles("/Oak")
    DownloadSpectralFiles("/stellar_spectra")

    # Setup spectral file
    print("Inserting stellar spectrum")
    StellarSpectrum.InsertStellarSpectrum(
        os.environ.get('FWL_DATA')+"/spectral_files/Oak/318/Oak.sf",
        os.environ.get('FWL_DATA')+"/spectral_files/stellar_spectra/Sun_t4_4Ga_claire_12.txt",
        dirs["output"]
    )
    print(" ")
    band_edges = ReadBandEdges(dirs["output"]+"star.sf")

    # Open config file
    cfg_file =  dirs["janus"]+"data/tests/config_runaway.toml"
    with open(cfg_file, 'r') as f:
          cfg = toml.load(f)

    # Planet
    time = { "planet": cfg['planet']['time'], "star": cfg['star']['time']}
    star_mass = cfg['star']['star_mass']
    mean_distance = cfg['star']['mean_distance']

    vol_mixing = {
                    "H2O" : 1.0,
                    "CO2" : 0.0,
                    "N2"  : 0.0
                }
    # Create atmosphere object
    atm = atmos.from_file(cfg_file, band_edges, vol_mixing=vol_mixing, vol_partial={})

    # Compute stellar heating
    atm.instellation = InterpolateStellarLuminosity(star_mass, time, mean_distance)

    #Run Janus
 
    # Run JANUS in a loop to generate runaway curve
    print("Running JANUS...")
    Ts_arr = np.linspace(200, 2800, 20)
    OLR_arr = []
    for i in range(20):
      print("T_surf = %d K" % Ts_arr[i])
      atmos.setSurfaceTemperature(atm, Ts_arr[i])

      _, atm_moist = RadConvEqm(dirs, time, atm, standalone=True, cp_dry=False, trppD=False, rscatter=False)

      OLR_arr.append(atm_moist.LW_flux_up[0])
      print(" ")

    OLR_arr = np.array(OLR_arr)
    Ts_arr  = np.array(Ts_arr)

    # Get literature data
    g2013 = np.loadtxt(dirs["janus"]+"data/comparison_data/Goldblatt13_data.txt",
                          dtype=float, skiprows=2, delimiter=',').T 
    k2013 = np.loadtxt(dirs["janus"]+"data/comparison_data/Kopparapu13_data.txt",
                          dtype=float, skiprows=2, delimiter=',').T 
    h2015 = np.loadtxt(dirs["janus"]+"data/comparison_data/Hamano15_data.txt",
                          dtype=float, skiprows=2, delimiter=',').T 
    s2023 = np.loadtxt(dirs["janus"]+"data/comparison_data/Selsis23_convective.txt",
                          dtype=float, skiprows=2, delimiter=',').T 

    # Setup plot
    print("Making plot")
    fig,ax = plt.subplots(1,1, figsize=(7,4))

    # Plot data
    lw = 2
    ax.plot(k2013[0], k2013[1], color='tab:red',   lw=lw, label='Kopparapu+2013')
    ax.plot(g2013[0], g2013[1], color='tab:green', lw=lw, label='Goldblatt+2013')
    ax.plot(h2015[0], h2015[1], color='tab:blue',  lw=lw, label='Hamano+2015')
    ax.plot(s2023[0], s2023[1], color='tab:orange',lw=lw, label='Selsis+2023')
    ax.plot(Ts_arr,   OLR_arr,  color='black',     lw=lw, label='JANUS')

    # Setup figure and save
    ax.legend(loc='upper left')

    ax.set_xlabel("Surface temperature [K]")
    ax.xaxis.set_minor_locator(MultipleLocator(100.0))
    ax.set_xlim(200.0,  2700.0)  

    ax.set_ylabel("OLR [W m$^{-2}$]")
    ax.set_ylim(np.amin(OLR_arr) - 10.0, 500.0)
    ax.yaxis.set_minor_locator(MultipleLocator(25.0))  

    fig.savefig(dirs["output"]+"runaway_demo.pdf", bbox_inches='tight')
    fig.savefig(dirs["output"]+"runaway_demo.png", bbox_inches='tight', dpi=190)
    print(" ")

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])

    # Done
    print("Done!")

