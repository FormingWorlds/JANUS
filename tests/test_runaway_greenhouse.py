#!/usr/bin/env python3

import os, shutil
import numpy as np
from matplotlib.ticker import MultipleLocator
from importlib.resources import files

from janus.modules.stellar_luminosity import InterpolateStellarLuminosity
from janus.modules.solve_pt import RadConvEqm
from janus.utils.socrates import CleanOutputDir

from janus.utils.atmosphere_column import atmos
import janus.utils.StellarSpectrum as StellarSpectrum
from janus.utils.ReadSpectralFile import ReadBandEdges


def run_once(T_surf, dirs, band_edges):

    # Planet 
    time = { "planet": 0., "star": 4.5e9 } # yr,
    star_mass     = 1.0                 # M_sun, mass of star
    mean_distance = 0.5                 # au, orbital distance
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

    # Tropopause calculation
    trppT = 0.0     # Fixed tropopause value if not calculated dynamically
    
    ##### Function calls

    # Create atmosphere object
    atm            = atmos(T_surf, P_surf, P_top, pl_radius, pl_mass, band_edges, vol_mixing=vol_mixing, trppT=trppT)

    # Compute stellar heating
    atm.instellation = InterpolateStellarLuminosity(star_mass, time, mean_distance)

    # Do rad trans
    _, atm_moist = RadConvEqm(dirs, time, atm, standalone=True, cp_dry=False, trppD=False, rscatter=rscatter) 

    return [atm_moist.LW_flux_up[0]]

def test_runaway_greenhouse():

    # Set up dirs
    if os.environ.get('RAD_DIR') == None:
        raise Exception("Socrates environment variables not set! Have you installed Socrates and sourced set_rad_env?")
    dirs = {
            "janus": str(files("janus"))+"/",
            "output": os.path.abspath(os.getcwd())+"/output/"
            }

    # Tidy directory
    if os.path.exists(dirs["output"]):
        shutil.rmtree(dirs["output"])
    os.mkdir(dirs["output"])

    # Setup spectral file
    print("Inserting stellar spectrum")
    StellarSpectrum.InsertStellarSpectrum(
        dirs["janus"]+"data/spectral_files/Oak/Oak.sf",
        dirs["janus"]+"data/spectral_files/stellar_spectra/Sun_t4_4Ga_claire_12.txt",
        dirs["output"]
    )
    band_edges = ReadBandEdges(dirs["output"]+"star.sf")

    #Get reference values
    OLR_ref = np.loadtxt(dirs["janus"]+"data/tests/data_runaway_greenhouse.csv",
                         dtype=float, skiprows=1, delimiter=',')

    #Run Janus
    Ts_arr = np.linspace(200, 2800, 20)
    for i in range(20):
      out = run_once(Ts_arr[i], dirs, band_edges)
      print("Output %s; Reference %s" % (out, OLR_ref[i][1]))
      np.testing.assert_allclose(out, OLR_ref[i][1], rtol=1e-5, atol=0)      

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])
