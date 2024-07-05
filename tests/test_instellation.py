#!/usr/bin/env python3

from importlib.resources import files
import os, shutil, toml
import numpy as np

from janus.modules.stellar_luminosity import InterpolateStellarLuminosity
from janus.modules.solve_pt import *
from janus.utils.socrates import CleanOutputDir
from janus.utils.atmosphere_column import atmos
import janus.utils.StellarSpectrum as StellarSpectrum
import janus.utils.phys as phys
from janus.utils.ReadSpectralFile import ReadBandEdges
from janus.utils.data import DownloadSpectralFiles

def test_instellation():

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
    band_edges = ReadBandEdges(dirs["output"]+"star.sf")

    # Open config file
    cfg_file =  dirs["janus"]+"data/tests/config_instellation.toml"
    with open(cfg_file, 'r') as f:
          cfg = toml.load(f)

    # Planet
    time = { "planet": cfg['planet']['time'], "star": cfg['star']['time']}
    star_mass = cfg['star']['star_mass']

    # Define volatiles by mole fractions
    vol_mixing = {
                    "H2O" : 1.0,
                    "CO2" : 0.0,
                    "N2"  : 0.0
                }

    #Get reference data
    ref = np.loadtxt(dirs["janus"]+"data/tests/data_instellation.csv",
                     dtype=float, skiprows=1, delimiter=',')

    # Create atmosphere object
    atm = atmos.from_file(cfg_file, band_edges, vol_mixing=vol_mixing, vol_partial={})
     
    r_arr = np.linspace(0.3, 1.4, 7) # orbital distance range [AU]
    for i in range(7):
      print("Orbital separation = %.2f AU" % r_arr[i])

      atm.instellation = InterpolateStellarLuminosity(star_mass, time, r_arr[i])
      atmos.setTropopauseTemperature(atm)

      atm = MCPA_CBL(dirs, atm, False, rscatter = True, T_surf_max=9.0e99, T_surf_guess = atm.trppT+100)

      out = [atm.SW_flux_down[0], atm.LW_flux_up[0], atm.net_flux[0], atm.ts, atm.trppT]
      print(out)
      print(ref[i][1:6])
      np.testing.assert_allclose(out, ref[i][1:6], rtol=1e-5, atol=0)

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])
