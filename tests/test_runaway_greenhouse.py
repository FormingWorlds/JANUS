#!/usr/bin/env python3

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

def test_runaway_greenhouse():

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

    #Get reference values
    OLR_ref = np.loadtxt(dirs["janus"]+"data/tests/data_runaway_greenhouse.csv",
                         dtype=float, skiprows=1, delimiter=',')

    # Create atmosphere object
    atm = atmos.from_file(cfg_file, band_edges, vol_mixing=vol_mixing, vol_partial={})

    # Compute stellar heating
    atm.instellation = InterpolateStellarLuminosity(star_mass, time, mean_distance)

    #Run Janus
    Ts_arr = np.linspace(200, 2800, 20)
    for i in range(20):

      atmos.setSurfaceTemperature(atm, Ts_arr[i])

      _, atm_moist = RadConvEqm(dirs, time, atm, standalone=True, cp_dry=False, trppD=False, rscatter=False)

      out = [atm_moist.LW_flux_up[0]]
      print("Output %s; Reference %s" % (out, OLR_ref[i][1]))
      np.testing.assert_allclose(out, OLR_ref[i][1], rtol=1e-5, atol=0)      

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])
