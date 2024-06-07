#!/usr/bin/env python3

from importlib.resources import files
import os, shutil
import numpy as np

from janus.modules.stellar_luminosity import InterpolateStellarLuminosity
from janus.modules.solve_pt import *
from janus.utils.socrates import CleanOutputDir
from janus.utils.atmosphere_column import atmos
import janus.utils.StellarSpectrum as StellarSpectrum
import janus.utils.phys as phys
from janus.utils.ReadSpectralFile import ReadBandEdges

def run_once(sep, dirs, T_magma, P_surf, skin_d, band_edges):

    # Planet 
    time = { "planet": 0., "star": 100e+6 } # yr,
    star_mass     = 1.0                 # M_sun, mass of star
    pl_radius     = 6.371e6             # m, planet radius
    pl_mass       = 5.972e24            # kg, planet mass

    # Boundary conditions for pressure & temperature
    P_top         = 1.0                  # Pa

    # Define volatiles by mole fractions
    vol_mixing = {
                    "H2O" : 1.0,
                    "CO2" : 0.0,
                    "N2"  : 0.0
                }
    
    rscatter = True
    A_B = 0.1  # bond albedo
    A_S = 0.1
    inst_sf = 3.0/8.0
    
    ##### Function calls

    S_0 = InterpolateStellarLuminosity(star_mass, time, sep) 

    zenith_angle = 48.19 # cronin+14 (also for scaling by a factor of 3/8 ^^)

    T_eqm = (S_0 * inst_sf * (1.0 - A_B) /phys.sigma)**(1.0/4.0)
    T_trpp = T_eqm * (0.5**0.25)  # radiative skin temperature

    # Create atmosphere object
    atm = atmos(T_magma,  P_surf * 1e5, P_top, pl_radius, pl_mass, band_edges,
                vol_mixing=vol_mixing, trppT=T_trpp, req_levels=150)
    atm.albedo_pl = A_B
    atm.albedo_s  = A_S
    atm.inst_sf = inst_sf
    atm.zenith_angle = zenith_angle
    atm.instellation = S_0
    atm.skin_d = skin_d
    atm.tmp_magma = T_magma

    # Do rad trans
    atm = MCPA_CBL(dirs, atm, False, rscatter, T_surf_max=9.0e99, T_surf_guess = T_trpp+100)

    return [atm.SW_flux_down[0], atm.LW_flux_up[0], atm.net_flux[0], atm.ts, T_trpp]

def test_instellation():

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

    # Set parameters
    P_surf  = 280.0   # surface pressure [bar]
    T_magma = 3000.0  # magma temperature [K]
    skin_d  = 1e-2  # conductive skin thickness [m]

    #Get reference data
    ref = np.loadtxt(dirs["janus"]+"data/tests/data_instellation.csv",
                     dtype=float, skiprows=1, delimiter=',')
     
    r_arr = np.linspace(0.3, 1.4, 7) # orbital distance range [AU]
    for i in range(7):
      print("Orbital separation = %.2f AU" % r_arr[i])
      out = run_once(r_arr[i], dirs, T_magma, P_surf, skin_d, band_edges)
      print(out)
      print(ref[i][1:6])
      np.testing.assert_allclose(out, ref[i][1:6], rtol=1e-5, atol=0)

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])
