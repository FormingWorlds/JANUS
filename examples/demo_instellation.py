#!/usr/bin/env python3

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams.update({'font.size': 12})
mpl.rcParams['axes.formatter.useoffset'] = False
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from importlib.resources import files

import os, shutil, toml
import numpy as np

from janus.utils.logs import setup_logger
log = setup_logger()

from janus.modules import MCPA_CBL
from janus.utils import atmos, CleanOutputDir, DownloadSpectralFiles, DownloadStellarSpectra, ReadBandEdges, StellarSpectrum

import mors

from janus.utils.data import FWL_DATA_DIR

if __name__=='__main__':

    log.info("Start")

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

    #Download required spectral files
    DownloadSpectralFiles("Oak")
    DownloadStellarSpectra()

    # Read spectrum
    spec = mors.Spectrum()
    spec.LoadTSV(str(FWL_DATA_DIR / 'stellar_spectra' / 'Named' / 'sun.txt'))

    # Convert to SOCRATES format
    socstar = os.path.join(dirs["output"], "socstar.txt")
    StellarSpectrum.PrepareStellarSpectrum(spec.wl, spec.fl, socstar)


    # Setup spectral file
    log.info("Inserting stellar spectrum")
    StellarSpectrum.InsertStellarSpectrum(
        str(FWL_DATA_DIR / 'spectral_files'/'Oak'/'318'/'Oak.sf'),
        socstar,
        dirs["output"]
    )
    band_edges = ReadBandEdges(dirs["output"]+"star.sf")

    # Open config file
    cfg_file =  dirs["janus"]+"data/tests/config_instellation.toml"
    with open(cfg_file, 'r') as f:
          cfg = toml.load(f)

    # Star luminosity
    time = { "planet": cfg['planet']['time'], "star": cfg['star']['time']}
    star_mass = cfg['star']['star_mass']
    mors.DownloadEvolutionTracks("Baraffe")
    baraffe = mors.BaraffeTrack(star_mass)

    # Define volatiles by mole fractions
    vol_mixing = {
                    "H2O" : 1.0,
                    "CO2" : 0.0,
                    "N2"  : 0.0
                }


    # PARAMETERS
    legend  = True    # make legend?
    dx_tick = 0.1     # x-tick spacing (set to 0 for automatic)
    # /PARAMETERS

    # Create atmosphere object
    atm = atmos.from_file(cfg_file, band_edges, vol_mixing=vol_mixing, vol_partial={})
    T_magma = atm.tmp_magma #get default value 3000 K but this could be a config option

    r_arr = np.linspace(0.3, 1.4, 7) # orbital distance range [AU]

    asf_arr  = []   # ASF
    OLR_arr = []    # OLR
    net_arr = []    # net flux at TOA
    ts_arr  = []    # surface temperature
    tr_arr  = []    # tropopause temperature
    for i in range(7):
      log.info("Orbital separation = %.2f AU" % r_arr[i])

      atm.instellation = baraffe.BaraffeSolarConstant(time['star'], r_arr[i])
      atmos.setTropopauseTemperature(atm)

      atm = MCPA_CBL(dirs, atm, False, rscatter = True, T_surf_max=9.0e99, T_surf_guess = atm.trppT+100)

      asf_arr.append(atm.SW_flux_down[0] * -1.0)  # directed downwards => minus sign
      OLR_arr.append(atm.LW_flux_up[0])
      net_arr.append(atm.net_flux[0])
      ts_arr.append(atm.ts)
      tr_arr.append(atm.trppT)

      # Plot case
      plt.ioff()
      fig,ax = plt.subplots(1,1)
      ax.plot(atm.tmpl, atm.pl, color='black', lw=2)
      ax.set_yscale("log"); ax.invert_yaxis()
      ax.set_ylabel("Pressure [Pa]")
      ax.set_xlabel("Temperature [K]")
      ax.set_title("a = %.2f AU" %r_arr[i])
      fig.savefig(dirs["output"]+"/profile%.2fAU.jpg"%r_arr[i],bbox_inches='tight', dpi=100)
      plt.close()

      # Save netcdf
      atm.write_ncdf(dirs["output"]+"/profile%.2fAU.nc"%r_arr[i])

      log.info(" ")

    save_arr = [r_arr, asf_arr, OLR_arr, net_arr, ts_arr, tr_arr]
    np.savetxt(dirs["output"]+"/data_%dK.csv"%T_magma,
               np.array(save_arr).T, fmt="%.5e", delimiter=",",
               header="r [AU],  S_0 [W m-2], OLR [W m-2], net [W m-2], ts [K], tr[K] ")
    
    log.info("Making plots")

    plt.ioff()
    fig,(ax1, ax2) = plt.subplots(2,1, figsize=(6,6.4))

    ax1.set_title("%d K" % T_magma)

    ax1.text(r_arr[0], +1.0, "Cooling", verticalalignment='center', color='seagreen', weight='bold', zorder=8).set_bbox(dict(facecolor='white', alpha=0.7, linewidth=0))
    ax1.text(r_arr[0], -1.0, "Warming", verticalalignment='center', color='seagreen', weight='bold', zorder=8).set_bbox(dict(facecolor='white', alpha=0.7, linewidth=0))
    ax1.axhline(y=0, color='black', lw=0.7, zorder=0)  # zero flux

    ax1.plot(r_arr, asf_arr, color='royalblue', lw=2.5, label='ASF')    # absorbed stellar flux (sw down)
    ax1.plot(r_arr, OLR_arr, color='crimson',   lw=2.5, label='OLR')    # outgoing longwave radiation (lw up)
    ax1.plot(r_arr, net_arr, color='seagreen' , lw=2.5, label='Net')    # net upward-directed flux

    ax1.legend(loc='center right', framealpha=1.0)
    ax1.set_yscale("symlog")
    ax1.set_ylabel("Upward flux [W m$^{-2}$]")
    ax1.set_xticklabels([])
    if dx_tick > 1.0e-10:
        ax1.xaxis.set_minor_locator(MultipleLocator(dx_tick))

    yticks = ax1.get_yticks()
    yticks = np.delete(yticks, np.argwhere( yticks == 0.0 ))
    ax1.set_yticks(yticks)

    # Planet data
    x_planets = {
        "Earth": 1.0,
        "Venus": 0.723,
        "Mercury": 0.387
    }
    c_planets = {
        "Earth": 'deepskyblue',
        "Venus": 'goldenrod',
        "Mercury": 'violet'
    }
    for p in x_planets.keys():
        for ax in (ax1,ax2):
            ax.axvline(x=x_planets[p], color=c_planets[p], lw=3, linestyle='dashed', label=p , zorder=1)

    arr_magma = np.ones(len(r_arr))*T_magma
    ax2.plot(r_arr, np.zeros(len(r_arr)), zorder=6, color='silver', lw=2.5, label=r"$\tilde{T}_s$") # magma temperature
    ax2.plot(r_arr, ts_arr - arr_magma,   zorder=6, color='black',  lw=2.5, label=r"$T_s$") # surface solution temperature

    ax2.legend(loc='center right', framealpha=1.0)
    ax2.set_ylabel(r"$T - \tilde{T_s}$ [K]")
    if dx_tick > 1.0e-10:
        ax2.xaxis.set_minor_locator(MultipleLocator(dx_tick))

    ax2.set_xlabel("Orbital separation [AU]")
    fig.subplots_adjust(hspace=0.08)
    fig.savefig(dirs["output"]+"inst_%dK.pdf"%T_magma, bbox_inches='tight')


    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs['output'])

    # Done
    log.info("Done!")

