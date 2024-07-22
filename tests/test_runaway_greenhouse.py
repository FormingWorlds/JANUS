import os
import shutil
import toml
import numpy as np
from importlib.resources import files

from janus.modules import RadConvEqm
from janus.utils import (
    atmos,
    CleanOutputDir,
    DownloadSpectralFiles,
    DownloadStellarSpectra,
    ReadBandEdges,
    StellarSpectrum,
)
import mors


def test_runaway_greenhouse():
    dirs = {
        "janus": str(files("janus")) + "/",
        "output": os.path.abspath(os.getcwd()) + "/output/",
    }

    # Tidy directory
    if os.path.exists(dirs["output"]):
        shutil.rmtree(dirs["output"])
    os.mkdir(dirs["output"])

    # Download required spectral files
    DownloadSpectralFiles("/Oak")
    DownloadStellarSpectra()

    # Read spectrum
    spec = mors.Spectrum()
    spec.LoadTSV(os.environ.get("FWL_DATA") + "/stellar_spectra/Named/sun.txt")

    # Convert to SOCRATES format
    socstar = os.path.join(dirs["output"], "socstar.txt")
    StellarSpectrum.PrepareStellarSpectrum(spec.wl, spec.fl, socstar)

    # Setup spectral file
    print("Inserting stellar spectrum")
    StellarSpectrum.InsertStellarSpectrum(
        os.environ.get("FWL_DATA") + "/spectral_files/Oak/318/Oak.sf",
        socstar,
        dirs["output"],
    )
    band_edges = ReadBandEdges(dirs["output"] + "star.sf")

    # Open config file
    cfg_file = dirs["janus"] + "data/tests/config_runaway.toml"
    with open(cfg_file) as f:
        cfg = toml.load(f)

    # Planet
    time = {"planet": cfg["planet"]["time"], "star": cfg["star"]["time"]}
    star_mass = cfg["star"]["star_mass"]
    mean_distance = cfg["star"]["mean_distance"]

    vol_mixing = {"H2O": 1.0, "CO2": 0.0, "N2": 0.0}

    # Get reference values
    OLR_ref = np.loadtxt(
        dirs["janus"] + "data/tests/data_runaway_greenhouse.csv",
        dtype=float,
        skiprows=1,
        delimiter=",",
    )

    # Create atmosphere object
    atm = atmos.from_file(cfg_file, band_edges, vol_mixing=vol_mixing, vol_partial={})

    # Compute stellar heating from Baraffe evolution tracks
    mors.DownloadEvolutionTracks("/Baraffe")
    baraffe = mors.BaraffeTrack(star_mass)
    atm.instellation = baraffe.BaraffeSolarConstant(time["star"], mean_distance)

    for Ts, expected in (
        (200, OLR_ref[0]),
        (2800, OLR_ref[19]),
    ):
        atmos.setSurfaceTemperature(atm, Ts)

        _, atm_moist = RadConvEqm(
            dirs, time, atm, standalone=True, cp_dry=False, trppD=False, rscatter=False
        )

        out = atm_moist.LW_flux_up[0]

        np.testing.assert_allclose(out, expected[1], rtol=1e-5, atol=0)

    CleanOutputDir(os.getcwd())
    CleanOutputDir(dirs["output"])
