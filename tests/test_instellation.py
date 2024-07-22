from importlib.resources import files
import os
import shutil
import toml
import numpy as np

from janus.modules import MCPA_CBL
from janus.utils import (
    atmos,
    CleanOutputDir,
    DownloadSpectralFiles,
    DownloadStellarSpectra,
    ReadBandEdges,
    StellarSpectrum,
)
import mors


def test_instellation(tmpdir):
    
    janus = str(files("janus")) + "/"
    output = str(tmpdir) + "/"

    # Download required spectral files
    DownloadSpectralFiles("/Oak")
    DownloadStellarSpectra()

    # Read spectrum
    spec = mors.Spectrum()
    spec.LoadTSV(os.environ.get("FWL_DATA") + "/stellar_spectra/Named/sun.txt")

    # Convert to SOCRATES format
    socstar = os.path.join(output, "socstar.txt")
    StellarSpectrum.PrepareStellarSpectrum(spec.wl, spec.fl, socstar)

    # Setup spectral file
    print("Inserting stellar spectrum")
    StellarSpectrum.InsertStellarSpectrum(
        os.environ.get("FWL_DATA") + "/spectral_files/Oak/318/Oak.sf",
        socstar,
        output,
    )
    band_edges = ReadBandEdges(output + "star.sf")

    # Open config file
    cfg_file = janus + "data/tests/config_instellation.toml"
    with open(cfg_file) as f:
        cfg = toml.load(f)

    # Star luminosity
    time = {"planet": cfg["planet"]["time"], "star": cfg["star"]["time"]}
    star_mass = cfg["star"]["star_mass"]
    mors.DownloadEvolutionTracks("/Baraffe")
    baraffe = mors.BaraffeTrack(star_mass)

    # Define volatiles by mole fractions
    vol_mixing = {"H2O": 1.0, "CO2": 0.0, "N2": 0.0}

    # Get reference data
    ref = np.loadtxt(
        janus + "data/tests/data_instellation.csv",
        dtype=float,
        skiprows=1,
        delimiter=",",
    )

    # Create atmosphere object
    atm = atmos.from_file(cfg_file, band_edges, vol_mixing=vol_mixing, vol_partial={})

    for r, expected in (
        (0.3, ref[0]),
        (1.4, ref[6]),
    ):
        atm.instellation = baraffe.BaraffeSolarConstant(time["star"], r)
        atm.setTropopauseTemperature()

        atm = MCPA_CBL(
            {'janus':janus, 'output':output},
            atm,
            False,
            rscatter=True,
            T_surf_max=9.0e99,
            T_surf_guess=atm.trppT + 100,
        )

        np.testing.assert_allclose(atm.SW_flux_down[0], expected[1], rtol=1e-5, atol=0)
        np.testing.assert_allclose(atm.LW_flux_up[0], expected[2], rtol=1e-5, atol=0)
        np.testing.assert_allclose(atm.net_flux[0], expected[3], rtol=1e-5, atol=0)
        np.testing.assert_allclose(atm.ts, expected[4], rtol=1e-5, atol=0)
        np.testing.assert_allclose(atm.trppT, expected[5], rtol=1e-5, atol=0)

    # Tidy
    CleanOutputDir(os.getcwd())
    CleanOutputDir(output)
