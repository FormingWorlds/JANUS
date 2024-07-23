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


DATA_DRC = files("janus") / "data" / "tests"


def get_spectrum_data(drc):
    DownloadSpectralFiles("/Oak")
    DownloadStellarSpectra()

    spec = mors.Spectrum()
    spec.LoadTSV(os.environ.get("FWL_DATA") + "/stellar_spectra/Named/sun.txt")
    socstar = os.path.join(drc, "socstar.txt")
    StellarSpectrum.PrepareStellarSpectrum(spec.wl, spec.fl, socstar)

    StellarSpectrum.InsertStellarSpectrum(
        os.environ.get("FWL_DATA") + "/spectral_files/Oak/318/Oak.sf",
        socstar,
        drc,
    )
    band_edges = ReadBandEdges(drc + "star.sf")

    return band_edges


def get_atmosphere_config(*, band_edges, cfg_name: str):
    cfg_file = DATA_DRC / cfg_name

    with open(cfg_file) as f:
        cfg = toml.load(f)

    star_time = cfg["star"]["time"]
    star_mass = cfg["star"]["star_mass"]
    distance = cfg["star"]["mean_distance"]

    vol_mixing = {"H2O": 1.0, "CO2": 0.0, "N2": 0.0}

    atm = atmos.from_file(cfg_file, band_edges, vol_mixing=vol_mixing, vol_partial={})

    mors.DownloadEvolutionTracks("/Baraffe")
    baraffe = mors.BaraffeTrack(star_mass)
    atm.instellation = baraffe.BaraffeSolarConstant(star_time, distance)

    return atm, cfg
