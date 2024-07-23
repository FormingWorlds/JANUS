import os
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

JANUS_DRC = str(files("janus")) + "/"


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


def get_atmosphere_config(band_edges):
    cfg_file = JANUS_DRC + "data/tests/config_runaway.toml"
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


def test_runaway_greenhouse(tmpdir):
    out_drc = str(tmpdir) + "/"

    band_edges = get_spectrum_data(out_drc)
    atm, cfg = get_atmosphere_config(band_edges=band_edges)

    time = {"planet": cfg["planet"]["time"], "star": cfg["star"]["time"]}

    Ts = 200
    expected = 9.07314e01

    atm.setSurfaceTemperature(Ts)

    _, atm_moist = RadConvEqm(
        {"janus": JANUS_DRC, "output": out_drc},
        time,
        atm,
        standalone=True,
        cp_dry=False,
        trppD=False,
        rscatter=False,
    )

    ret = atm_moist.LW_flux_up[0]

    np.testing.assert_allclose(ret, expected, rtol=1e-5, atol=0)

    CleanOutputDir(os.getcwd())
    CleanOutputDir(out_drc)
