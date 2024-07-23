import os
from contextlib import contextmanager
from importlib.resources import files
from pathlib import Path

import mors
import toml
from janus.utils import (
    DownloadSpectralFiles,
    DownloadStellarSpectra,
    ReadBandEdges,
    StellarSpectrum,
    atmos,
)

DATA_DRC = files('janus') / 'data' / 'tests'
FWL_DATA = os.environ.get('FWL_DATA')


@contextmanager
def work_directory(path: Path | str):
    """Changes working directory and returns to previous on exit.

    Parameters
    ----------
    path : Path | str
        Temporarily change to this directory.
    """
    prev_cwd = Path.cwd().resolve()
    try:
        os.chdir(path)
        yield
    finally:  # In any case, no matter what happens, go back eventually
        os.chdir(prev_cwd)


def get_spectrum_data(drc):
    DownloadSpectralFiles('/Oak')
    DownloadStellarSpectra()

    spec = mors.Spectrum()
    spec.LoadTSV(FWL_DATA + '/stellar_spectra/Named/sun.txt')
    socstar = os.path.join(drc, 'socstar.txt')
    StellarSpectrum.PrepareStellarSpectrum(spec.wl, spec.fl, socstar)

    StellarSpectrum.InsertStellarSpectrum(
        FWL_DATA + '/spectral_files/Oak/318/Oak.sf',
        socstar,
        drc,
    )
    band_edges = ReadBandEdges(drc + 'star.sf')

    return band_edges


def get_atmosphere_config(*, band_edges, cfg_name: str):
    cfg_file = DATA_DRC / cfg_name

    with open(cfg_file) as f:
        cfg = toml.load(f)

    star_time = cfg['star']['time']
    star_mass = cfg['star']['star_mass']
    distance = cfg['star']['mean_distance']

    vol_mixing = {'H2O': 1.0, 'CO2': 0.0, 'N2': 0.0}

    atm = atmos.from_file(cfg_file, band_edges, vol_mixing=vol_mixing, vol_partial={})

    mors.DownloadEvolutionTracks('/Baraffe')
    baraffe = mors.BaraffeTrack(star_mass)
    atm.instellation = baraffe.BaraffeSolarConstant(star_time, distance)

    return atm, cfg
