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

from helpers import get_spectrum_data, get_atmosphere_config


JANUS_DRC = str(files("janus")) + "/"


def test_runaway_greenhouse(tmpdir):
    out_drc = str(tmpdir) + "/"

    band_edges = get_spectrum_data(out_drc)
    atm, cfg = get_atmosphere_config(
        band_edges=band_edges, cfg_name="config_runaway.toml"
    )

    time = {"planet": cfg["planet"]["time"], "star": cfg["star"]["time"]}

    Ts = 200
    expected = 9.07314e01

    atm.setSurfaceTemperature(Ts)

    _, atm_moist = RadConvEqm(
        {"output": out_drc},
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
