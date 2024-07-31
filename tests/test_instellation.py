from numpy.testing import assert_allclose
from helpers import get_atmosphere_config, get_spectrum_data, work_directory
from janus.modules import MCPA_CBL
import pytest


TEST_DATA = (
    (0.3, (2.34297e03, 1.94397e03, -4.67185e01, 3.00023e03, 4.19568e02)),
    (1.4, (1.07585e02, 4.15768e02, 3.24365e02, 2.99838e03, 1.94222e02)),
)


@pytest.mark.parametrize("inp,expected", TEST_DATA)
def test_instellation(tmpdir, inp, expected):
    out_drc = str(tmpdir) + "/"

    band_edges = get_spectrum_data(out_drc)
    atm, cfg = get_atmosphere_config(
        band_edges=band_edges,
        distance=inp,
        cfg_name="config_instellation.toml",
    )
    atm.setTropopauseTemperature()

    with work_directory(tmpdir):
        atm = MCPA_CBL(
            {"output": out_drc},
            atm,
            False,
            rscatter=True,
            T_surf_max=9.0e99,
            T_surf_guess=atm.trppT + 100,
        )

    ret = (atm.SW_flux_down[0], atm.LW_flux_up[0], atm.net_flux[0], atm.ts, atm.trppT)

    assert_allclose(ret, expected, rtol=1e-5, atol=0)
