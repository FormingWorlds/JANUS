import pytest
from helpers import get_atmosphere_config, get_spectrum_data, work_directory
from janus.modules import MCPA_CBL
from numpy.testing import assert_allclose

TEST_DATA = (
    (0.3, (2.34297e03, 1.94397e03, -4.67185e01, 3.00023e03, 4.19568e02)),
#    (0.483333, (9.02641e02, 8.43916e02, 7.70272e01, 2.99961e03, 3.30552e02)),
#    (0.666667, (4.74451e02, 5.01731e02, 9.86425e01, 2.99951e03, 2.81454e02)),
#    (0.85, (2.91857e02, 4.21028e02, 1.73069e02, 2.99913e03, 2.49260e02)),
#    (1.03333, (1.97482e02, 4.16319e02, 2.48540e02, 2.99876e03, 2.26070e02)),
#    (1.21667, (1.42451e02, 4.15915e02, 2.94890e02, 2.99853e03, 2.08342e02)),
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
