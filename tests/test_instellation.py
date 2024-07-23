import numpy as np
from helpers import get_atmosphere_config, get_spectrum_data, work_directory
from janus.modules import MCPA_CBL


def test_instellation(tmpdir):
    out_drc = str(tmpdir) + '/'

    band_edges = get_spectrum_data(out_drc)
    atm, cfg = get_atmosphere_config(band_edges=band_edges, cfg_name='config_instellation.toml')
    atm.setTropopauseTemperature()

    with work_directory(tmpdir):
        atm = MCPA_CBL(
            {'output': out_drc},
            atm,
            False,
            rscatter=True,
            T_surf_max=9.0e99,
            T_surf_guess=atm.trppT + 100,
        )

    for ret, expected in (
        (atm.SW_flux_down[0], 2.34297e03),
        (atm.LW_flux_up[0], 1.94397e03),
        (atm.net_flux[0], -4.67185e01),
        (atm.ts, 3.00023e03),
        (atm.trppT, 4.19568e02),
    ):
        np.testing.assert_allclose(ret, expected, rtol=1e-5, atol=0)
