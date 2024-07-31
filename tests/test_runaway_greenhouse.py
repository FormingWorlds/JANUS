import numpy as np
from helpers import get_atmosphere_config, get_spectrum_data, work_directory
from janus.modules import RadConvEqm
import pytest


TEST_DATA = (
    (200.00000000,9.07314e+01),
    (1705.26315789,2.78884e+02),
    (2800.00000000,6.58174e+03),
)


@pytest.mark.parametrize("inp,expected", TEST_DATA)
def test_runaway_greenhouse(tmpdir, inp, expected):
    out_drc = str(tmpdir) + '/'

    band_edges = get_spectrum_data(out_drc)
    atm, cfg = get_atmosphere_config(band_edges=band_edges, distance=0.3, cfg_name='config_runaway.toml')

    time = {'planet': cfg['planet']['time'], 'star': cfg['star']['time']}

    atm.setSurfaceTemperature(inp)

    with work_directory(tmpdir):
        _, atm_moist = RadConvEqm(
            {'output': out_drc},
            time,
            atm,
            standalone=True,
            cp_dry=False,
            trppD=False,
            rscatter=False,
        )

    ret = atm_moist.LW_flux_up[0]
    np.testing.assert_allclose(ret, expected, rtol=1e-5, atol=0)
