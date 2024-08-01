import numpy as np
from helpers import get_atmosphere_config, get_spectrum_data, work_directory
from janus.modules import RadConvEqm
import pytest


TEST_DATA = (
    (200.00000000,9.07314e+01),
#    (336.84210526,2.77416e+02),
#    (473.68421053,2.77349e+02),
#    (610.52631579,2.77138e+02),
#    (747.36842105,2.77217e+02),
#    (884.21052632,2.77049e+02),
#    (1021.05263158,2.77241e+02),
#    (1157.89473684,2.76997e+02),
#    (1294.73684211,2.77013e+02),
#    (1431.57894737,2.77166e+02),
#    (1568.42105263,2.77504e+02),
    (1705.26315789,2.78884e+02),
#    (1842.10526316,2.85329e+02),
#    (1978.94736842,3.07788e+02),
#    (2115.78947368,3.74184e+02),
#    (2252.63157895,5.47014e+02),
#    (2389.47368421,9.50390e+02),
#    (2526.31578947,1.80287e+03),
#    (2663.15789474,3.49537e+03),
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
