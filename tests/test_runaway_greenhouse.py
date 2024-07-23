import numpy as np
from helpers import get_atmosphere_config, get_spectrum_data, work_directory
from janus.modules import RadConvEqm


def test_runaway_greenhouse(tmpdir):
    out_drc = str(tmpdir) + '/'

    band_edges = get_spectrum_data(out_drc)
    atm, cfg = get_atmosphere_config(band_edges=band_edges, cfg_name='config_runaway.toml')

    time = {'planet': cfg['planet']['time'], 'star': cfg['star']['time']}

    atm.setSurfaceTemperature(200)

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
    np.testing.assert_allclose(ret, 9.07314e01, rtol=1e-5, atol=0)
