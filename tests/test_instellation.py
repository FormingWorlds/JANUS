"""Pipeline reference test for the MCPA_CBL surface energy balance.

Runs the full JANUS single-column radiative-convective solver (real SOCRATES
binary) at fixed instellation and pins the returned top-of-atmosphere shortwave
and longwave fluxes, net flux, surface temperature, and tropopause temperature
against reference values. This is the instellation flux anchor in the
reference_pinned inventory for the solve_pt pipeline. See docs/How-to/test.md.
The active reference values come from the SOCRATES version in helpers.SOCRATES_PINNED_VERSION.
"""

import numpy as np
import pytest

pytest.importorskip('mors')

from helpers import (  # noqa: E402
    SOCRATES_PIN_NOTE,
    get_atmosphere_config,
    get_spectrum_data,
    work_directory,
)

from janus.modules import MCPA_CBL  # noqa: E402

pytestmark = [pytest.mark.integration, pytest.mark.timeout(300)]

TEST_DATA = (
    (0.3, (2.3429661e03, 1.9457766e03, -4.4875244e01, 3.0002244e03, 4.1957014e02)),
    #    (0.483333, (9.02641e02, 8.43916e02, 7.70272e01, 2.99961e03, 3.30552e02)),
    #    (0.666667, (4.74451e02, 5.01731e02, 9.86425e01, 2.99951e03, 2.81454e02)),
    #    (0.85, (2.91857e02, 4.21028e02, 1.73069e02, 2.99913e03, 2.49260e02)),
    #    (1.03333, (1.97482e02, 4.16319e02, 2.48540e02, 2.99876e03, 2.26070e02)),
    #    (1.21667, (1.42451e02, 4.15915e02, 2.94890e02, 2.99853e03, 2.08342e02)),
    (1.4, (1.0758520e02, 4.1712155e02, 3.2571954e02, 2.9983714e03, 1.9422323e02)),
)
CASE_IDS = {0.3: 'inner_0p3au', 1.4: 'outer_1p4au'}


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
@pytest.mark.parametrize('inp,expected', TEST_DATA, ids=[CASE_IDS[x] for x, _ in TEST_DATA])
def test_instellation(tmpdir, inp, expected, validation_plot):
    """Full-column MCPA_CBL fluxes match pinned references at fixed instellation.

    For a pure water-vapour column at two orbital distances the converged
    single-column solver must reproduce the reference top-of-atmosphere
    downward shortwave flux, upward longwave flux, net flux, surface temperature, and
    tropopause temperature. Before the pin the physical sanity of the solve is
    checked: surface and tropopause temperatures are positive and every
    returned flux is finite, ruling out a diverged or NaN-producing solve that
    could otherwise coincidentally satisfy a loose tolerance.
    """
    out_drc = str(tmpdir) + '/'

    band_edges = get_spectrum_data(out_drc)
    atm, cfg = get_atmosphere_config(
        band_edges=band_edges,
        distance=inp,
        cfg_name='config_instellation.toml',
    )
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

    ret = (atm.SW_flux_down[0], atm.LW_flux_up[0], atm.net_flux[0], atm.ts, atm.trppT)
    # Physical sanity before the pin: positive temperatures, finite fluxes.
    assert atm.ts > 0.0
    assert atm.trppT > 0.0
    assert np.all(np.isfinite(ret))
    np.testing.assert_allclose(ret, expected, rtol=1e-5, atol=0, err_msg=SOCRATES_PIN_NOTE)

    def figure():
        import matplotlib.pyplot as plt

        names = ['SW down', 'LW up', 'net', 'Ts', 'T trop']
        fig, ax = plt.subplots()
        ax.axhline(1.0, ls=':', color='0.5')
        ax.plot(names, np.asarray(ret) / np.asarray(expected), 'o', label='code / stored value')
        ax.set(ylabel='code value / stored value', ylim=(0.99, 1.01))
        ax.legend(fontsize='small')
        fig.tight_layout()
        return fig

    validation_plot(figure, f'janus_instellation_{inp:g}au'.replace('.', 'p'))
