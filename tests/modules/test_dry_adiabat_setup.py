"""Tests for src/janus/modules/dry_adiabat_setup.py.

dry_adiabat_setup.py provides ``dry_adiabat_atm``, which fills the staggered
temperature profile of a dry (Poisson) adiabat T = ts*(p/ps)**Rcp with
Rcp = R_universal/cp_mix. This file exercises:

* Potential-temperature conservation: theta = T*(ps/p)**Rcp is constant and
  equal to the surface temperature ts along the whole column.
* The Poisson slope: temperature rises monotonically toward the surface, the
  surface boundary condition T(ps)=ts holds, and Rcp lands in the physical
  R/cp band for common gases.

Invariant families exercised: conservation (constant potential temperature),
positivity/boundedness, monotonicity, and pinned-value-with-discrimination-
guard. See docs/How-to/test.md.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from janus.modules.dry_adiabat_setup import dry_adiabat_atm

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def _make_atm(vol_list, ts, nlev=20, ps=1.0e5):
    """Build a minimal atmosphere namespace for the dry-adiabat routine.

    pl is ascending (top of atmosphere -> surface) so pl[-1]=ps is the surface
    node; p sits at the geometric cell centres between staggered levels.
    """
    pl = np.logspace(3.0, np.log10(ps), nlev)  # 1e3 Pa -> ps, ascending
    p = np.sqrt(pl[:-1] * pl[1:])  # geometric cell centres
    return SimpleNamespace(
        vol_list=dict(vol_list),
        ts=ts,
        ps=ps,
        pl=pl,
        tmpl=np.zeros(nlev),
        p=p,
    )


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_dry_adiabat_conserves_potential_temperature():
    """The Poisson adiabat holds potential temperature constant at ts.

    For a dry adiabat T = ts*(p/ps)**Rcp the potential temperature
    theta = T*(ps/p)**Rcp is invariant with height and equals the surface
    temperature ts (the analytical Poisson identity). The surface node p=ps
    must return exactly ts. An isothermal column (a "forgot the adiabat" bug)
    would instead give theta = ts*(ps/p)**Rcp, which departs from ts by a
    large factor at the top of the column, discriminating the pin.
    """
    ts = 1500.0
    atm = _make_atm({'N2': 1.0}, ts=ts)
    dry_adiabat_atm(atm)

    # Potential temperature reconstructed from the filled profile is constant.
    theta = atm.tmpl * (atm.ps / atm.pl) ** atm.Rcp
    np.testing.assert_allclose(theta, ts, rtol=1e-10)
    # Surface boundary condition: at p=ps the exponent base is 1, so T=ts.
    assert atm.tmpl[-1] == pytest.approx(ts, rel=1e-12)
    # Rcp for a diatomic-dominated column sits near R/cp ~ 2/7.
    assert 0.0 < atm.Rcp < 0.4

    # Discrimination: an isothermal profile (tmpl == ts) would yield a
    # potential temperature far above ts at the low-pressure top level.
    theta_if_isothermal = ts * (atm.ps / atm.pl) ** atm.Rcp
    assert theta_if_isothermal[0] > 1.5 * ts


@pytest.mark.physics_invariant
def test_dry_adiabat_temperature_rises_monotonically_toward_surface():
    """Dry-adiabat temperature increases with pressure and stays positive.

    On a Poisson adiabat with Rcp>0 the temperature rises monotonically from
    the cold, low-pressure top toward the hot surface, remains strictly
    positive, and the interpolated cell-centre profile stays bounded by the
    staggered end points. Rcp for a realistic volatile mixture must land in the
    physical R/cp band (0, 0.4), and the surface-to-top temperature ratio must
    exceed unity by a resolvable margin, ruling out an isothermal profile.
    """
    ts = 1200.0
    atm = _make_atm({'H2O': 0.5, 'CO2': 0.3, 'N2': 0.2}, ts=ts)
    dry_adiabat_atm(atm)

    # Monotone increase with pressure (pl ascends toward the surface).
    assert np.all(np.diff(atm.tmpl) > 0.0)
    assert np.all(atm.tmp > 0.0)
    # Surface boundary condition holds for the mixture too.
    assert atm.tmpl[-1] == pytest.approx(ts, rel=1e-12)
    # Interpolated cell-centre temperatures stay within the staggered bounds.
    assert atm.tmp.min() >= atm.tmpl.min()
    assert atm.tmp.max() <= atm.tmpl.max()
    # Physical R/cp band and a resolvable surface-to-top contrast.
    assert 0.0 < atm.Rcp < 0.4
    assert atm.tmpl[-1] / atm.tmpl[0] > 1.5
