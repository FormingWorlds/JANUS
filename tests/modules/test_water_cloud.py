"""Tests for src/janus/modules/water_cloud.py.

water_cloud.simple_cloud fills the three SOCRATES cloud input arrays (effective
radius, liquid water mass fraction, cloud fraction) wherever water condenses in
the column. It has two regimes selected by whether a condensation scheme has
already supplied a condensate field:

* No scheme (x_cond all zero): a level clouds when its temperature is at or below
  the local H2O dew point, taking the user cloud settings; otherwise the cloud
  arrays are zeroed.
* Scheme active (x_cond present): a level clouds where its condensate fraction is
  positive, and the liquid-water mass fraction is set to that condensate fraction.

Both regimes skip levels whose H2O partial pressure is negligible and never touch
the deepest level. This file exercises the saturated (clouding) case, the dry
(no-op) case, and the scheme-active case.

Invariant families exercised: boundedness of the cloud fraction and liquid-water
mass fraction in [0, 1], the dry-column no-op limit, and preservation of levels
the routine is contracted to skip. See docs/How-to/test.md.
"""

from types import SimpleNamespace

import numpy as np
import pytest

import janus.utils.GeneralAdiabat as ga
from janus.modules.water_cloud import simple_cloud

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Sentinel value on the cloud arrays before the call, so a level the routine is
# contracted to skip (negligible H2O, or the deepest level) is distinguishable
# from a level it processes and zeroes.
_SENTINEL = 99.0


def _cloud_atm(tmp, pp_h2o, x_cond, arr_init):
    """Build a minimal atmosphere stand-in carrying only what simple_cloud reads.

    arr_init seeds the re / lwm / clfr arrays so a skipped level keeps the seed
    while a processed level is overwritten with either the cloud settings or zero.
    """
    n = len(tmp)
    return SimpleNamespace(
        p=np.zeros(n),
        tmp=np.array(tmp, dtype=float),
        p_vol={'H2O': np.array(pp_h2o, dtype=float)},
        x_cond={'H2O': np.array(x_cond, dtype=float)},
        re=np.full(n, arr_init),
        lwm=np.full(n, arr_init),
        clfr=np.full(n, arr_init),
        effective_radius=1.0e-5,
        liquid_water_fraction=0.8,
        cloud_fraction=0.5,
    )


@pytest.mark.physics_invariant
def test_simple_cloud_saturated_column_activates_bounded_cloud():
    """A cold saturated column clouds every resolved level with valid settings.

    At 200 K every level with non-negligible H2O partial pressure sits below the
    local dew point (verified explicitly against ga.Tdew), so the routine writes
    the user cloud settings. The activated cloud fraction and liquid-water mass
    fraction must lie in [0, 1], and the effective radius must equal the supplied
    value. The second level carries a negligible H2O partial pressure, so the
    routine must skip it and leave its seeded sentinel intact, discriminating a
    broken guard that would cloud a water-free level.
    """
    n = 5
    pp = [1.0e4, 0.0, 1.0e4, 1.0e4, 1.0e4]  # index 1 is negligible -> skipped
    atm = _cloud_atm(np.full(n, 200.0), pp, np.zeros(n), _SENTINEL)

    simple_cloud(atm)

    # Precondition: the whole clouding column really is below the dew point.
    assert 200.0 < ga.Tdew('H2O', 1.0e4)
    # Level 0 is resolved and saturated -> user settings, bounded fractions.
    assert atm.clfr[0] == pytest.approx(0.5, rel=1e-12)
    assert 0.0 <= atm.clfr[0] <= 1.0
    assert 0.0 <= atm.lwm[0] <= 1.0
    assert atm.lwm[0] == pytest.approx(0.8, rel=1e-12)
    assert atm.re[0] == pytest.approx(1.0e-5, rel=1e-12)
    # Negligible-H2O level is skipped: the seeded sentinel survives untouched.
    assert atm.clfr[1] == pytest.approx(_SENTINEL, rel=1e-12)
    assert atm.re[1] == pytest.approx(_SENTINEL, rel=1e-12)


@pytest.mark.physics_invariant
def test_simple_cloud_dry_column_leaves_no_condensate():
    """A warm subsaturated column produces no cloud anywhere it is resolved.

    At 600 K every level is far above the local dew point, so the no-scheme
    branch zeroes the cloud arrays on all resolved levels. This is the dry-column
    no-op limit: the routine must not activate a cloud when nothing condenses. The
    deepest level is never visited by the loop, so its seeded value is left
    unchanged, an edge the pinned-zero levels are checked against.
    """
    n = 5
    atm = _cloud_atm(np.full(n, 600.0), np.full(n, 1.0e4), np.zeros(n), 0.0)
    atm.re[-1] = _SENTINEL  # deepest level is outside the loop range

    simple_cloud(atm)

    # Precondition: the column is genuinely subsaturated everywhere.
    assert 600.0 > ga.Tdew('H2O', 1.0e4)
    # Every resolved level is zeroed: no spurious condensate.
    assert np.all(atm.clfr[:-1] == pytest.approx(0.0, abs=1e-30))
    assert np.all(atm.lwm[:-1] == pytest.approx(0.0, abs=1e-30))
    # The unvisited deepest level keeps its seed, proving the loop bound.
    assert atm.re[-1] == pytest.approx(_SENTINEL, rel=1e-12)


@pytest.mark.physics_invariant
def test_simple_cloud_uses_condensate_fraction_when_scheme_active():
    """With a condensation scheme, cloud water follows the condensate fraction.

    When x_cond carries a non-trivial condensate field the routine takes its
    scheme-active branch: each resolved level with a positive condensate fraction
    is clouded with its liquid-water mass fraction set to that fraction, and each
    level with zero condensate is zeroed. The liquid-water mass fractions must
    equal the input condensate fractions and stay bounded in [0, 1], while the
    cloud fraction takes the fixed user setting. A level with negligible H2O
    partial pressure is skipped even when its condensate fraction is positive (the
    guard limit), and a level with zero condensate is zeroed rather than clouded.
    """
    n = 5
    # Index 1 has a positive condensate fraction but negligible H2O partial
    # pressure, so the routine must skip it; index 3 is a zeroed resolved level.
    x_cond = [0.3, 0.5, 0.5, 0.0, 0.2]
    pp = [1.0e4, 0.0, 1.0e4, 1.0e4, 1.0e4]
    atm = _cloud_atm(np.full(n, 250.0), pp, x_cond, 0.0)

    simple_cloud(atm)

    # Condensate-bearing resolved levels adopt their condensate fraction as lwm.
    assert atm.lwm[0] == pytest.approx(0.3, rel=1e-12)
    assert atm.lwm[2] == pytest.approx(0.5, rel=1e-12)
    assert np.all((atm.lwm[:-1] >= 0.0) & (atm.lwm[:-1] <= 1.0))
    assert atm.clfr[0] == pytest.approx(0.5, rel=1e-12)
    # Negligible-H2O level is skipped despite a positive condensate fraction.
    assert atm.lwm[1] == pytest.approx(0.0, abs=1e-30)
    assert atm.clfr[1] == pytest.approx(0.0, abs=1e-30)
    # A zero-condensate resolved level is zeroed, not clouded.
    assert atm.lwm[3] == pytest.approx(0.0, abs=1e-30)
    assert atm.clfr[3] == pytest.approx(0.0, abs=1e-30)
