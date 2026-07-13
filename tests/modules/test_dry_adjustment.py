"""Tests for src/janus/modules/dry_adjustment.py.

dry_adjustment.py provides ``DryAdj``, a single-sweep dry convective adjustment
that relaxes super-adiabatic level pairs toward the dry adiabat with equal layer
masses. Each pairwise swap conserves the layer-summed temperature. The routine
runs a downward sweep followed by an upward sweep, because a downward fix can
raise a level and destabilise the already-processed pair above it. This file
exercises:

* Layer-sum conservation across the adjustment (the equal-mass swap invariant).
* The upward sweep: it restabilises the top pair that the downward sweep alone
  leaves super-adiabatic, so the two-sweep result differs from a downward-only
  result and lowers the residual super-adiabatic excess.
* The no-op contract on an already dry-stable and on a dry-neutral column.

Invariant families exercised: conservation (layer-summed temperature),
positivity, monotonicity of the stability excess, and
pinned-value-with-discrimination-guard. See docs/How-to/test.md.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from janus.modules.dry_adjustment import DryAdj

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

RCP = 2.0 / 7.0  # diatomic R/cp


def _downward_only(temp, press, rcp):
    """Replicate a single downward adjustment sweep (no upward sweep).

    This is the incomplete algorithm the upward sweep exists to complete; the
    two-sweep ``DryAdj`` output is compared against it as the discrimination
    guard.
    """
    temp = temp.copy()
    for i in range(len(temp) - 1):
        pfact = (press[i] / press[i + 1]) ** rcp
        if temp[i] < temp[i + 1] * pfact:
            tbar = 0.5 * (temp[i] + temp[i + 1])
            temp[i + 1] = 2.0 * tbar / (1.0 + pfact)
            temp[i] = temp[i + 1] * pfact
    return temp


def _superadiabatic_excess(temp, press, rcp, i):
    """Super-adiabatic excess of pair (i, i+1); positive means unstable."""
    pfact = (press[i] / press[i + 1]) ** rcp
    return temp[i + 1] * pfact - temp[i]


@pytest.mark.physics_invariant
def test_dry_adjustment_conserves_layer_sum_and_upward_pass_restabilises_top():
    """The upward sweep restabilises the top pair while conserving the layer sum.

    The initial column is dry-stable at the top pair but super-adiabatic at the
    bottom pair. The downward sweep fixes the bottom pair, which raises the
    shared middle level and drives the top pair super-adiabatic; only the upward
    sweep can then restabilise it. The two-sweep result therefore differs from a
    downward-only sweep, raises the top level the downward sweep leaves
    untouched, and lowers the top-pair super-adiabatic excess. Every pairwise
    swap uses equal layer masses, so the layer-summed temperature is conserved.
    """
    press = np.array([1.0e4, 3.16e4, 1.0e5])  # ascending toward the surface
    temp = np.array([150.0, 200.0, 300.0])  # stable (0,1); unstable (1,2)
    assert _superadiabatic_excess(temp, press, RCP, 0) < 0.0  # top pair stable
    assert _superadiabatic_excess(temp, press, RCP, 1) > 0.0  # bottom pair unstable

    downward = _downward_only(temp, press, RCP)
    atm = SimpleNamespace(tmp=temp.copy(), p=press.copy(), Rcp=RCP)
    DryAdj(atm)
    both = atm.tmp

    # Conservation: equal-mass swaps preserve the layer-summed temperature.
    assert both.sum() == pytest.approx(temp.sum(), rel=1e-12)
    assert np.all(both > 0.0)
    # The downward sweep leaves the top level untouched at its initial 150 K and
    # the top pair super-adiabatic; the upward sweep raises the top level and
    # neutralises that excess.
    assert downward[0] == pytest.approx(150.0, rel=1e-12)
    assert _superadiabatic_excess(downward, press, RCP, 0) > 0.1
    assert both[0] > downward[0] + 0.1
    assert _superadiabatic_excess(both, press, RCP, 0) < _superadiabatic_excess(
        downward, press, RCP, 0
    )
    # Discrimination: dropping the upward sweep changes the top-pair result, so
    # the two profiles are not equal.
    assert not np.allclose(both, downward)


@pytest.mark.physics_invariant
def test_dry_adjustment_leaves_stable_and_neutral_columns_unchanged():
    """An already dry-stable or dry-neutral column is returned unchanged.

    ``DryAdj`` only touches pairs that are super-adiabatic. A column whose
    temperature rises toward the surface more slowly than the dry adiabat (every
    pair sub-adiabatic) must pass through untouched, and a column exactly on the
    adiabat (zero excess everywhere, the neutral edge case) must stay on it. A
    spurious "always adjust" path would perturb both.
    """
    press = np.array([1.0e4, 3.16e4, 1.0e5])

    # Strongly stable column: every pair is sub-adiabatic (negative excess).
    stable = np.array([280.0, 285.0, 290.0])
    for i in range(len(stable) - 1):
        assert _superadiabatic_excess(stable, press, RCP, i) < 0.0
    atm_stable = SimpleNamespace(tmp=stable.copy(), p=press.copy(), Rcp=RCP)
    DryAdj(atm_stable)
    np.testing.assert_allclose(atm_stable.tmp, stable, rtol=0, atol=0)

    # Neutral column: built exactly on the adiabat, so the excess is zero and no
    # adjustment triggers.
    ts = 300.0
    pf12 = (press[1] / press[2]) ** RCP
    pf01 = (press[0] / press[1]) ** RCP
    neutral = np.array([ts * pf12 * pf01, ts * pf12, ts])
    atm_neutral = SimpleNamespace(tmp=neutral.copy(), p=press.copy(), Rcp=RCP)
    DryAdj(atm_neutral)
    np.testing.assert_allclose(atm_neutral.tmp, neutral, rtol=1e-12, atol=0)
