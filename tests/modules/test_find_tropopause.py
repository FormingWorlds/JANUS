"""Tests for src/janus/modules/find_tropopause.py.

find_tropopause locates the tropopause of a moist column by one of two criteria:
a dynamic search driven by the net radiative heating profile, or a fixed
temperature threshold. This file exercises:

* The dynamic branch when a clear heating maximum sits above a stratospheric
  sign change: the tropopause index is placed where the heating first drops
  below the significance threshold.
* The dynamic guard paths that reject a candidate that is too shallow (closer to
  the model top than the stratosphere depth) or too weak (mean heating above it
  below the significance floor); both leave the column untouched.
* The fixed-temperature branch: the cold-point of a stratospheric inversion, the
  make-safe clamp when the threshold exceeds the surface temperature, and the
  not-triggered limit that returns the sentinel (idx -1, P 1e-10).

Invariant families exercised: positivity/boundedness of the returned tropopause
pressure and temperature, index-within-bounds after the clamp, and the degenerate
limits that leave the column intact. See docs/How-to/test.md.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from janus.modules.find_tropopause import find_tropopause

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Sentinel values placed on the tropopause fields before each call so an
# untouched column is distinguishable from one where the routine wrote a result.
_IDX_SENTINEL = -999
_P_SENTINEL = -999.0


def _dynamic_atm(net_heating, trppT=180.0):
    """Build a minimal moist-atmosphere stand-in for the dynamic branch.

    Only the attributes read by the dynamic path are populated. The pressure and
    temperature edge grids descend with index so a written trppP / trppT are
    strictly positive and easy to pin.
    """
    n = len(net_heating)
    pl = np.linspace(1.0e5, 1.0, n)
    tmpl = np.linspace(300.0, 150.0, n)
    return SimpleNamespace(
        net_heating=np.array(net_heating, dtype=float),
        p=np.zeros(n),
        pl=pl,
        tmpl=tmpl,
        tmp=tmpl.copy(),
        trppT=trppT,
        trppidx=_IDX_SENTINEL,
        trppP=_P_SENTINEL,
    )


def _fixed_atm(tmpl, trppT):
    """Build a minimal moist-atmosphere stand-in for the fixed-threshold branch.

    tmp is the cell-centre array (one shorter than the edge array tmpl), so the
    make-safe clamp against ``len(tmp) - 1`` has a shorter target than the edge
    grid and can actually reduce the found index.
    """
    tmpl = np.array(tmpl, dtype=float)
    pl = np.logspace(5.0, 2.0, len(tmpl))
    tmp = (tmpl[:-1] + tmpl[1:]) / 2.0
    return SimpleNamespace(
        pl=pl,
        tmpl=tmpl,
        tmp=tmp,
        trppT=trppT,
        trppidx=_IDX_SENTINEL,
        trppP=_P_SENTINEL,
        net_heating=np.zeros(len(tmp)),
        p=np.zeros(len(tmp)),
    )


@pytest.mark.physics_invariant
def test_dynamic_tropopause_located_at_heating_dropoff():
    """A strong heating maximum above a sign change fixes the tropopause index.

    The heating profile holds 60 K/day over the three deepest levels, drops to
    5 K/day, then turns negative near the top (the sign change that qualifies the
    search). The routine walks down from the heating maximum while heating stays
    above the 30 K/day significance threshold, so the tropopause lands at the
    first level where heating falls to or below 30. The exact index is pinned and
    cross-checked against the physical criterion (heating just above it exceeds
    30, heating at it does not), and the written trppP / trppT are strictly
    positive, ruling out an off-by-one that would grab a supra-threshold level.
    """
    net_heating = [60.0, 60.0, 60.0] + [5.0] * 90 + [-10.0] * 7
    atm = _dynamic_atm(net_heating)

    find_tropopause(atm, dynamic=True)

    assert atm.trppidx == 3
    # The index is exactly where heating first drops out of the significant band.
    assert net_heating[atm.trppidx - 1] > 30.0
    assert net_heating[atm.trppidx] <= 30.0
    # Positivity: the tropopause pressure and temperature come from the descending
    # edge grids, so both must be strictly positive at the located index.
    assert atm.trppP == pytest.approx(atm.pl[3], rel=1e-12)
    assert atm.trppT == pytest.approx(atm.tmpl[3], rel=1e-12)
    assert atm.trppP > 0.0
    assert atm.trppT > 0.0


def test_dynamic_no_significant_heating_leaves_column_intact():
    """Sub-threshold heating everywhere writes no tropopause result.

    With the whole profile at a uniform 5 K/day the peak heating never exceeds
    the 50 K/day significance ceiling, so the search sets no index and the
    tropopause fields keep their pre-call sentinel values. This is the degenerate
    limit where no tropopause exists; a regression that unconditionally wrote a
    result would overwrite the sentinels.
    """
    atm = _dynamic_atm([5.0] * 100)

    find_tropopause(atm, dynamic=True)

    assert atm.trppidx == _IDX_SENTINEL
    assert atm.trppP == pytest.approx(_P_SENTINEL, rel=1e-12)


@pytest.mark.parametrize(
    'net_heating',
    [
        # Candidate one level below the top: shallower than the stratosphere
        # depth (2 levels for a 100-level column), so it is rejected.
        [60.0] + [5.0] * 92 + [-10.0] * 7,
        # Candidate deep enough, but the mean heating above it is strongly
        # negative and falls below the 10 K/day floor, so it is rejected.
        [-100.0, -100.0, 60.0, 60.0] + [5.0] * 89 + [-10.0] * 7,
    ],
    ids=['too-shallow', 'weak-mean-heating'],
)
def test_dynamic_tropopause_rejected_when_shallow_or_weak(net_heating):
    """Shallow or weak candidates are rejected and leave the column intact.

    Both profiles clear the significance ceiling and contain a sign change, so
    the inner search runs, but each candidate is then discarded: one is too close
    to the model top, the other has insufficient mean heating above it. In both
    cases the routine must fall back to writing nothing, so the sentinels survive.
    """
    atm = _dynamic_atm(net_heating)

    find_tropopause(atm, dynamic=True)

    assert atm.trppidx == _IDX_SENTINEL
    assert atm.trppP == pytest.approx(_P_SENTINEL, rel=1e-12)


@pytest.mark.physics_invariant
def test_fixed_temperature_tropopause_at_inversion_coldpoint():
    """The fixed threshold marks the cold-point of a stratospheric inversion.

    The edge temperature profile falls to a 150 K minimum then warms again above
    it (a stratospheric inversion). With a 180 K threshold only the cold-point
    level lies below it, so the upward scan stops there. The located level sits
    below the threshold while the level just above it is warmer than the
    threshold, confirming the routine found the inversion minimum rather than an
    arbitrary sub-threshold level; the returned pressure stays strictly positive.
    """
    atm = _fixed_atm([300.0, 250.0, 200.0, 150.0, 200.0, 250.0, 300.0], trppT=180.0)

    find_tropopause(atm, dynamic=False)

    assert atm.trppidx == 3
    assert atm.tmpl[atm.trppidx] < atm.trppT
    # The level above the cold-point is warmer than the threshold: the inversion.
    assert atm.tmpl[atm.trppidx + 1] >= atm.trppT
    assert atm.trppP == pytest.approx(atm.pl[3], rel=1e-12)
    assert atm.trppP > 0.0


@pytest.mark.physics_invariant
def test_fixed_temperature_threshold_above_surface_clamps_index():
    """A threshold above the surface temperature clamps the index in-bounds.

    When the fixed threshold exceeds every column temperature the upward scan
    matches at the model-top edge (index ``len(tmpl) - 1``). The make-safe clamp
    limits the index to ``len(tmp) - 1`` (the shorter cell-centre grid), so the
    returned index stays a valid centre index and the returned pressure remains
    strictly positive. Without the clamp the index would over-run the centre
    arrays.
    """
    atm = _fixed_atm([300.0, 250.0, 200.0, 150.0, 120.0, 100.0, 90.0], trppT=500.0)

    find_tropopause(atm, dynamic=False)

    assert atm.trppidx == len(atm.tmp) - 1
    assert atm.trppidx <= len(atm.tmpl) - 1
    assert atm.trppP == pytest.approx(atm.pl[atm.trppidx], rel=1e-12)
    assert atm.trppP > 0.0


def test_fixed_temperature_threshold_touched_but_not_crossed():
    """A threshold met only by equality triggers the search but sets no new level.

    When an edge temperature equals the threshold exactly, the trigger condition
    (temperature at or below the threshold) is satisfied, so the fixed branch runs
    its upward scan, but no level lies strictly below the threshold, so the scan
    completes without selecting a new level. The make-safe clamp then leaves the
    pre-existing index in place and reads back a consistent pressure. This is the
    boundary case that separates "reached the threshold" from "dropped below it".
    """
    # Index 0 is the default index carried in from construction.
    atm = _fixed_atm([300.0, 250.0, 300.0], trppT=250.0)
    atm.trppidx = 0

    find_tropopause(atm, dynamic=False)

    assert atm.trppidx == 0
    assert atm.trppP == pytest.approx(atm.pl[0], rel=1e-12)


def test_fixed_temperature_no_crossing_sets_untriggered_sentinel():
    """No level below the threshold returns the untriggered sentinel result.

    With a 50 K threshold that no level reaches, the fixed branch takes its
    not-triggered path and writes the documented sentinel: index -1 and a tiny
    positive pressure floor of 1e-10 Pa. This is the error-contract limit that
    signals "no tropopause" to the caller.
    """
    atm = _fixed_atm([300.0, 250.0, 200.0, 150.0, 120.0, 100.0, 90.0], trppT=50.0)

    find_tropopause(atm, dynamic=False)

    assert atm.trppidx == -1
    assert atm.trppP == pytest.approx(1.0e-10, rel=1e-6)
