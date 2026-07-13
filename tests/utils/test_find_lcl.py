"""Tests for src/janus/utils/find_lcl.py.

find_lcl.py provides ``find_intersection``, which locates the lifting
condensation level (LCL) of a column as the level where the water partial
pressure meets the saturation curve. This file exercises:

* The saturated-surface limit, where the partial-pressure and saturation
  curves already coincide at the surface, so the LCL is the surface itself.
* The monotonic relationship between dryness and LCL altitude: a drier parcel
  reaches saturation only higher up, at a lower pressure.
* The multiple-crossing tie-break: with more than one crossing the routine returns
  the lowest-altitude (last-index) one, which the pinned test records despite the
  source comment naming it the "first occurrence".
* The no-intersection error contract, where no level is within tolerance and the
  routine reports the failure and returns None.

Invariant families exercised: monotonicity of the LCL with parcel dryness,
boundedness of the returned index, and pinned-value discrimination on the
returned saturation pressure. See docs/How-to/test.md.
"""

import numpy as np
import pytest

from janus.utils.find_lcl import find_intersection

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Column pressure grid, ascending from the top of atmosphere (index 0) to the
# surface (index -1), in Pa. The saturation curve rises toward the warm surface.
PRESSURE = np.array([1.0e3, 1.0e4, 5.0e4, 1.0e5])
P_SAT = np.array([1.0e2, 1.0e3, 2.0e4, 8.0e4])


def test_saturated_surface_parcel_condenses_at_the_surface():
    """A parcel already saturated at the surface has its LCL at the surface.

    When the water partial pressure equals the saturation pressure only at the
    surface, ``find_intersection`` returns the surface level and its saturation
    pressure. The returned index must be the last (surface) node and the returned
    value must equal the surface saturation pressure, discriminating a routine
    that mistakenly reports the top-of-atmosphere crossing.
    """
    # Partial pressure meets the saturation curve only at the surface node.
    partial = np.array([5.0e3, 5.0e3, 5.0e3, 8.0e4])
    idx, value = find_intersection(P_SAT, partial, tolerance=1.0)
    # The LCL sits at the surface (the highest-pressure, last node).
    assert idx == len(P_SAT) - 1
    # The returned value is the surface saturation pressure, pinned exactly.
    assert value == pytest.approx(P_SAT[-1], rel=1e-12)
    # Discrimination: the surface value is far from the top-of-atmosphere
    # saturation pressure, so an index slip would be caught.
    assert value > P_SAT[0] * 100.0


@pytest.mark.physics_invariant
def test_drier_parcel_reaches_its_lifting_condensation_level_higher_up():
    """A drier parcel condenses at a lower pressure (higher altitude) than a moist one.

    The LCL rises as a parcel dries out: with less water it must be lifted and
    cooled further before it saturates. Both parcels here have a single unambiguous
    crossing of the saturation curve (the moist one at the surface, index 3; the
    drier one aloft, index 1), so the comparison is not confounded by the tie-break
    rule. The drier parcel's LCL must sit at a strictly lower index and lower
    pressure.
    """
    # Each parcel touches the saturation curve at exactly one level.
    partial_moist = np.array([5.0e3, 5.0e3, 5.0e3, 8.0e4])  # single crossing at surface idx 3
    partial_dry = np.array([5.0e3, 1.0e3, 5.0e3, 5.0e3])  # single crossing aloft at idx 1
    idx_moist, _ = find_intersection(P_SAT, partial_moist, tolerance=1.0)
    idx_dry, val_dry = find_intersection(P_SAT, partial_dry, tolerance=1.0)
    # Drier air saturates higher up: lower index, lower pressure.
    assert idx_dry == 1 and idx_moist == 3  # single, unambiguous crossings
    assert PRESSURE[idx_dry] < PRESSURE[idx_moist]
    # The drier LCL pressure is well below the surface value, not a rounding-scale
    # difference; pin the recovered saturation pressure there.
    assert val_dry == pytest.approx(P_SAT[idx_dry], rel=1e-12)
    assert val_dry > 0.0


def test_multiple_crossings_return_the_lowest_level_crossing():
    """With two crossings the routine returns the lowest-altitude (last-index) one.

    A parcel can meet the saturation curve at more than one level. The source scans
    the whole column and returns the highest-index (lowest-altitude, closest to the
    surface) crossing: it takes ``intersection_indices[-1]``. Note that the inline
    source comment and the ``first_intersection_index`` variable name say "first
    occurrence", which contradicts the actual last-index behavior pinned here. This
    test fixes the observed behavior so a later change to the tie-break is caught.
    """
    # Partial pressure meets the saturation curve at index 0 and index 2.
    partial_two = np.array([1.0e2, 5.0e3, 2.0e4, 5.0e3])
    idx, value = find_intersection(P_SAT, partial_two, tolerance=1.0)
    # The returned crossing is the last index (2), not the first (0).
    assert idx == 2
    assert idx != 0  # discriminates the actual behavior from the "first occurrence" comment
    # The value is the saturation pressure at that lower level, pinned exactly.
    assert value == pytest.approx(P_SAT[2], rel=1e-12)


def test_no_condensation_level_reports_failure_and_returns_none(caplog):
    """With no level within tolerance the routine logs an error and returns None.

    If the partial-pressure curve never approaches the saturation curve, there is
    no lifting condensation level. The routine must report this rather than
    silently returning a spurious index: it logs an error and returns None. This
    is the error contract, exercised by an everywhere-subsaturated parcel.
    """
    # Partial pressure is one megapascal below saturation at every level.
    partial_none = P_SAT - 1.0e6
    with caplog.at_level('ERROR'):
        result = find_intersection(P_SAT, partial_none, tolerance=1.0)
    assert result is None
    # The failure is surfaced in the log, not swallowed.
    assert any('No LCL found' in record.message for record in caplog.records)
