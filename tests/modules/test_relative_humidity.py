"""Tests for src/janus/modules/relative_humidity.py.

relative_humidity.py provides ``compute_Rh``, which returns the relative humidity
profile (water partial pressure over saturation pressure), and
``update_for_constant_RH``, which rewrites the water partial-pressure profile so
the column holds a prescribed relative humidity, then rebuilds the total and
surface pressures by Dalton summation. This file exercises:

* ``compute_Rh``: the partial-over-saturation ratio, its unity value at
  saturation, and the current unclamped divergence when saturation vanishes.
* The constant-relative-humidity rewrite: the new water partial pressure equals
  RH times the saturation pressure at every level.
* Dalton's law closure: the total pressure equals the sum of the component
  partial pressures, and the surface pressure equals that sum at the base.
* The dry and saturated relative-humidity limits (RH = 0 and RH = 1).
* The current unguarded crash when no lifting condensation level exists (a known
  defect, recorded rather than designed).

Invariant families exercised: Dalton-summation conservation,
positivity/boundedness of the rebuilt pressures, and pinned-value discrimination
against the saturation curve. See docs/How-to/test.md.
"""

from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from janus.modules.relative_humidity import compute_Rh, update_for_constant_RH

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# A saturation pressure that scales linearly with temperature keeps the pinned
# values easy to derive by hand while still discriminating a p vs p_sat mix-up.
_PSAT_SLOPE = 50.0  # Pa per K


def _fake_p_sat(_molecule, temperature):
    """Deterministic stand-in for ga.p_sat, linear in temperature."""
    return _PSAT_SLOPE * float(temperature)


def _make_atm(rh_partial):
    """Build a synthetic three-level column whose water curve meets saturation at the base."""
    tmp = np.array([200.0, 250.0, 300.0])
    p_sat = _PSAT_SLOPE * tmp  # [10000, 12500, 15000] Pa
    # Water partial pressure meets the saturation curve at the surface so the
    # lifting-condensation-level search finds an intersection there.
    h2o = np.array([5000.0, 5000.0, p_sat[-1]]) if rh_partial else np.array([1.0, 1.0, 1.0])
    co2 = np.array([100.0, 200.0, 300.0])
    return SimpleNamespace(
        tmp=tmp,
        p=h2o + co2,
        p_vol={'H2O': h2o.copy(), 'CO2': co2.copy()},
        vol_list={'H2O': 0.5, 'CO2': 0.5},
        ps=float((h2o + co2)[-1]),
    )


@pytest.mark.physics_invariant
def test_compute_rh_is_partial_pressure_over_saturation():
    """compute_Rh returns the water partial pressure divided by its saturation pressure.

    Relative humidity is the ratio of the water partial pressure (x_gas['H2O'] * p)
    to the local saturation pressure. The values are chosen so the wrong-order
    division (p_sat over partial pressure) lands far outside tolerance, and a
    subsaturated column returns a humidity strictly below one.
    """
    tmp = np.array([200.0, 250.0, 300.0])
    p = np.array([1.0e4, 2.0e4, 4.0e4])
    x_h2o = np.array([0.1, 0.2, 0.25])  # partial_p = x_gas * p = [1000, 4000, 10000] Pa
    p_sat_map = {200.0: 4000.0, 250.0: 8000.0, 300.0: 40000.0}
    atm = SimpleNamespace(tmp=tmp, p=p, x_gas={'H2O': x_h2o})
    with patch(
        'janus.modules.relative_humidity.ga.p_sat',
        side_effect=lambda _mol, t: p_sat_map[float(t)],
    ):
        r_h = compute_Rh(atm)
    # R_h = partial_p / p_sat = [0.25, 0.5, 0.25].
    expected = np.array([1000.0, 4000.0, 10000.0]) / np.array([4000.0, 8000.0, 40000.0])
    np.testing.assert_allclose(r_h, expected, rtol=1e-12)
    # Wrong-order guard: p_sat / partial_p would give 4.0 at the surface level,
    # sixteen times the correct 0.25 and far outside any rounding tolerance.
    assert r_h[0] == pytest.approx(0.25, rel=1e-12)
    assert abs(r_h[0] - (p_sat_map[200.0] / 1000.0)) > 1.0
    # Subsaturated air everywhere: relative humidity strictly below one.
    assert np.all(r_h < 1.0)


@pytest.mark.physics_invariant
def test_compute_rh_reaches_unity_at_saturation():
    """compute_Rh returns exactly one when the water partial pressure equals saturation.

    At saturation the water partial pressure equals the saturation pressure, so the
    relative humidity is one at every level. This is the physical boundary between
    subsaturated and supersaturated air; a wrong normalisation would miss unity.
    """
    tmp = np.array([220.0, 260.0, 300.0])
    p = np.array([1.0e4, 2.0e4, 5.0e4])
    x_h2o = np.array([0.3, 0.3, 0.3])
    partial = x_h2o * p  # [3000, 6000, 15000] Pa, all non-trivial
    p_sat_map = {220.0: 3000.0, 260.0: 6000.0, 300.0: 15000.0}
    atm = SimpleNamespace(tmp=tmp, p=p, x_gas={'H2O': x_h2o})
    with patch(
        'janus.modules.relative_humidity.ga.p_sat',
        side_effect=lambda _mol, t: p_sat_map[float(t)],
    ):
        r_h = compute_Rh(atm)
    np.testing.assert_allclose(r_h, np.ones(3), rtol=1e-12)
    # The partial pressures are genuinely positive, so unity is a real saturation
    # result and not a degenerate zero-over-zero.
    assert np.all(partial > 0.0)


def test_compute_rh_diverges_as_saturation_vanishes():
    """compute_Rh diverges to infinity as the saturation pressure goes to zero.

    This records the current behavior, not a designed contract: with no clamp, a
    vanishing saturation pressure at a very cold level makes the ratio partial over
    saturation divide by zero, so compute_Rh returns positive infinity there. The
    warmer level with a finite saturation pressure stays finite, isolating the
    divergence to the cold node.
    """
    tmp = np.array([100.0, 300.0])
    p = np.array([1.0e4, 4.0e4])
    x_h2o = np.array([0.2, 0.2])  # partial_p = [2000, 8000] Pa
    p_sat_map = {100.0: 0.0, 300.0: 40000.0}  # cold-level saturation collapses to zero
    atm = SimpleNamespace(tmp=tmp, p=p, x_gas={'H2O': x_h2o})
    with patch(
        'janus.modules.relative_humidity.ga.p_sat',
        side_effect=lambda _mol, t: p_sat_map[float(t)],
    ):
        with np.errstate(divide='ignore'):  # the unclamped divide-by-zero is the point
            r_h = compute_Rh(atm)
    # Cold level: partial_p / 0 -> +inf (unclamped, current behavior).
    assert np.isinf(r_h[0])
    assert r_h[0] > 0.0  # positive infinity, from a positive partial pressure
    # The warmer level with finite saturation stays finite and bounded.
    assert np.isfinite(r_h[1])
    assert r_h[1] == pytest.approx(8000.0 / 40000.0, rel=1e-12)


@pytest.mark.physics_invariant
def test_constant_relative_humidity_rewrites_water_and_closes_dalton_sum():
    """A constant-RH update sets water to RH*p_sat and rebuilds pressures by Dalton's law.

    ``update_for_constant_RH`` overwrites the water partial pressure with
    RH*p_sat at every level, then rebuilds the total pressure as the sum of all
    component partial pressures and the surface pressure as that sum at the base.
    The rewrite is pinned against the saturation curve (discriminating a routine
    that scales the total pressure instead), and the total pressure must equal the
    component sum exactly (Dalton conservation).
    """
    atm = _make_atm(rh_partial=True)
    p_sat = _PSAT_SLOPE * atm.tmp
    co2 = atm.p_vol['CO2'].copy()
    r_h = np.full(3, 0.5)  # uniform 50 percent relative humidity
    with patch('janus.modules.relative_humidity.ga.p_sat', side_effect=_fake_p_sat):
        out = update_for_constant_RH(atm, r_h)
    # The water profile is RH times the saturation curve, not RH times total p.
    np.testing.assert_allclose(out.p_vol['H2O'], 0.5 * p_sat, rtol=1e-12)
    # Surface pin: 0.5 * 15000 Pa, well away from a total-pressure scaling.
    assert out.p_vol['H2O'][-1] == pytest.approx(7500.0, rel=1e-12)
    # Dalton's law: total pressure is the sum of the component partial pressures.
    np.testing.assert_allclose(out.p, out.p_vol['H2O'] + co2, rtol=1e-12)
    # Surface pressure is the base-level Dalton sum and is a positive scalar.
    assert out.ps == pytest.approx(float((0.5 * p_sat + co2)[-1]), rel=1e-12)
    assert isinstance(out.ps, float)
    assert np.all(out.p > 0.0)


@pytest.mark.physics_invariant
def test_relative_humidity_dry_and_saturated_limits():
    """The RH = 0 and RH = 1 limits give a dry and a fully saturated water profile.

    At zero relative humidity the water partial pressure collapses to zero and the
    total pressure reduces to the non-condensable background; at unit relative
    humidity the water partial pressure equals the saturation curve exactly. These
    two boundary inputs frame the physical range and confirm the rebuilt pressures
    stay positive.
    """
    p_sat_expected = _PSAT_SLOPE * np.array([200.0, 250.0, 300.0])

    atm_dry = _make_atm(rh_partial=True)
    co2 = atm_dry.p_vol['CO2'].copy()
    with patch('janus.modules.relative_humidity.ga.p_sat', side_effect=_fake_p_sat):
        dry = update_for_constant_RH(atm_dry, np.zeros(3))
    # No water anywhere: total pressure is the background alone.
    np.testing.assert_allclose(dry.p_vol['H2O'], np.zeros(3), atol=1e-30)
    np.testing.assert_allclose(dry.p, co2, rtol=1e-12)
    assert dry.ps == pytest.approx(float(co2[-1]), rel=1e-12)

    atm_sat = _make_atm(rh_partial=True)
    with patch('janus.modules.relative_humidity.ga.p_sat', side_effect=_fake_p_sat):
        sat = update_for_constant_RH(atm_sat, np.ones(3))
    # Fully saturated: water sits exactly on the saturation curve.
    np.testing.assert_allclose(sat.p_vol['H2O'], p_sat_expected, rtol=1e-12)
    assert np.all(sat.p_vol['H2O'] > dry.p_vol['H2O'])
    assert np.all(sat.p > 0.0)


def test_update_crashes_when_no_lifting_condensation_level_exists():
    """Document the current unguarded crash when no lifting condensation level exists.

    This is not a designed contract: it records a known defect. The rewrite locates
    the lifting condensation level by unpacking ``find_intersection``, which returns
    None (a single value, not a pair) when no level is within tolerance. There is no
    guard for the missing-LCL case, so on an everywhere-subsaturated column the
    unpack raises TypeError. The crash happens before the water profile is
    overwritten, so the input column is left untouched.
    """
    atm = _make_atm(rh_partial=False)  # water pressures far below saturation
    original_h2o = atm.p_vol['H2O'].copy()
    with patch('janus.modules.relative_humidity.ga.p_sat', side_effect=_fake_p_sat):
        with pytest.raises(TypeError):
            update_for_constant_RH(atm, np.full(3, 0.5))
    # The failure occurs before the water profile is overwritten, so the input
    # column is left untouched.
    np.testing.assert_array_equal(atm.p_vol['H2O'], original_h2o)
