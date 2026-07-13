"""Tests for src/janus/utils/water_tables.py.

Exercises the tabulated water thermodynamics lookup: pinned triple-point and
critical-point values, Clausius-Clapeyron consistency between the tabulated
latent heat and the saturation-curve slope, monotonicity, and the
beyond-range clamp. See docs/How-to/test.md.
"""

import numpy as np
import pytest

from janus.utils import water_tables as wt

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_lookup_pins_triple_point_and_critical_point():
    """The saturation curve pins the water triple point and critical point.

    The IAPWS reference values are p_tp = 611.657 Pa at 273.16 K and
    p_c = 22.064 MPa at 647.096 K (Wagner and Pruss 2002, J. Phys. Chem.
    Ref. Data 31, 387). Both are table endpoints, so an index-arithmetic
    slip in the linear-spacing lookup shifts them far outside tolerance.
    """
    p_triple = wt.lookup('psat', wt.T_tp)
    assert p_triple == pytest.approx(611.657, rel=1e-4)
    p_crit = wt.lookup('psat', wt.T_c)
    assert p_crit == pytest.approx(2.2064e7, rel=1e-3)
    # Scale guard: the two pins differ by more than four orders of
    # magnitude, so a table/units mix-up cannot satisfy both.
    assert p_crit / p_triple > 1e4
    # Sign guard: saturation pressure is strictly positive.
    assert p_triple > 0.0
    # Latent heat vanishes at the critical point and is ~2.5 MJ/kg at the
    # triple point; a column swap would exchange these magnitudes.
    assert wt.lookup('L_vap', wt.T_c) == pytest.approx(0.0, abs=1e3)
    assert wt.lookup('L_vap', wt.T_tp) == pytest.approx(2.5005e6, rel=1e-3)


@pytest.mark.physics_invariant
def test_lookup_monotonicity_and_clausius_clapeyron_consistency():
    """psat rises and L_vap falls with temperature; the tabulated slope
    matches L / (R_v T).

    The phase_grad column is documented as dln(psat)/dln(T), which for an
    ideal condensable equals L_vap / (R_v T); checking at 350 K ties the
    three tables together, so a swapped or rescaled column fails.
    """
    temps = np.linspace(wt.T_tp, wt.T_c, 25)
    psat_vals = wt.lookup('psat', temps)
    lvap_vals = wt.lookup('L_vap', temps)
    assert np.all(np.diff(psat_vals) > 0.0)
    assert np.all(np.diff(lvap_vals) < 0.0)
    assert np.all(np.isfinite(psat_vals))

    t_mid = 350.0
    rv = 461.52  # J kg-1 K-1, specific gas constant of water vapour
    expected_grad = wt.lookup('L_vap', t_mid) / (rv * t_mid)
    assert wt.lookup('phase_grad', t_mid) == pytest.approx(expected_grad, rel=2e-2)


def test_lookup_array_input_and_beyond_range_clamp():
    """Array temperatures vectorise and post-critical input extrapolates
    from the last interval without failing.

    Beyond-range behaviour is the documented no-fail contract of the
    lookup; the result must stay finite and continue the final segment's
    linear trend.
    """
    arr = np.array([300.0, 400.0, 500.0])
    vals = wt.lookup('psat', arr)
    assert vals.shape == (3,)
    assert np.all(np.diff(vals) > 0.0)

    beyond = wt.lookup('psat', 700.0)
    assert np.isfinite(beyond)
    # Linear extrapolation of the last interval must exceed the critical
    # pressure since the final segment slope is positive.
    assert beyond > wt.lookup('psat', wt.T_c)

    with pytest.raises(KeyError):
        wt.lookup('not_a_table', 300.0)
