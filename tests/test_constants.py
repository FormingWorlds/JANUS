"""Cross-cutting sanity checks on JANUS physical constants and simple laws.

This cross-cutting file (a documented exception to the mirror-the-source rule)
spot-checks three quantities read directly from the constants tables: the
inverse-square surface-gravity law, the tabulated Earth reference constants in
``janus.utils.planets``, and the constant-mode heat capacity ``cpv``. The
per-source deep tests live in tests/utils/test_height.py and
tests/utils/test_cp_funcs.py; this file guards the values a caller reads
straight from the tables. See docs/How-to/test.md.
"""

import numpy as np
import pytest

from janus.utils.cp_funcs import cpv
from janus.utils.height import gravity
from janus.utils.planets import Earth

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.physics_invariant
def test_gravity_decreases_with_radius():
    """Surface gravity falls as 1/r^2 with planetary radius at fixed mass.

    ``gravity`` is the Newtonian inverse-square law, so at fixed mass a larger
    radius gives weaker surface gravity, and doubling the radius must quarter g
    exactly, not merely lower it. The exact-quarter check discriminates an
    inverse-linear 1/r implementation, which would only halve g.
    """
    m = 5.972e24
    r_earth = 6.371e6
    g_surface = gravity(m, r_earth)
    g_double = gravity(m, 2.0 * r_earth)
    # Monotone decrease with radius.
    assert g_surface > g_double
    # Inverse-square, not inverse-linear: 2x radius -> exactly g/4.
    assert g_double == pytest.approx(g_surface / 4.0, rel=1e-12)
    assert g_double < g_surface / 3.0  # rules out the inverse-linear g/2 result
    # Positivity.
    assert g_surface > 0.0


def test_earth_constants_match_database():
    """The tabulated Earth constants hold their documented reference values.

    ``janus.utils.planets.Earth`` pins surface gravity, Bond albedo, and mean
    surface temperature. Each is pinned with a tolerance rather than exact
    float equality and bracketed so a units slip (fraction vs percent albedo,
    degrees Celsius vs kelvin, m/s^2 vs a 10x error) would fall outside the
    guard band.
    """
    # Surface gravity ~9.8 m/s^2, not 0.98 or 98.
    assert Earth.g == pytest.approx(9.798, rel=1e-6)
    assert 9.0 < Earth.g < 10.0
    # Bond albedo is a dimensionless fraction in (0, 1), ~0.306, not 30.6.
    assert Earth.albedo == pytest.approx(0.306, rel=1e-6)
    assert 0.0 < Earth.albedo < 1.0
    # Mean surface temperature in kelvin (~288 K), not degrees Celsius (~15).
    assert Earth.Tsbar == pytest.approx(288.0, rel=1e-6)
    assert Earth.Tsbar > 273.0


@pytest.mark.physics_invariant
def test_cpv_constant_positive():
    """Constant-mode heat capacity is finite and positive for every species.

    In constant mode ``cpv`` returns the per-species heat capacity from the
    coefficient table; a physical heat capacity must be finite and strictly
    positive for all supported volatiles. An unknown species falls through to
    zero (the guard path), which discriminates the positive-return contract
    for a real gas from the empty-table default.
    """
    for vol in ['H2O', 'CO2', 'N2', 'H2', 'CH4', 'O2']:
        cp = cpv(vol, 300.0, cp_mode='constant')
        assert np.isfinite(cp)
        assert cp > 0.0
    # Guard path: an unknown species returns exactly zero, not a positive cp.
    assert cpv('Xx_not_a_gas', 300.0, cp_mode='constant') == pytest.approx(0.0, abs=1e-12)
