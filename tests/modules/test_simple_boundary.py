"""Tests for src/janus/modules/simple_boundary.py.

simple_boundary.py provides ``simple_boundary_tend``, the surface-atmosphere
thermal exchange used by the standalone timestepping driver. It relaxes the
bottom atmospheric level toward the surface temperature and removes the
equal-and-opposite tendency from the surface. This file exercises:

* The antisymmetry of the exchange: with equal coefficients the atmosphere
  gains exactly what the surface loses (energy is exchanged, not created).
* The relaxation direction, magnitude, and coefficient scaling, plus the
  zero-contrast equilibrium in which no exchange occurs.

Invariant families exercised: conservation (energy-exchange antisymmetry),
monotonicity/sign of the relaxation, and pinned-value-with-discrimination-guard.
See docs/How-to/test.md.
"""

import numpy as np
import pytest

from janus.modules.simple_boundary import simple_boundary_tend

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.physics_invariant
def test_surface_exchange_is_antisymmetric_and_relaxes_toward_surface():
    """The surface-atmosphere exchange conserves energy and relaxes the base layer.

    With equal atmosphere and surface exchange coefficients the tendency added
    to the bottom atmospheric level equals the tendency removed from the surface
    (energy is exchanged, not created). The base-layer tendency is the pure
    relaxation (ts - T_base)/mix_coeff: a hotter surface warms the base layer
    with no additive offset, so the tendency stays of order the small relaxation
    rate rather than of order the surface temperature.
    """
    npz = -1
    n = 5
    t_base = 250.0
    ts = 300.0
    mix = 1.0e6
    pt = np.full(n, t_base)

    t_dt_out, ts_dt_out = simple_boundary_tend(
        npz, pt, ts, np.zeros(n), 0.0, mix_coeff_atmos=mix, mix_coeff_surf=mix
    )

    expected = (ts - t_base) / mix  # pure relaxation term
    # Base-layer heating equals the closed-form relaxation, no spurious offset.
    assert t_dt_out[npz] == pytest.approx(expected, rel=1e-12)
    # Antisymmetry: equal coefficients mean the surface loses what the air gains.
    assert ts_dt_out == pytest.approx(-t_dt_out[npz], rel=1e-12)
    # Sign: a hotter surface warms the base layer.
    assert t_dt_out[npz] > 0.0
    # Discrimination: an additive surface temperature would push the tendency to
    # order ts (hundreds of K), dwarfing the ~5e-5 relaxation term.
    assert t_dt_out[npz] < 1.0e-3
    assert abs((ts + expected) - t_dt_out[npz]) > 1.0


@pytest.mark.physics_invariant
def test_surface_exchange_reverses_sign_and_scales_with_coefficient():
    """A colder surface cools the base layer, and larger coefficients slow it.

    Reversing the temperature contrast (surface colder than the base layer)
    flips the sign of both tendencies, and doubling the exchange coefficient
    halves the relaxation rate. The zero-contrast case, surface at the
    base-layer temperature, produces no exchange at all: the equilibrium edge.
    """
    npz = -1
    n = 4
    t_base = 320.0
    ts = 280.0  # surface colder than the base layer

    pt = np.full(n, t_base)
    t_dt_out, ts_dt_out = simple_boundary_tend(
        npz, pt, ts, np.zeros(n), 0.0, mix_coeff_atmos=1.0e6, mix_coeff_surf=1.0e6
    )
    # Colder surface: the base layer cools and the surface warms.
    assert t_dt_out[npz] < 0.0
    assert ts_dt_out > 0.0

    # Coefficient scaling: doubling mix_coeff halves the tendency magnitude.
    t_dt_slow, _ = simple_boundary_tend(
        npz, pt, ts, np.zeros(n), 0.0, mix_coeff_atmos=2.0e6, mix_coeff_surf=2.0e6
    )
    assert t_dt_slow[npz] == pytest.approx(t_dt_out[npz] / 2.0, rel=1e-12)

    # Equilibrium edge case: no temperature contrast means no exchange.
    t_dt_eq, ts_dt_eq = simple_boundary_tend(npz, np.full(n, ts), ts, np.zeros(n), 0.0)
    assert t_dt_eq[npz] == pytest.approx(0.0, abs=1e-15)
    assert ts_dt_eq == pytest.approx(0.0, abs=1e-15)
