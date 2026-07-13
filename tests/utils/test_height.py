"""Tests for src/janus/utils/height.py.

height.py provides the Newtonian surface-gravity law ``gravity`` and the
hydrostatic altitude integrator ``integrate_heights``. This file exercises:

* The inverse-square gravity law pinned at Earth's mass and radius.
* Hydrostatic integration of an isothermal column against the analytic
  scale-height integral H*ln(p_surface/p).
* The blow-up guard that flags a diverging integration and returns a
  bounded dummy grid.
* The domain guard that rejects a non-positive radius rather than dividing
  by zero.

Invariant families exercised: positivity/boundedness, monotonicity, and
pinned-value-with-discrimination-guard. See docs/How-to/test.md.
"""

import math
from types import SimpleNamespace

import numpy as np
import pytest

import janus.utils.height as height
import janus.utils.phys as phys
from janus.utils.height import integrate_heights

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Earth reference values (mean radius, so g sits slightly above 9.81).
M_EARTH = 5.972e24
R_EARTH = 6.371e6


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_surface_gravity_matches_earth_inverse_square():
    """Newtonian surface gravity g=GM/r^2 recovers Earth's ~9.82 m/s^2.

    ``gravity`` is the inverse-square Newtonian law. Evaluated at Earth's
    mass and mean radius it returns ~9.82 m/s^2 (just above the standard
    9.81 because the mean radius is smaller than the equatorial radius).
    Doubling the radius must quarter g, which discriminates an inverse-linear
    bug; doubling the mass must double g.
    """
    grav = height.gravity(M_EARTH, R_EARTH)
    assert grav == pytest.approx(9.82, abs=0.02)
    # Inverse-square guard: 2x radius -> g/4 exactly, not g/2 (inverse-linear).
    grav_double_r = height.gravity(M_EARTH, 2.0 * R_EARTH)
    assert grav_double_r == pytest.approx(grav / 4.0, rel=1e-12)
    assert grav_double_r < grav / 3.0  # rules out the inverse-linear g/2 result
    # Mass linearity and positivity.
    assert height.gravity(2.0 * M_EARTH, R_EARTH) == pytest.approx(2.0 * grav, rel=1e-12)
    assert grav > 0.0


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_integrate_heights_recovers_isothermal_scale_height():
    """Hydrostatic integration of an isothermal N2 column matches H*ln(p_s/p).

    For an isothermal ideal-gas column the geopotential height grows as
    H*ln(p_surface/p) with scale height H=R_gas*T/(M*g) (analytic
    hydrostatic limit). ``integrate_heights`` must reproduce this: altitude
    increases monotonically as pressure falls, the total column height
    matches the scale-height integral to a few percent on a fine grid, and
    the input pressure array is restored (the routine reverses internally).
    """
    temp = 288.0
    nlev = 200
    # Input convention: index 0 = top of atmosphere (low p), index -1 = surface.
    press = np.logspace(3, 5, nlev)
    atm = SimpleNamespace(
        p=press.copy(),
        tmp=np.full(nlev, temp),
        vol_list={'N2': 1.0},
        x_gas={'N2': np.ones(nlev)},
    )
    p_input = atm.p.copy()
    zcen, zedge, err = integrate_heights(atm, M_EARTH, R_EARTH)

    assert err is False
    assert len(zcen) == nlev
    assert len(zedge) == nlev + 1
    # Altitude and pressure are anti-monotone: sorted by ascending pressure,
    # altitude must be strictly descending.
    order = np.argsort(atm.p)
    assert np.all(np.diff(zcen[order]) < 0.0)

    grav = height.gravity(M_EARTH, R_EARTH)
    scale_h = phys.R_gas * temp / (phys.molar_mass['N2'] * grav)
    analytic_span = scale_h * np.log(press.max() / press.min())
    z_span = zcen.max() - zcen.min()
    # Rectangle-rule hydrostatic sum matches the analytic log integral to a
    # few percent on this grid; a missing 1/p factor or wrong molar mass would
    # throw the scale height far outside 5%.
    assert z_span == pytest.approx(analytic_span, rel=0.05)
    # Autonomy: the input pressure array is unchanged after the internal flip.
    np.testing.assert_allclose(atm.p, p_input, rtol=0, atol=0)


@pytest.mark.physics_invariant
def test_integrate_heights_flags_blowup_and_returns_bounded_dummy():
    """A runaway hydrostatic integration is flagged and returns bounded dummy z.

    When the hydrostatic step diverges (an unphysically hot column drives dz
    or z past 1e8 m), ``integrate_heights`` must set the height_error flag and
    replace the altitude grid with a bounded 0..1000 m linear ramp rather than
    propagating the blow-up. This is the routine's error contract.
    """
    nlev = 20
    press = np.logspace(3, 5, nlev)
    atm = SimpleNamespace(
        p=press.copy(),
        tmp=np.full(nlev, 1.0e12),  # absurd temperature forces divergence
        vol_list={'N2': 1.0},
        x_gas={'N2': np.ones(nlev)},
    )
    zcen, zedge, err = integrate_heights(atm, M_EARTH, R_EARTH)
    assert err is True
    # Dummy grid is the bounded linear ramp, not a runaway value.
    assert zcen.min() >= 0.0
    assert zcen.max() == pytest.approx(1000.0, rel=1e-12)
    assert len(zcen) == nlev


@pytest.mark.physics_invariant
def test_gravity_rejects_nonpositive_radius():
    """Surface gravity raises on a non-positive radius instead of dividing by zero.

    ``gravity`` evaluates GM/r^2, singular at r=0 and unphysical for r<0. The
    routine guards the domain and raises ValueError rather than returning an
    infinite or negative gravity; this is the error contract. Just inside the
    valid domain gravity stays finite and positive, and shrinking the radius
    tenfold multiplies gravity by exactly 100 (inverse-square growth), the
    physical edge that discriminates a missing or inverse-linear guard.
    """
    with pytest.raises(ValueError):
        height.gravity(M_EARTH, 0.0)
    with pytest.raises(ValueError):
        height.gravity(M_EARTH, -1.0e6)
    # Just inside the valid domain: finite, positive, and larger than at the
    # surface because the radius is smaller.
    g_close = height.gravity(M_EARTH, R_EARTH / 10.0)
    assert math.isfinite(g_close)
    assert g_close > 0.0
    assert g_close == pytest.approx(height.gravity(M_EARTH, R_EARTH) * 100.0, rel=1e-9)
