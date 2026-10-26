"""Tests for src/janus/utils/height.py.

height.py provides the Newtonian surface-gravity law ``gravity`` and the
hydrostatic altitude integrator ``integrate_heights``. This file exercises:

* The inverse-square gravity law pinned at Earth's mass and radius.
* Hydrostatic integration of an isothermal column against the analytic
  scale-height integral H*ln(p_surface/p).
* The blow-up guard that flags a diverging integration and returns a
  bounded dummy grid.
* The local gravity profile left on the column, which must be the
  inverse-square law evaluated on the returned heights.
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
def test_surface_gravity_matches_earth_inverse_square(validation_plot):
    """Newtonian surface gravity g=GM/r^2 recovers Earth's ~9.82 m/s^2.

    ``gravity`` is the inverse-square Newtonian law. Evaluated at Earth's
    mass and mean radius it returns ~9.82 m/s^2 (just above the standard
    9.81 because the mean radius is smaller than the equatorial radius).
    Doubling the radius must quarter g, which discriminates an inverse-linear
    bug; doubling the mass must double g.
    """
    g_ref, g_tol, r_far = 9.82, 0.02, 2.0
    grav = height.gravity(M_EARTH, R_EARTH)
    assert grav == pytest.approx(g_ref, abs=g_tol)
    # Inverse-square guard: 2x radius -> g/4 exactly, not g/2 (inverse-linear).
    grav_double_r = height.gravity(M_EARTH, r_far * R_EARTH)
    assert grav_double_r == pytest.approx(grav / r_far**2, rel=1e-12)
    assert grav_double_r < grav / 3.0  # rules out the inverse-linear g/2 result
    # Mass linearity and positivity.
    assert height.gravity(2.0 * M_EARTH, R_EARTH) == pytest.approx(2.0 * grav, rel=1e-12)
    assert grav > 0.0

    def figure():
        import matplotlib.pyplot as plt

        radii = np.linspace(1.0, 3.0, 41)
        fig, ax = plt.subplots()
        ax.plot(
            radii, [height.gravity(M_EARTH, r * R_EARTH) for r in radii], label='code, gravity'
        )
        ax.errorbar(
            [1.0], [g_ref], yerr=g_tol, fmt='s', mfc='none', capsize=4, label=r'9.82 $\pm$ 0.02'
        )
        ax.plot([r_far], [grav / r_far**2], 's', mfc='none', ms=10, label='g/4')
        ax.plot([1.0, r_far], [grav, grav_double_r], 'o', label='code at R and 2R')
        ax.plot([r_far], [grav / r_far], 'x', label='g/2 (inverse-linear)')
        ax.set(xlabel=r'r [$R_\oplus$]', ylabel=r'g [m s$^{-2}$]')
        ax.legend(fontsize='small')
        fig.tight_layout()
        return fig

    validation_plot(figure, 'janus_height_gravity')


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_integrate_heights_recovers_isothermal_scale_height(validation_plot):
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

    def figure():
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.plot(scale_h * np.log(press.max() / press) / 1e3, press, label=r'$H \ln(p_s/p)$')
        keep = np.unique(np.r_[0:nlev:5, np.argmin(zcen), np.argmax(zcen)])
        height_km = (zcen - zcen.min()) / 1e3
        ax.plot(height_km[keep], atm.p[keep], 'o', ms=3, label='code, every 5th level')
        ax.set(xlabel='height above the lowest level [km]', ylabel='p [Pa]', yscale='log')
        ax.invert_yaxis()
        ax.set_title(
            f'span: code {z_span / 1e3:.2f} km, analytic {analytic_span / 1e3:.2f} km',
            fontsize='small',
        )
        ax.legend(fontsize='small')
        fig.tight_layout()
        return fig

    validation_plot(figure, 'janus_height_isothermal_column')


@pytest.mark.physics_invariant
def test_integrate_heights_flags_blowup_and_returns_bounded_dummy():
    """A runaway hydrostatic integration is flagged and returns bounded dummy z.

    When the hydrostatic step diverges (an unphysically hot column drives dz
    or z past 1e8 m), ``integrate_heights`` must set the height_error flag and
    replace the altitude grid with a bounded 0..1000 m linear ramp rather than
    propagating the blow-up. The local gravity profile must then be evaluated
    on that dummy grid. This is the routine's error contract.
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
    # Local gravity follows the dummy grid rather than the diverged heights or
    # a partially filled array from the interrupted loop: finite, bounded by
    # its values at 0 and 1000 m, and inverse-square in the returned heights.
    grav_s = height.gravity(M_EARTH, R_EARTH)
    assert len(atm.grav_z) == nlev
    assert np.all(np.isfinite(atm.grav_z))
    assert np.all(atm.grav_z <= grav_s * (1.0 + 1e-12))
    assert np.all(atm.grav_z >= height.gravity(M_EARTH, R_EARTH + 1000.0) * (1.0 - 1e-12))
    np.testing.assert_allclose(
        atm.grav_z, grav_s * R_EARTH**2 / (R_EARTH + zcen) ** 2, rtol=1e-12
    )


@pytest.mark.physics_invariant
def test_integrate_heights_writes_local_gravity_on_returned_heights():
    """Local gravity is stored on the same heights the integrator returns.

    ``integrate_heights`` evaluates g(z) = g_s R^2 / (R + z)^2 at each level as
    it steps up the column and must leave that profile on ``atm.grav_z``, in
    the same level order as the returned heights. A column seeded with a stale
    constant gravity must come back with the inverse-square profile, equal to
    g_s at the surface (z = 0, the edge case) and strictly decreasing upward.
    The isothermal 288 K N2 column from 1e5 to 1e3 Pa is about 40 km tall, so
    gravity at its top sits about 1.2% below g_s: far outside tolerance for a
    stale constant array, and the reversed profile fails the relation too.
    Each height step must also close the hydrostatic relation with the stored
    gravity, so the profile is the one the integration actually used.
    """
    nlev = 200
    press = np.logspace(3, 5, nlev)  # index 0 = top of atmosphere
    grav_s = height.gravity(M_EARTH, R_EARTH)
    atm = SimpleNamespace(
        p=press.copy(),
        tmp=np.full(nlev, 288.0),
        vol_list={'N2': 1.0},
        x_gas={'N2': np.ones(nlev)},
        grav_z=np.full(nlev, grav_s),  # stale profile left by an earlier solve
    )
    zcen, _, err = integrate_heights(atm, M_EARTH, R_EARTH)
    assert err is False

    expected = grav_s * R_EARTH**2 / (R_EARTH + zcen) ** 2
    assert len(atm.grav_z) == nlev
    np.testing.assert_allclose(atm.grav_z, expected, rtol=1e-12)
    # Edge case: the surface level (highest pressure, last index) has z = 0.
    assert zcen[-1] == pytest.approx(0.0, abs=1e-9)
    assert atm.grav_z[-1] == pytest.approx(grav_s, rel=1e-12)
    # Gravity falls strictly as height rises.
    order = np.argsort(zcen)
    assert np.all(np.diff(atm.grav_z[order]) < 0.0)
    # Stale-constant guard: the top level is well below surface gravity.
    assert atm.grav_z[0] < 0.99 * grav_s
    # Orientation guard: the level-reversed profile does not match the heights.
    assert not np.allclose(atm.grav_z[::-1], expected, rtol=1e-6, atol=0.0)
    # Hydrostatic closure per step, surface first: dz_n = R T dp / (M g_n p_n).
    # Stepping with constant g_s instead of g_n shifts dz by up to ~1.2%.
    p_up, z_up, g_up = press[::-1], zcen[::-1], atm.grav_z[::-1]
    dz_expected = (
        phys.R_gas
        * 288.0
        / (phys.molar_mass['N2'] * g_up[:-1] * p_up[:-1])
        * (p_up[:-1] - p_up[1:])
    )
    np.testing.assert_allclose(np.diff(z_up), dz_expected, rtol=1e-9)


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
