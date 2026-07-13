"""Tests for src/janus/modules/dry_adiabat_timestep.py.

dry_adiabat_timestep.py provides ``compute_dry_adiabat``, the deprecated
radiative time-stepping loop: each step it computes radiative heating, optionally
does a surface energy balance and a pure-steam adjustment, applies the heating,
runs the dry convective adjustment, and clamps temperatures to a floor. This
file exercises:

* Positive radiative heating raising temperatures and the dry adjustment running
  its full inner sweep every step.
* Strong radiative cooling being clamped at the 10 K temperature floor.
* The surface time-stepping branch heating the surface and calling the boundary
  scheme.
* The pure-steam adjustment adding its moist tendency each step.
* The radiative-solver failure path breaking the loop without corrupting the
  profile.

The heavy internals (general adiabat, dry-adiabat setup, SOCRATES, dry and moist
adjustment, cloud, surface boundary) are mocked at the narrowest scope with
physically plausible heating rates. See docs/How-to/test.md.
"""

import copy  # noqa: F401  (compute_dry_adiabat deep-copies the input atmosphere)
from contextlib import ExitStack
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

import janus.modules.dry_adiabat_timestep as dat
from janus.modules.dry_adiabat_timestep import compute_dry_adiabat

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

N = 5  # number of vertical levels in the synthetic column
CONV_STEPS = 30  # dry convective adjustment passes per radiation step (source constant)


def _install(stack, atm_ga, atm_dry, rad_side, moist_side=None, boundary_side=None):
    """Patch every heavy internal of the time-stepping loop and return the mocks."""

    def _identity(atm):
        return atm

    def _default_moist(atm, tau):
        return np.zeros(N)

    def _default_boundary(idx, tmp, ts, d_atm, d_surf, mix_a, mix_s):
        return d_atm, d_surf

    stack.enter_context(patch.object(dat.ga, 'general_adiabat', return_value=atm_ga))
    stack.enter_context(patch.object(dat, 'dry_adiabat_atm', return_value=atm_dry))
    rad = stack.enter_context(patch.object(dat.socrates, 'radCompSoc', side_effect=rad_side))
    dryadj = stack.enter_context(patch.object(dat, 'DryAdj', side_effect=_identity))
    cloud = stack.enter_context(patch.object(dat, 'simple_cloud', side_effect=_identity))
    moist = stack.enter_context(
        patch.object(dat, 'moist_adj', side_effect=moist_side or _default_moist)
    )
    boundary = stack.enter_context(
        patch.object(
            dat, 'simple_boundary_tend', side_effect=boundary_side or _default_boundary
        )
    )
    return dict(rad=rad, dryadj=dryadj, cloud=cloud, moist=moist, boundary=boundary)


@pytest.mark.physics_invariant
def test_positive_heating_raises_temperatures_and_applies_dry_adjustment():
    """Positive net heating warms every level and the dry adjustment runs each step.

    With a uniform positive radiative heating rate, each step raises the profile
    by net_heating*dt (below the 20 K per-step cap), so the final temperatures
    exceed the initial ones at every level and the array length is preserved. The
    dry convective adjustment must run its full inner sweep (conv_steps passes)
    per radiation step; the accumulated heating is pinned to the number of steps,
    discriminating a loop that drops or double-applies the tendency.
    """
    atm_dry = SimpleNamespace(
        tmp=np.array([200.0, 220.0, 240.0, 260.0, 280.0]), ts=300.0, dt=0.5
    )
    initial = atm_dry.tmp.copy()
    atm_ga = SimpleNamespace(tmp=np.zeros(N), do_cloud=True)

    def rad(atm, dirs, recalc=False, rscatter=False):
        atm.net_heating = np.full(N, 4.0)  # K/day, uniform radiative heating
        atm.LW_flux_up = np.full(N, 300.0)  # constant OLR so the loop converges
        return atm

    with ExitStack() as stack:
        mocks = _install(stack, atm_ga, atm_dry, rad)
        out = compute_dry_adiabat(SimpleNamespace(), dirs={}, standalone=False)

    step = 4.0 * 0.5  # net_heating * dt, uniform and below the per-step cap
    assert len(out.tmp) == N  # column length preserved through the loop
    assert np.all(out.tmp > initial)  # warmed in the direction of positive heating
    np.testing.assert_allclose(out.tmp - initial, mocks['rad'].call_count * step, rtol=1e-9)
    # The dry adjustment runs conv_steps passes for each radiation step.
    assert mocks['dryadj'].call_count == mocks['rad'].call_count * CONV_STEPS
    # The break guard requires i>5, so the loop advances past the sixth step.
    assert mocks['rad'].call_count >= 6


@pytest.mark.physics_invariant
def test_strong_cooling_is_clamped_at_the_temperature_floor():
    """Strong radiative cooling is clamped so no level falls below the 10 K floor.

    A large negative heating rate would drive temperatures negative, but the loop
    clamps any level below the 10 K floor back up to it. Starting from a cold
    profile, the final temperatures must all sit at or above the floor, at least
    one level must be pinned exactly to the floor (the clamp fired), and every
    level must be cooler than it started. This is the positivity guard of the
    time-stepping loop.
    """
    atm_dry = SimpleNamespace(tmp=np.array([15.0, 20.0, 25.0, 30.0, 35.0]), ts=300.0, dt=0.5)
    initial = atm_dry.tmp.copy()
    atm_ga = SimpleNamespace(tmp=np.zeros(N), do_cloud=True)

    def rad(atm, dirs, recalc=False, rscatter=False):
        atm.net_heating = np.full(N, -100.0)  # strong radiative cooling
        atm.LW_flux_up = np.full(N, 300.0)
        return atm

    with ExitStack() as stack:
        _install(stack, atm_ga, atm_dry, rad)
        out = compute_dry_adiabat(SimpleNamespace(), dirs={}, standalone=False)

    floor = 10.0
    assert np.all(out.tmp >= floor)  # floor enforced everywhere
    assert np.any(np.isclose(out.tmp, floor))  # at least one level pinned to the floor
    assert np.all(out.tmp < initial)  # the whole column cooled


@pytest.mark.physics_invariant
def test_surface_timestepping_heats_surface_and_runs_boundary_scheme():
    """The surface time-stepping branch warms the surface under a net downward flux.

    With surface time-stepping enabled and a net downward surface flux (absorbed
    shortwave and downwelling longwave exceeding the upwelling terms), the surface
    temperature must rise from its initial value, and the boundary tendency scheme
    must be invoked once per radiation step. The column length is preserved.
    """
    atm_dry = SimpleNamespace(
        tmp=np.array([200.0, 220.0, 240.0, 260.0, 280.0]), ts=300.0, dt=0.5
    )
    ts_initial = atm_dry.ts
    atm_ga = SimpleNamespace(tmp=np.zeros(N), do_cloud=True)

    def rad(atm, dirs, recalc=False, rscatter=False):
        atm.net_heating = np.full(N, 2.0)
        # Edge arrays span the N+1 level interfaces; the surface term is index N.
        atm.LW_flux_up = np.full(N + 1, 100.0)  # OLR and upwelling surface LW
        atm.SW_flux_down = np.full(N + 1, 300.0)
        atm.LW_flux_down = np.full(N + 1, 200.0)
        atm.SW_flux_up = np.full(N + 1, 50.0)
        return atm

    with ExitStack() as stack:
        mocks = _install(stack, atm_ga, atm_dry, rad)
        out = compute_dry_adiabat(SimpleNamespace(), dirs={}, standalone=True, surf_dt=True)

    # Net surface flux 300+200-100-50 = 350 W/m^2 downward warms the surface.
    assert out.ts > ts_initial
    assert mocks['boundary'].call_count == mocks['rad'].call_count
    assert len(out.tmp) == N


@pytest.mark.physics_invariant
def test_pure_steam_adjustment_adds_moist_tendency_each_step():
    """The pure-steam adjustment contributes its moist tendency on every step.

    When pure-steam adjustment is enabled, each step adds the moist-adjustment
    tendency on top of the radiative heating. With a uniform radiative heating and
    a small uniform moist tendency, the accumulated warming per step is their sum,
    pinned against the step count; the moist routine must be called once per
    radiation step. This discriminates a loop that skips the moist contribution.
    """
    atm_dry = SimpleNamespace(
        tmp=np.array([200.0, 220.0, 240.0, 260.0, 280.0]), ts=300.0, dt=0.5
    )
    initial = atm_dry.tmp.copy()
    atm_ga = SimpleNamespace(tmp=np.zeros(N), do_cloud=True)

    def rad(atm, dirs, recalc=False, rscatter=False):
        atm.net_heating = np.full(N, 4.0)
        atm.LW_flux_up = np.full(N, 300.0)
        return atm

    def moist(atm, tau):
        return np.full(N, 0.5)  # small positive moist-adjustment tendency

    with ExitStack() as stack:
        mocks = _install(stack, atm_ga, atm_dry, rad, moist_side=moist)
        out = compute_dry_adiabat(
            SimpleNamespace(), dirs={}, standalone=False, pure_steam_adj=True
        )

    step = 4.0 * 0.5 + 0.5  # radiative heating plus the moist tendency per step
    assert mocks['moist'].call_count == mocks['rad'].call_count
    np.testing.assert_allclose(out.tmp - initial, mocks['rad'].call_count * step, rtol=1e-9)
    assert np.all(out.tmp > initial)


def test_radiation_failure_breaks_loop_without_corrupting_profile():
    """A radiative-solver failure breaks the loop and leaves the profile untouched.

    The loop wraps each step in a guard: if the radiative transfer raises, the
    step is abandoned and the loop breaks. Because the failure happens before the
    heating is applied, the returned profile must equal the input profile, and the
    solver must have been called exactly once. This is the error-contract path for
    a solver that cannot run.
    """
    tmp0 = np.array([200.0, 220.0, 240.0, 260.0, 280.0])
    atm_dry = SimpleNamespace(tmp=tmp0.copy(), ts=300.0, dt=0.5)
    atm_ga = SimpleNamespace(tmp=np.zeros(N), do_cloud=True)

    def rad(atm, dirs, recalc=False, rscatter=False):
        raise RuntimeError('radiative transfer failed')

    with ExitStack() as stack:
        mocks = _install(stack, atm_ga, atm_dry, rad)
        out = compute_dry_adiabat(SimpleNamespace(), dirs={}, standalone=True)

    assert mocks['rad'].call_count == 1  # failed on the first radiation call
    np.testing.assert_array_equal(out.tmp, tmp0)  # profile untouched by the failed step


def test_clear_sky_call_performs_no_radiative_stepping():
    """Document a known defect: a clear-sky call runs no radiation and no stepping.

    In this deprecated module the radiative increment is computed only inside the
    cloud branch, so a clear-sky atmosphere (clouds off, no surface stepping) never
    assigns it. Applying the increment then raises an internal NameError that the
    bare except swallows, breaking the loop on the first pass. Radiation is never
    run and the temperature profile is returned unchanged. This records the current
    behavior, not a designed contract: radiation should not depend on the cloud
    flag.
    """
    tmp0 = np.array([200.0, 220.0, 240.0, 260.0, 280.0])
    atm_dry = SimpleNamespace(tmp=tmp0.copy(), ts=300.0, dt=0.5)
    atm_ga = SimpleNamespace(tmp=np.zeros(N), do_cloud=False)  # clear sky

    def rad(atm, dirs, recalc=False, rscatter=False):
        atm.net_heating = np.full(N, 5.0)
        atm.LW_flux_up = np.full(N, 300.0)
        return atm

    with ExitStack() as stack:
        mocks = _install(stack, atm_ga, atm_dry, rad)
        out = compute_dry_adiabat(SimpleNamespace(), dirs={}, standalone=False)

    # The radiative solve sits inside the cloud branch, so a clear-sky call skips it.
    assert mocks['rad'].call_count == 0
    # The dry adjustment is never reached; the profile is returned untouched.
    assert mocks['dryadj'].call_count == 0
    np.testing.assert_array_equal(out.tmp, tmp0)
