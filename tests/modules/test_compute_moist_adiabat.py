"""Tests for src/janus/modules/compute_moist_adiabat.py.

compute_moist_adiabat.py orchestrates one moist-adiabat atmosphere solve: it
builds the general adiabat, sets up clouds when requested, runs the radiative
transfer, and (when a tropopause is active) prescribes a stratosphere and
recomputes the fluxes. This file exercises:

* The cloudy path with a dynamically calculated tropopause, where the cloud is
  set up both before radiation and again after the stratosphere is prescribed,
  and the radiative transfer runs twice (initial and recalculation).
* The clear-sky limit with no active tropopause, where the stratosphere step is
  skipped and radiation runs exactly once.

The heavy internals (general adiabat, relative humidity, cloud, SOCRATES,
tropopause, stratosphere) are mocked at the narrowest scope with physically
plausible fluxes; this file checks the control flow and the flux threading, not
the physics of those internals. See docs/How-to/test.md.
"""

from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

import janus.modules.compute_moist_adiabat as cma
from janus.modules.compute_moist_adiabat import compute_moist_adiabat

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def _identity(atm, *args, **kwargs):
    """Return the atmosphere unchanged (stand-in for a structure-only step)."""
    return atm


def _fake_rad(atm, dirs, recalc=False, rscatter=False):
    """Stand-in for SOCRATES: attach physically plausible fluxes and return."""
    atm.net_flux = np.array([120.0, 60.0, -10.0])  # W/m^2, net down-positive
    atm.LW_flux_up = np.array([250.0, 240.0, 230.0])  # W/m^2, outgoing longwave
    return atm


def _make_atm(do_cloud, trppT, minT):
    return SimpleNamespace(do_cloud=do_cloud, trppT=trppT, minT=minT)


def test_cloudy_column_with_tropopause_runs_radiation_twice_and_reclouds():
    """A cloudy column with a dynamic tropopause clouds twice and recomputes fluxes.

    With clouds enabled and a dynamically calculated tropopause, the routine sets
    up the cloud before the first radiative solve and again after the stratosphere
    is prescribed, then runs SOCRATES a second time to recompute the fluxes. This
    verifies both cloud set-ups fire and that the second radiation call is a
    recalculation, discriminating a path that skips the post-stratosphere update.
    """
    atm = _make_atm(do_cloud=True, trppT=250.0, minT=180.0)
    rh_profile = np.array([0.4, 0.5, 0.6])
    with (
        patch.object(cma.ga, 'general_adiabat', side_effect=_identity),
        patch.object(cma, 'compute_Rh', return_value=rh_profile),
        patch.object(cma, 'simple_cloud', side_effect=_identity) as cloud,
        patch.object(cma.socrates, 'radCompSoc', side_effect=_fake_rad) as rad,
        patch.object(cma, 'find_tropopause', side_effect=_identity) as trop,
        patch.object(cma, 'set_stratosphere', side_effect=_identity) as strat,
    ):
        out = compute_moist_adiabat(atm, dirs={}, standalone=True, trppD=True, rscatter=False)

    # Radiation runs twice: the initial solve and the post-stratosphere recompute.
    assert rad.call_count == 2
    assert rad.call_args_list[0].kwargs['recalc'] is False
    assert rad.call_args_list[1].kwargs['recalc'] is True
    # The cloud is built before radiation and again after the stratosphere step.
    assert cloud.call_count == 2
    assert trop.call_count == 1 and strat.call_count == 1
    # The relative-humidity profile is stored on the returned atmosphere.
    np.testing.assert_array_equal(out.rh, rh_profile)
    # Plausible fluxes are threaded through to the caller.
    assert out.net_flux[0] == pytest.approx(120.0, rel=1e-12)
    assert np.all(np.isfinite(out.net_flux))


def test_clear_sky_without_tropopause_skips_stratosphere_and_radiates_once():
    """A clear-sky column with no active tropopause radiates once and skips the stratosphere.

    When clouds are off, the tropopause is not calculated dynamically, and the
    fixed tropopause temperature does not exceed the floor, the stratosphere step
    must be skipped entirely and radiation must run exactly once. This is the
    degenerate limit of the solve: no cloud set-up and no stratosphere
    recalculation. The tropopause temperature equal to the floor is the boundary
    input that keeps the stratosphere branch closed.
    """
    atm = _make_atm(do_cloud=False, trppT=180.0, minT=180.0)  # trppT == minT, not above
    with (
        patch.object(cma.ga, 'general_adiabat', side_effect=_identity),
        patch.object(cma, 'compute_Rh', return_value=np.array([0.3, 0.3, 0.3])),
        patch.object(cma, 'simple_cloud', side_effect=_identity) as cloud,
        patch.object(cma.socrates, 'radCompSoc', side_effect=_fake_rad) as rad,
        patch.object(cma, 'find_tropopause', side_effect=_identity) as trop,
        patch.object(cma, 'set_stratosphere', side_effect=_identity) as strat,
    ):
        out = compute_moist_adiabat(atm, dirs={}, standalone=False, trppD=False, rscatter=False)

    # Exactly one radiative solve; the stratosphere branch never fires.
    assert rad.call_count == 1
    assert trop.call_count == 0 and strat.call_count == 0
    # No cloud is built in the clear-sky case.
    assert cloud.call_count == 0
    # The single solve still returns the plausible fluxes.
    assert out.net_flux[-1] == pytest.approx(-10.0, rel=1e-12)
    assert out.LW_flux_up[0] > 0.0


def test_fixed_tropopause_enters_stratosphere_without_cloud_or_standalone_logging():
    """A fixed tropopause above the floor enters the stratosphere with clouds and logging off.

    The stratosphere branch fires whenever the fixed tropopause temperature exceeds
    the floor, even without dynamic tropopause calculation. With clouds off and
    standalone off, the routine must still find the tropopause, prescribe the
    stratosphere, and recompute the fluxes, but skip both the cloud set-up and the
    standalone flux logging. This exercises the cloud-off and quiet siblings of the
    cloudy standalone path. The tropopause temperature just above the floor is the
    boundary input that opens the branch.
    """
    atm = _make_atm(do_cloud=False, trppT=250.0, minT=180.0)  # trppT > minT opens the branch
    with (
        patch.object(cma.ga, 'general_adiabat', side_effect=_identity),
        patch.object(cma, 'compute_Rh', return_value=np.array([0.3, 0.3, 0.3])),
        patch.object(cma, 'simple_cloud', side_effect=_identity) as cloud,
        patch.object(cma.socrates, 'radCompSoc', side_effect=_fake_rad) as rad,
        patch.object(cma, 'find_tropopause', side_effect=_identity) as trop,
        patch.object(cma, 'set_stratosphere', side_effect=_identity) as strat,
    ):
        out = compute_moist_adiabat(atm, dirs={}, standalone=False, trppD=False, rscatter=False)

    # The stratosphere branch runs on the fixed tropopause, recomputing the fluxes.
    assert rad.call_count == 2
    assert trop.call_count == 1 and strat.call_count == 1
    # Clouds are off, so no cloud set-up happens on either radiation pass.
    assert cloud.call_count == 0
    assert out.net_flux[-1] == pytest.approx(-10.0, rel=1e-12)
