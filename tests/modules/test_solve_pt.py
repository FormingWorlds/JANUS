"""Tests for src/janus/modules/solve_pt.py.

solve_pt.py exposes the atmosphere-structure pipeline entries ``RadConvEqm``,
``MCPA`` and ``MCPA_CBL``. This file exercises the guard and fallback branches of
those entries:

* ``RadConvEqm`` with the deprecated dry-adiabat time-stepping branch enabled,
  and its empty-result fallback when the dry case is disabled.
* ``MCPA`` returning the moist-adiabat state directly.
* ``MCPA_CBL`` finding the surface temperature that closes the conductive-lid
  energy balance: the top-of-atmosphere and surface flux boundary conditions, the
  secant default second guess, the bracketing (brentq) method, the surface-ceiling
  clamp, the invalid-method error contract, and the non-convergence warning.

The moist-adiabat and dry-adiabat solvers are mocked at the narrowest scope with
plausible fluxes; the atmosphere objects passed to ``MCPA_CBL`` are real ``atmos``
instances built from synthetic band edges, never downloaded data. See
docs/How-to/test.md.
"""

from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

import janus.modules.solve_pt as solve_pt
from janus.utils.atmosphere_column import atmos

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Synthetic ascending band edges in nm; three bands is enough for a valid column.
BAND_EDGES = [1.0, 100.0, 1000.0, 10000.0]


def _flux_atm(net_surface, olr):
    """Small atmosphere-like object carrying only the fluxes the log lines read."""
    return SimpleNamespace(
        net_flux=np.array([net_surface, 0.5 * net_surface]),
        LW_flux_up=np.array([olr, 0.9 * olr]),
    )


def _make_atm_inp(t_surf=1000.0):
    """Build a real clear-sky atmos with synthetic bands for the MCPA_CBL solver."""
    return atmos(
        t_surf,
        1.0e5,
        1.0,
        6.371e6,
        5.972e24,
        BAND_EDGES,
        vol_mixing={'H2O': 0.3, 'CO2': 0.4, 'N2': 0.3},
        req_levels=15,
        trppT=250.0,
        minT=1.0,
        maxT=9000.0,
    )


def _fake_moist(atm_in, dirs, standalone, trppD, rscatter=False):
    """Mock moist-adiabat solve: attach plausible fluxes with a strong top/bottom contrast.

    The net flux runs +200 W/m^2 at the top to -200 W/m^2 at the base. Against the
    skin flux 200*(3000 - Ts), the top end balances at Ts = 2999 K and the bottom
    end at Ts = 3001 K, so the two boundary-condition roots are 2 K apart, far
    beyond the 0.1 K pin tolerance used below.
    """
    n = atm_in.nlev_save
    atm_in.net_flux = np.linspace(200.0, -200.0, n)  # W/m^2, sign flips top to bottom
    atm_in.LW_flux_up = np.linspace(280.0, 240.0, n)  # W/m^2 outgoing longwave
    return atm_in


def test_radconveqm_dry_branch_returns_both_states(caplog):
    """RadConvEqm with the dry branch on returns both the dry and moist states.

    With the dry-adiabat case enabled, ``RadConvEqm`` runs both the moist and the
    dry solve and returns the pair, logging their fluxes in standalone mode. Both
    solvers must be invoked exactly once and their returned states threaded back
    unchanged with finite fluxes.
    """
    atm_moist = _flux_atm(120.0, 250.0)
    atm_dry = _flux_atm(90.0, 240.0)
    with (
        patch('janus.modules.solve_pt.compute_moist_adiabat', return_value=atm_moist) as moist,
        patch('janus.modules.solve_pt.compute_dry_adiabat', return_value=atm_dry) as dry,
        caplog.at_level('INFO'),
    ):
        out_dry, out_moist = solve_pt.RadConvEqm(
            {}, {}, SimpleNamespace(), True, True, False, False
        )
    assert out_moist is atm_moist
    assert out_dry is atm_dry
    assert moist.call_count == 1 and dry.call_count == 1
    assert np.isfinite(out_dry.net_flux[0]) and np.isfinite(out_moist.net_flux[0])


def test_radconveqm_without_dry_returns_empty_mapping():
    """RadConvEqm with the dry case disabled returns an empty dry result.

    When the dry-adiabat case is turned off, the dry solver must not run and the
    dry return slot falls back to an empty mapping while the moist state is still
    returned. This is the degenerate branch of the entry point.
    """
    atm_moist = _flux_atm(120.0, 250.0)
    with (
        patch('janus.modules.solve_pt.compute_moist_adiabat', return_value=atm_moist),
        patch('janus.modules.solve_pt.compute_dry_adiabat') as dry,
    ):
        out_dry, out_moist = solve_pt.RadConvEqm(
            {}, {}, SimpleNamespace(), False, False, False, False
        )
    assert out_dry == {}
    assert out_moist is atm_moist
    assert dry.call_count == 0  # dry solver skipped entirely


def test_mcpa_returns_moist_adiabat_state():
    """MCPA returns the moist-adiabat state from a single solve.

    ``MCPA`` is the thin entry that prescribes a stratosphere and returns the
    moist-adiabat solution directly. It must call the moist solver once and return
    its result with finite fluxes.
    """
    atm_moist = _flux_atm(120.0, 250.0)
    with patch('janus.modules.solve_pt.compute_moist_adiabat', return_value=atm_moist) as moist:
        out = solve_pt.MCPA({}, SimpleNamespace(), False, True, False)
    assert out is atm_moist
    assert moist.call_count == 1
    assert np.isfinite(out.net_flux[0])


@pytest.mark.physics_invariant
def test_mcpa_cbl_secant_solves_toa_energy_balance():
    """MCPA_CBL (secant, TOA) finds a bounded surface temperature closing the lid balance.

    With the top-of-atmosphere boundary condition and the secant method, MCPA_CBL
    steps the surface temperature until the conductive-lid skin flux matches the
    top-of-column net flux. The default disabled surface guess triggers the second
    initial guess at 0.8*T_magma. The residual must balance against the top net
    flux (+200 W/m^2), giving Ts = 2999 K, which is 2 K from the surface-BC root
    (3001 K) and so discriminates the boundary-condition branch. The returned
    surface temperature lies inside [minT, maxT] and is strictly positive.
    """
    atm_inp = _make_atm_inp()
    with patch(
        'janus.modules.solve_pt.compute_moist_adiabat', side_effect=_fake_moist
    ) as moist:
        atm = solve_pt.MCPA_CBL({}, atm_inp, trppD=False, rscatter=False, atm_bc=0)
    # Skin flux 200*(3000 - Ts) balances the top net flux +200 W/m^2 at 2999 K; the
    # surface-BC root 3001 K is 2 K away, 20x the 0.1 K tolerance, so this pin
    # discriminates the top-of-atmosphere branch from the surface branch.
    assert atm.ts == pytest.approx(2999.0, abs=0.1)
    assert atm_inp.minT <= atm.ts <= atm_inp.maxT  # bounded to the allowed window
    assert atm.ts > 0.0
    # The residual balanced the skin flux against the TOP of the column: the
    # reconstructed skin flux at the solution matches net_flux[0], not net_flux[-1].
    skin = atm.skin_k / atm.skin_d * (atm.tmp_magma - atm.ts)
    assert skin == pytest.approx(atm.net_flux[0], abs=1.0)
    assert abs(skin - atm.net_flux[-1]) > 100.0  # not the bottom flux
    assert moist.call_count >= 2  # residual evaluations plus the final solve


@pytest.mark.physics_invariant
def test_mcpa_cbl_surface_boundary_condition_uses_bottom_flux():
    """MCPA_CBL with the surface boundary condition reads the bottom-of-column net flux.

    Selecting the surface boundary condition makes both the residual and the final
    diagnostic use the bottom-of-column net flux rather than the top. The residual
    balances against the bottom net flux (-200 W/m^2), giving Ts = 3001 K, which is
    2 K from the top-of-atmosphere root (2999 K); the reconstructed skin flux at the
    solution matches the bottom net flux and not the top, so the assertion depends on
    which end fed the residual.
    """
    atm_inp = _make_atm_inp()
    with patch('janus.modules.solve_pt.compute_moist_adiabat', side_effect=_fake_moist):
        atm = solve_pt.MCPA_CBL({}, atm_inp, trppD=False, rscatter=False, atm_bc=1)
    # Skin flux 200*(3000 - Ts) balances the bottom net flux -200 W/m^2 at 3001 K; the
    # TOA root 2999 K is 2 K away, 20x the 0.1 K tolerance.
    assert atm.ts == pytest.approx(3001.0, abs=0.1)
    # The residual balanced the skin flux against the BOTTOM of the column: the
    # reconstructed skin flux at the solution matches net_flux[-1], not net_flux[0].
    skin = atm.skin_k / atm.skin_d * (atm.tmp_magma - atm.ts)
    assert skin == pytest.approx(atm.net_flux[-1], abs=1.0)
    assert abs(skin - atm.net_flux[0]) > 100.0  # not the top flux
    assert atm_inp.minT <= atm.ts <= atm_inp.maxT


@pytest.mark.physics_invariant
def test_mcpa_cbl_surface_ceiling_clamps_solution():
    """MCPA_CBL clamps the solved surface temperature to the supplied ceiling.

    When a surface-temperature ceiling below the unconstrained solution is given,
    the returned surface temperature must be clamped to that ceiling rather than
    the ~3000 K balance point. The clamp is the boundedness guard on the solver
    output.
    """
    atm_inp = _make_atm_inp()
    with patch('janus.modules.solve_pt.compute_moist_adiabat', side_effect=_fake_moist):
        atm = solve_pt.MCPA_CBL(
            {}, atm_inp, trppD=False, rscatter=False, atm_bc=0, T_surf_max=1500.0
        )
    assert atm.ts == pytest.approx(1500.0, rel=1e-9)  # pinned at the ceiling
    assert atm.ts < 2000.0  # well below the ~3000 K unconstrained solution


@pytest.mark.physics_invariant
def test_mcpa_cbl_brentq_method_brackets_solution():
    """MCPA_CBL with the bracketing method finds the solution inside the bracket.

    The brentq method brackets the surface temperature between 800 K and the
    ceiling (here capped below the maximum), then narrows to the balance point.
    The returned surface temperature must lie inside that bracket and match the
    analytic root of the mock residual.
    """
    atm_inp = _make_atm_inp()
    with patch('janus.modules.solve_pt.compute_moist_adiabat', side_effect=_fake_moist):
        atm = solve_pt.MCPA_CBL(
            {}, atm_inp, trppD=False, rscatter=False, atm_bc=0, T_surf_max=5000.0, method=1
        )
    assert 800.0 < atm.ts < 5000.0  # inside the brentq bracket
    assert atm.ts == pytest.approx(2999.0, abs=0.1)


def test_mcpa_cbl_invalid_method_raises():
    """MCPA_CBL rejects an unknown root-finding method.

    Only the secant (0) and brentq (1) methods are defined; any other method
    selector must raise rather than silently returning a bogus surface
    temperature. This is the error contract of the solver dispatch.
    """
    atm_inp = _make_atm_inp()
    with patch('janus.modules.solve_pt.compute_moist_adiabat', side_effect=_fake_moist):
        with pytest.raises(Exception, match='Invalid solution method'):
            solve_pt.MCPA_CBL({}, atm_inp, trppD=False, rscatter=False, method=2)


def test_mcpa_cbl_reports_non_convergence(caplog):
    """MCPA_CBL warns and still returns a state when the root solve does not converge.

    If the root finder reports non-convergence, the solver must surface a warning
    rather than claim success, and still return the atmosphere state evaluated at
    the (clamped) best-estimate surface temperature. The non-converged flag is the
    error-contract path.
    """
    atm_inp = _make_atm_inp()
    fake_result = SimpleNamespace(root=1234.0, converged=False)
    with (
        patch('janus.modules.solve_pt.compute_moist_adiabat', side_effect=_fake_moist),
        patch('janus.modules.solve_pt.optimise.root_scalar', return_value=fake_result),
        caplog.at_level('WARNING'),
    ):
        atm = solve_pt.MCPA_CBL({}, atm_inp, trppD=False, rscatter=False, atm_bc=0)
    assert any('Did not find solution' in record.message for record in caplog.records)
    assert atm.ts == pytest.approx(1234.0, rel=1e-9)  # best estimate, inside [minT, maxT]
