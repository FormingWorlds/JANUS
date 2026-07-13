"""Tests for src/janus/utils/cp_funcs.py.

cp_funcs.py provides the per-species molar heat capacity ``cpv`` with two
branches (a constant value and a temperature-dependent NIST Shomate
polynomial) and the condensate heat capacity ``cp_cond``. This file exercises:

* The Shomate branch for water vapour pinned against the NIST-JANAF value
  at 500 K.
* The unknown-species guard return (cp=0) in both branches.
* Positivity of a known species' heat capacity and rough agreement between
  the two branches.
* The liquid-hydrogen condensate ``cp_cond`` staying positive outside its
  cubic fit window through the temperature clamp.

Invariant families exercised: positivity/boundedness and
pinned-value-with-discrimination-guard. See docs/How-to/test.md.
"""

import pytest

import janus.utils.cp_funcs as cp_funcs

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_water_heat_capacity_matches_nist_shomate_at_500K():
    """cpv(H2O) reproduces the NIST Shomate water-vapour value at 500 K.

    The T-dependent branch evaluates the NIST Shomate polynomial. At 500 K
    the tabulated molar heat capacity of gaseous water is 35.22 J/mol/K
    (NIST-JANAF; Chase 1998); the implementation returns 35.218. Heat
    capacity must be positive and rise with temperature for water vapour as
    vibrational modes activate, which discriminates the Shomate branch from
    the flat constant branch.
    """
    cp500 = cp_funcs.cpv('H2O', 500.0, 'T-dependent')
    assert cp500 == pytest.approx(35.22, rel=2e-3)  # NIST Shomate H2O(g), 500 K
    # Positivity guard.
    assert cp500 > 0.0
    # Scale guard: order 3.5e1 J/mol/K, not 3.5e0 (per-gram slip) or 3.5e2.
    assert 25.0 < cp500 < 55.0
    # Temperature dependence: cp of water vapour rises with T.
    cp1500 = cp_funcs.cpv('H2O', 1500.0, 'T-dependent')
    assert cp1500 > cp500


@pytest.mark.physics_invariant
def test_unknown_species_returns_zero_and_known_species_positive():
    """Unknown volatiles fall through to cp=0; known species stay positive.

    Both branches look up per-gas coefficients and return 0.0 for a species
    with no table entry (the guard path / error contract). A known species
    returns a positive heat capacity in both branches, and the two branches
    agree to order of magnitude for CO2.
    """
    assert cp_funcs.cpv('Xx_not_a_gas', 500.0, 'T-dependent') == pytest.approx(0.0, abs=1e-12)
    assert cp_funcs.cpv('Xx_not_a_gas', 500.0, 'constant') == pytest.approx(0.0, abs=1e-12)
    cp_const = cp_funcs.cpv('CO2', 500.0, 'constant')
    cp_tdep = cp_funcs.cpv('CO2', 500.0, 'T-dependent')
    assert cp_const > 0.0
    assert cp_tdep > 0.0
    # Both branches land in the same 30-55 J/mol/K band for CO2, discriminating
    # a branch that silently returned the wrong unit scale.
    assert cp_tdep == pytest.approx(cp_const, rel=0.35)


@pytest.mark.physics_invariant
def test_condensate_hydrogen_heat_capacity_stays_positive_outside_fit_window():
    """Condensate cp for hydrogen stays positive beyond its liquid fit window.

    ``cp_cond`` fits liquid-hydrogen heat capacity with a cubic valid only
    between the triple point (13.8 K) and the critical point (33 K). Outside
    that range the bare cubic turns negative, which is unphysical; the routine
    clamps the evaluation temperature into the window so the returned molar heat
    capacity stays positive. The clamp is the error contract; a temperature far
    above the critical point is the edge case, and the raw cubic there is the
    discriminating wrong value.
    """
    # Far above the critical point the clamp pins the 33 K endpoint value.
    cp_hot = cp_funcs.cp_cond('H2', 300.0)
    cp_at_ceiling = cp_funcs.cp_cond('H2', 33.0)
    assert cp_hot > 0.0
    assert cp_hot == pytest.approx(cp_at_ceiling, rel=1e-12)
    # Below the triple point the clamp holds the 13.8 K endpoint, also positive.
    cp_cold = cp_funcs.cp_cond('H2', 5.0)
    cp_at_floor = cp_funcs.cp_cond('H2', 13.8)
    assert cp_cold > 0.0
    assert cp_cold == pytest.approx(cp_at_floor, rel=1e-12)
    # Discrimination: the unclamped cubic evaluated at 300 K is strongly
    # negative, so the clamp changes the sign of the result, not its magnitude.
    raw_cubic_300 = (14.43877 - 1.691 * 300.0 + 0.10687 * 300.0**2 - 0.00174 * 300.0**3) * 2.02
    assert raw_cubic_300 < 0.0
    assert cp_hot > 0.0
