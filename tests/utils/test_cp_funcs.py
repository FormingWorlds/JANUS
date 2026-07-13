"""Tests for src/janus/utils/cp_funcs.py.

cp_funcs.py provides the per-species molar heat capacity ``cpv`` with two
branches (a constant value and a temperature-dependent NIST Shomate
polynomial) and the condensate heat capacity ``cp_cond``. This file exercises:

* The Shomate branch for water vapour pinned against the NIST-JANAF value
  at 500 K.
* The Shomate branch for every tabulated volatile (H2O, CO2, H2, N2, CH4,
  CO, O2, He, NH3) across each of its temperature sub-ranges, pinned against
  the NIST Chemistry WebBook molar heat capacities.
* The below-fit-window temperature clamp of the Shomate branch.
* The mode dispatch of ``cpv``: the constant branch for every species, the
  divergence between the constant and T-dependent branches at high T, and
  the None return for an unrecognised ``cp_mode``.
* The unknown-species guard return (cp=0) in both branches.
* Positivity of a known species' heat capacity and rough agreement between
  the two branches.
* The condensate ``cp_cond`` for every species: the constant-mode
  representative values, the saturated-liquid interpolation with its
  triple-point / critical-point clamps, the fixed-temperature liquid-hydrogen
  branch, the liquid-oxygen literature value, and the unknown-species guard.

Invariant families exercised: positivity/boundedness, monotonicity, and
pinned-value-with-discrimination-guard. See docs/How-to/test.md.
"""

import math

import pytest

import janus.utils.cp_funcs as cp_funcs

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Molar gas constant (CODATA 2018), J mol-1 K-1. The monatomic ideal-gas heat
# capacity is 5/2 R, the analytical limit the helium Shomate fit reproduces.
R_GAS = 8.314462618


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


# NIST Chemistry WebBook Shomate molar heat capacities (J mol-1 K-1) of the gas
# phase near the low end of each species' fit window. The source reproduces
# these published values from the Shomate coefficients; each tuple is
# (species, temperature_K, nist_cp).
_NIST_CPV_LOW_T = [
    ('H2O', 500.0, 35.22),
    ('CO2', 298.0, 37.13),
    ('H2', 298.0, 28.84),
    ('N2', 300.0, 29.12),
    ('CH4', 298.0, 35.64),
    ('CO', 298.0, 29.15),
    ('O2', 298.0, 29.38),
    ('NH3', 298.0, 35.65),
]


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    'vol,tmp,nist_cp',
    _NIST_CPV_LOW_T,
    ids=[
        'water-vapour-500K',
        'carbon-dioxide-298K',
        'hydrogen-298K',
        'nitrogen-300K',
        'methane-298K',
        'carbon-monoxide-298K',
        'oxygen-298K',
        'ammonia-298K',
    ],
)
def test_shomate_cpv_reproduces_nist_low_temperature_reference(vol, tmp, nist_cp):
    """Each volatile's Shomate cp matches the NIST WebBook value at low T.

    The T-dependent branch selects the lowest temperature sub-range and
    evaluates the NIST Shomate polynomial. Pinning against the published NIST
    Chemistry WebBook molar heat capacity for every tabulated species confirms
    the per-species coefficient block and the lowest-index selection. The value
    at 3000 K, which routes through the highest sub-range, must exceed the
    low-temperature value because internal modes activate with temperature;
    that monotone rise discriminates a wrong-index selection that would flatten
    the curve.
    """
    cp = cp_funcs.cpv(vol, tmp, 'T-dependent')
    # Primary pin: matches the NIST tabulated molar cp to about 1%.
    assert cp == pytest.approx(nist_cp, rel=1.2e-2)
    # Sign / positivity guard: a heat capacity is strictly positive.
    assert cp > 0.0
    # Scale / unit guard: molar cp of a small-molecule gas sits in 15 to 95
    # J/mol/K. A per-gram slip (dividing by molar mass) would land near 1 to 3;
    # a factor-1000 slip would land far above this band.
    assert 15.0 < cp < 95.0
    # Exponent / index guard: the highest sub-range value at 3000 K exceeds the
    # low-temperature value for every one of these polyatomic / diatomic gases.
    cp_hot = cp_funcs.cpv(vol, 3000.0, 'T-dependent')
    assert cp_hot > cp


# NIST WebBook molar heat capacities in the MIDDLE Shomate sub-range of the
# three-range species; (species, low_T, mid_T, mid_nist_cp, high_T).
_NIST_CPV_MID_RANGE = [
    ('H2', 298.0, 2000.0, 34.28, 3000.0),
    ('N2', 300.0, 1000.0, 32.69, 3000.0),
    ('O2', 298.0, 1000.0, 34.86, 3000.0),
]


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    'vol,low_t,mid_t,mid_cp,high_t',
    _NIST_CPV_MID_RANGE,
    ids=['hydrogen-mid-range', 'nitrogen-mid-range', 'oxygen-mid-range'],
)
def test_shomate_cpv_selects_middle_temperature_subrange(vol, low_t, mid_t, mid_cp, high_t):
    """The Shomate branch selects the middle sub-range for H2, N2, O2.

    Hydrogen, nitrogen and oxygen each have three Shomate coefficient sets. A
    temperature inside the middle window must select the middle set and
    reproduce the NIST WebBook cp there. The middle value must lie strictly
    between the low-range and high-range values (heat capacity rises
    monotonically with temperature for these gases), which discriminates a
    regression that selected the lowest or highest set instead of the middle
    one.
    """
    cp_mid = cp_funcs.cpv(vol, mid_t, 'T-dependent')
    assert cp_mid == pytest.approx(mid_cp, rel=1.2e-2)
    assert cp_mid > 0.0
    # Monotone ordering across the three sub-ranges is the index-selection guard.
    cp_low = cp_funcs.cpv(vol, low_t, 'T-dependent')
    cp_high = cp_funcs.cpv(vol, high_t, 'T-dependent')
    assert cp_low < cp_mid < cp_high


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_shomate_cpv_helium_matches_monatomic_ideal_limit():
    """Helium's Shomate cp equals the monatomic ideal-gas limit 5/2 R.

    Helium has a single Shomate set whose higher-order coefficients are
    negligible, so the T-dependent branch returns the monatomic ideal-gas heat
    capacity 5/2 R = 20.79 J/mol/K at every temperature. The value is
    essentially temperature-independent, which is the analytical limit and the
    edge case: a temperature below the 298 K fit floor clamps to the floor and
    still returns the same value. The near-perfect flatness discriminates the
    monatomic case from a diatomic gas (cp near 29 J/mol/K).
    """
    cp_warm = cp_funcs.cpv('He', 300.0, 'T-dependent')
    cp_hot = cp_funcs.cpv('He', 6000.0, 'T-dependent')
    cp_clamped = cp_funcs.cpv('He', 100.0, 'T-dependent')  # below the 298 K floor
    # Analytical limit: 5/2 R.
    assert cp_warm == pytest.approx(2.5 * R_GAS, rel=1e-3)
    # Temperature independence: monatomic cp has no rotational / vibrational
    # temperature dependence, so the 6000 K value matches the 300 K value.
    assert cp_hot == pytest.approx(cp_warm, rel=1e-6)
    # Clamp edge: below the fit floor the value is unchanged, still positive.
    assert cp_clamped == pytest.approx(cp_warm, rel=1e-6)
    assert cp_warm > 0.0
    # Discrimination from a diatomic gas: helium's cp is far below 29 J/mol/K.
    assert cp_warm < 24.0


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    'vol,cold_t,floor_t',
    [('H2O', 400.0, 500.0), ('CO2', 200.0, 298.0), ('N2', 50.0, 100.0)],
    ids=['water-below-500K-floor', 'co2-below-298K-floor', 'n2-below-100K-floor'],
)
def test_shomate_cpv_clamps_below_fit_window_floor(vol, cold_t, floor_t):
    """Below its fit window the Shomate branch clamps to the floor temperature.

    Each species restricts its Shomate fit to a validity window and clamps the
    evaluation temperature to the window floor. A temperature below the floor
    (the limit-input / error-contract path) must return exactly the value at
    the floor, and that value must be a positive, physically sized molar heat
    capacity. Equality with the floor value discriminates a missing clamp that
    would extrapolate the polynomial to an unphysical value below the window.
    """
    cp_cold = cp_funcs.cpv(vol, cold_t, 'T-dependent')
    cp_floor = cp_funcs.cpv(vol, floor_t, 'T-dependent')
    # Clamp identity: the sub-floor temperature evaluates at the floor.
    assert cp_cold == pytest.approx(cp_floor, rel=1e-12)
    assert cp_cold > 0.0
    assert 15.0 < cp_cold < 95.0


@pytest.mark.physics_invariant
@pytest.mark.parametrize('vol', ['H2O', 'CH4', 'CO2'])
def test_shomate_and_constant_branches_diverge_at_high_temperature(vol):
    """The constant and T-dependent branches separate at high temperature.

    At room temperature the constant heat capacity and the Shomate polynomial
    nearly agree, but by 2000 K the Shomate curve has risen well above the
    fixed constant value. Pinning the divergence enforces that ``cp_mode``
    actually dispatches between the two branches: a regression that silently
    routed the T-dependent request to the constant branch (or the reverse)
    would collapse the gap and fail this test.
    """
    cp_const = cp_funcs.cpv(vol, 2000.0, 'constant')
    cp_tdep = cp_funcs.cpv(vol, 2000.0, 'T-dependent')
    assert cp_const > 0.0
    assert cp_tdep > 0.0
    # The Shomate value rises with T; the constant value does not. The gap at
    # 2000 K is several J/mol/K, far outside any rounding tolerance.
    assert cp_tdep > cp_const + 3.0


# Constant-branch molar heat capacities are derived from phys.py per-species
# specific heats; every entry is positive and in the small-molecule band.
_CONSTANT_MODE_SPECIES = ['H2O', 'CH4', 'CO2', 'CO', 'N2', 'O2', 'H2', 'He', 'NH3']


@pytest.mark.physics_invariant
@pytest.mark.parametrize('vol', _CONSTANT_MODE_SPECIES)
def test_cpv_constant_branch_positive_and_physically_ordered(vol):
    """The constant branch returns a positive, physically sized cp per species.

    The constant branch multiplies each species' specific heat from phys.py by
    its molecular weight. Every species must return a strictly positive molar
    heat capacity in the small-molecule band. Helium, being monatomic, must
    return the smallest value of the set (fewest degrees of freedom); a unit or
    lookup error for any species would break that ordering or the band.
    """
    cp = cp_funcs.cpv(vol, 298.0, 'constant')
    assert cp > 0.0
    assert 15.0 < cp < 95.0
    # Physical ordering: the monatomic species has the lowest molar cp.
    cp_helium = cp_funcs.cpv('He', 298.0, 'constant')
    assert cp_helium <= cp
    # Helium sits near the monatomic 5/2 R limit, well below the diatomics.
    assert cp_helium == pytest.approx(2.5 * R_GAS, rel=5e-3)


@pytest.mark.physics_invariant
def test_cpv_unknown_mode_returns_none_while_valid_mode_returns_value():
    """An unrecognised cp_mode returns None; a valid mode returns a cp.

    ``cpv`` handles only the ``"T-dependent"`` and ``"constant"`` modes. Any
    other mode string matches neither branch and the function falls through to
    an implicit None return (the mode contract / error path). A recognised mode
    on the same species still returns a positive molar heat capacity, so the
    None is specific to the bad mode and not a broken species lookup.
    """
    assert cp_funcs.cpv('H2O', 500.0, 'nonsense-mode') is None
    cp_valid = cp_funcs.cpv('H2O', 500.0, 'T-dependent')
    assert cp_valid == pytest.approx(35.22, rel=1.2e-2)
    assert cp_valid > 0.0


# Constant-mode condensate heat capacities (J mol-1 K-1) hard-coded in the
# source as representative saturated-liquid values; (species, constant_cp).
_CONDENSATE_CONSTANT = [
    ('H2O', 76.0),
    ('CO2', 100.0),
    ('N2', 58.0),
    ('CH4', 55.0),
    ('CO', 60.0),
    ('NH3', 81.0),
]


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    'vol,const_cp',
    _CONDENSATE_CONSTANT,
    ids=['water', 'co2', 'nitrogen', 'methane', 'carbon-monoxide', 'ammonia'],
)
def test_condensate_constant_mode_returns_representative_value(vol, const_cp):
    """Constant-mode ``cp_cond`` returns the fixed representative liquid value.

    In constant mode each condensate returns a single representative
    saturated-liquid molar heat capacity. Pinning the value enforces the mode
    dispatch; the constant value must also differ from the temperature-resolved
    value at the same temperature, which discriminates a regression that
    ignored ``cp_mode`` and always interpolated. Every value is a positive,
    physically sized liquid molar heat capacity.
    """
    cp_const = float(cp_funcs.cp_cond(vol, 250.0, 'constant'))
    assert cp_const == pytest.approx(const_cp, rel=1e-9)
    assert cp_const > 0.0
    # Mode guard: the interpolated value at the same T differs from the constant.
    cp_interp = float(cp_funcs.cp_cond(vol, 250.0, 'T-dependent'))
    assert cp_interp != pytest.approx(cp_const, rel=1e-6)


@pytest.mark.physics_invariant
def test_condensate_hydrogen_constant_mode_uses_fixed_evaluation_temperature():
    """Constant-mode liquid-hydrogen ``cp_cond`` evaluates its cubic at 10 K.

    The liquid-hydrogen branch has no interpolation table; in constant mode it
    evaluates the cubic specific-heat fit at a fixed 10 K rather than the
    requested temperature. The returned molar heat capacity must be positive
    and must differ from the temperature-resolved value at the triple-to-critical
    window, which discriminates a regression that ignored the constant mode and
    evaluated the cubic at the requested temperature instead.
    """
    cp_const = cp_funcs.cp_cond('H2', 30.0, 'constant')
    cp_tdep = cp_funcs.cp_cond('H2', 30.0, 'T-dependent')
    assert cp_const > 0.0
    # Fixed 10 K evaluation: 14.43877 - 1.691*10 + 0.10687*100 - 0.00174*1000,
    # times the 2.02 g/mol molar mass, is about 13.1 J/mol/K.
    assert cp_const == pytest.approx(13.08, rel=1e-3)
    # Mode guard: at 30 K the T-dependent cubic gives a distinctly larger value.
    assert cp_tdep > cp_const + 5.0


# Saturated-liquid condensate endpoints from the NIST WebBook fluid tables;
# (species, below_triple_T, triple_cp, triple_T, above_critical_T, critical_cp).
_CONDENSATE_TABLE = [
    ('H2O', 270.0, 75.97, 274.0, 700.0, 3685.6),
    ('CO2', 210.0, 86.0, 217.0, 310.0, 694.8),
    ('N2', 50.0, 56.07, 64.0, 200.0, 271.2),
    ('CH4', 80.0, 54.05, 91.0, 250.0, 1508.2),
    ('CO', 60.0, 60.33, 69.0, 200.0, 672.65),
    ('NH3', 180.0, 71.61, 196.0, 500.0, 5907.5),
]


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    'vol,cold_t,triple_cp,triple_t,hot_t,crit_cp',
    _CONDENSATE_TABLE,
    ids=['water', 'co2', 'nitrogen', 'methane', 'carbon-monoxide', 'ammonia'],
)
def test_condensate_temperature_dependent_clamps_and_diverges(
    vol, cold_t, triple_cp, triple_t, hot_t, crit_cp
):
    """T-dependent ``cp_cond`` clamps to the table ends and diverges near it.

    The saturated-liquid heat capacity is interpolated over a NIST fluid table.
    Below the triple point the evaluation temperature clamps to the first table
    entry; above the critical point it clamps to the last. A temperature below
    the triple point (the edge case / clamp contract) must return the
    triple-point value, and a temperature above the critical point must return
    the near-critical value, which is finite and far larger than the
    triple-point value as the liquid heat capacity diverges toward the critical
    point. The clamp identities and the divergence together discriminate a
    missing clamp (which would extrapolate off the table) from the correct
    behaviour.
    """
    cp_cold = float(cp_funcs.cp_cond(vol, cold_t, 'T-dependent'))
    cp_triple = float(cp_funcs.cp_cond(vol, triple_t, 'T-dependent'))
    cp_hot = float(cp_funcs.cp_cond(vol, hot_t, 'T-dependent'))
    # Floor clamp: below the triple point returns the first table entry.
    assert cp_cold == pytest.approx(triple_cp, rel=1e-3)
    assert cp_cold == pytest.approx(cp_triple, rel=1e-12)
    # Ceiling clamp: above the critical point returns the last table entry.
    assert cp_hot == pytest.approx(crit_cp, rel=1e-3)
    assert math.isfinite(cp_hot)
    # Physical divergence: the near-critical heat capacity exceeds twice the
    # triple-point value, and both are positive.
    assert cp_hot > 2.0 * cp_cold
    assert cp_cold > 0.0
    # Mid-window interpolation stays bounded between the extremes and positive.
    cp_mid = float(cp_funcs.cp_cond(vol, 0.5 * (triple_t + hot_t), 'T-dependent'))
    assert cp_cold <= cp_mid <= cp_hot


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_condensate_oxygen_returns_liquid_reference_value():
    """Liquid-oxygen ``cp_cond`` returns the Giauque and Johnston 1929 value.

    The liquid-oxygen branch returns a single literature molar heat capacity,
    4.184 * 10 = 41.84 J/mol/K (Giauque and Johnston 1929), independent of both
    temperature and mode. Requesting different temperatures and both modes must
    return the same positive value; a mode-dependent or temperature-dependent
    result would signal a regression that added spurious dispatch to this
    branch.
    """
    cp_cold = cp_funcs.cp_cond('O2', 60.0, 'T-dependent')
    cp_warm = cp_funcs.cp_cond('O2', 150.0, 'constant')
    assert cp_cold == pytest.approx(41.84, rel=1e-6)
    assert cp_warm == pytest.approx(41.84, rel=1e-6)
    assert cp_cold > 0.0


@pytest.mark.physics_invariant
def test_condensate_unknown_species_returns_zero():
    """Unknown condensate species fall through to the cp=0 guard.

    ``cp_cond`` returns 0.0 for any species without a table entry, in both
    modes (the guard / error path). A known species returns a positive value at
    the same temperature, so the zero is specific to the missing entry and not
    a collapsed interpolation for everything.
    """
    assert cp_funcs.cp_cond('Xx_not_a_liquid', 250.0, 'T-dependent') == pytest.approx(
        0.0, abs=1e-12
    )
    assert cp_funcs.cp_cond('Xx_not_a_liquid', 250.0, 'constant') == pytest.approx(
        0.0, abs=1e-12
    )
    cp_known = float(cp_funcs.cp_cond('H2O', 300.0, 'constant'))
    assert cp_known > 0.0
