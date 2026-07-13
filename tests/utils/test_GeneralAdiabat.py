"""Tests for src/janus/utils/GeneralAdiabat.py.

GeneralAdiabat builds a self-consistent multi-species moist pseudoadiabat from
the ground up. This file exercises the closed-form thermodynamic helpers and the
column-construction machinery that the integrator drives:

* ``p_sat`` and ``Tdew``: per-volatile saturation vapour pressure and its dew
  point inverse, including the dry-forced (H2) and unknown-species sentinels.
* ``get_beta``, ``get_T_crit``, ``L_heat``: the Clausius-Clapeyron phase-gradient,
  the per-volatile critical temperature, and the latent-heat dispatch across the
  sublimation, vaporization, and super-critical regimes.
* ``slopeRay``: the single-condensible (water/air) moist adiabat slope and its
  dry limit.
* ``dry_adiabat`` / ``dry_adiabat_pressure``: the Poisson dry adiabat and its
  pressure inverse.
* ``atm_z``: the hydrostatic height step and its non-monotone-height guard.
* ``condensation`` and ``general_adiabat``: surface mixing-ratio renormalization,
  wet/dry classification, and the assembled column profile.
* ``plot_adiabats``: the diagnostic figure (rendered on the Agg backend).

Invariant families exercised: positivity/boundedness of saturation pressures,
latent heats, heat-capacity-driven slopes and mixing ratios; monotonicity of
saturation pressure and dew point with temperature and of the dry adiabat with
pressure; conservation of the surface volatile inventory under renormalization;
and pinned analytical limits (the triple-point identity, the dew-point round trip,
the dry adiabat lapse R/cp). See docs/How-to/test.md.
"""

import math

import numpy as np
import pytest

import janus.utils.GeneralAdiabat as ga
import janus.utils.phys as phys
from janus.utils.atmosphere_column import atmos

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Volatiles whose p_sat / Tdew / L_heat match-cases are exercised here. H2 is
# handled separately because it is force-dried in p_sat.
_CONDENSABLES = ['CH4', 'CO', 'O2', 'O3', 'He', 'NH3']

# Hardcoded dew-point reference temperatures from the source (the arbitrary
# coexistence-curve anchor point for each volatile). Paired with the reference
# pressure p_sat(vol, Tref) they form the identity Tdew(vol, p_sat(vol, Tref)) = Tref.
_TDEW_REF = {
    'CH4': 148.15,
    'CO': 100.0,
    'N2': 98.15,
    'O2': 123.15,
    'O3': 161.85,
    'H2': 23.15,
    'NH3': 273.15,
    'CO2': 253.0,  # p > p_triple branch anchor
}


def _make_atm(vol_mixing, T_surf=300.0, P_surf=1.0e5, P_top=1.0e4, req_levels=15):
    """Construct a bare atmosphere with synthetic (empty) band edges, no data I/O."""
    return atmos(
        T_surf,
        P_surf,
        P_top,
        6.371e6,
        5.972e24,
        [],
        vol_mixing=vol_mixing,
        req_levels=req_levels,
    )


def _expected_L_heat(vol, T):
    """Reproduce the source latent-heat dispatch for the discrimination guard."""
    g = getattr(phys, vol)
    L_sub = g.L_sublimation
    if vol in ('H2', 'He'):
        # Source substitutes the vaporization value; these gases have no
        # tabulated sublimation latent heat.
        L_sub = g.L_vaporization
    L_vap = g.L_vaporization
    mw = g.MolecularWeight
    if T <= g.TriplePointT:
        return L_sub * mw * 1e-3
    if T <= g.CriticalPointT:
        return L_vap * mw * 1e-3
    return 0.0


# ---------------------------------------------------------------------------
# p_sat
# ---------------------------------------------------------------------------


@pytest.mark.physics_invariant
@pytest.mark.parametrize('vol', _CONDENSABLES)
def test_p_sat_positive_and_monotonic_with_temperature(vol):
    """Saturation vapour pressure is positive and rises with temperature.

    For each condensable the Clausius-Clapeyron saturation pressure must be
    strictly positive and finite, and it must increase with temperature. Two
    temperatures at 0.6 and 0.9 of the critical temperature keep the value away
    from float underflow while resolving the exponential slope, so a flipped
    temperature ratio (which would make p_sat fall with T) fails the monotonicity
    check. The colder point is the boundary case nearest underflow.
    """
    t_crit = ga.get_T_crit(vol)
    t_low = 0.6 * t_crit
    t_high = 0.9 * t_crit
    p_low = ga.p_sat(vol, t_low)
    p_high = ga.p_sat(vol, t_high)
    assert np.isfinite(p_low) and np.isfinite(p_high)
    assert p_low > 0.0
    # Monotone increasing: the warmer level has the higher saturation pressure.
    assert p_high > p_low


def test_p_sat_dry_forced_and_unknown_return_large_sentinel():
    """H2 and unknown species report an effectively infinite saturation pressure.

    H2 is force-dried in the source, so p_sat('H2', T) short-circuits to the 1e30
    sentinel that keeps hydrogen non-condensing at any temperature. An unknown
    species hits the match default and returns the same sentinel. Both are the
    contract paths that flag a volatile as never saturating; a real condensable
    (water) must instead return a finite laboratory-scale pressure, discriminating
    a regression that force-dried everything.
    """
    assert ga.p_sat('H2', 300.0) == pytest.approx(1.0e30, rel=1e-12)
    assert ga.p_sat('NotAGas', 300.0) == pytest.approx(1.0e30, rel=1e-12)
    # Contrast with a genuine condensable: finite and far below the sentinel.
    p_water = ga.p_sat('H2O', 300.0)
    assert p_water < 1.0e6
    assert p_water > 0.0


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
@pytest.mark.parametrize('vol', ['CO2', 'CH4', 'N2'])
def test_p_sat_triple_point_identity(vol):
    """At the triple point the saturation pressure equals the triple pressure.

    The Clausius-Clapeyron exponential collapses to unity at the reference
    temperature, so p_sat(vol, T_triple) must reproduce the tabulated triple-point
    pressure exactly (the analytical fixed point). This pins the pre-exponential
    scale e0; a unit-conversion slip would move the value off the tabulated Pa
    figure. A point 30 K above the triple temperature is strictly larger, so the
    identity is not a degenerate any-scale pass.
    """
    g = getattr(phys, vol)
    val = ga.p_sat(vol, g.TriplePointT)
    assert val == pytest.approx(g.TriplePointP, rel=1e-9)
    assert val > 0.0
    # Discrimination: 30 K warmer must exceed the triple-point pressure.
    assert ga.p_sat(vol, g.TriplePointT + 30.0) > val


# ---------------------------------------------------------------------------
# Tdew
# ---------------------------------------------------------------------------


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
@pytest.mark.parametrize('vol', list(_TDEW_REF.keys()) + ['He'])
def test_tdew_round_trip_recovers_reference_temperature(vol):
    """Tdew is the closed-form inverse of the saturation curve at the anchor.

    Each volatile's dew point is defined so that Tdew(vol, p_sat(vol, Tref)) = Tref,
    the analytical inverse identity of the constant-latent-heat Clausius-Clapeyron
    relation. Feeding the reference saturation pressure back through Tdew must
    recover the anchor temperature to machine precision. Helium uses the source's
    fixed 1e5 Pa reference. A pressure well above the anchor yields a strictly
    warmer dew point, so the identity is not degenerate.
    """
    if vol == 'He':
        tref = 4.22
        pref = 1.0e5
    else:
        tref = _TDEW_REF[vol]
        pref = ga.p_sat(vol, tref)
    val = ga.Tdew(vol, pref)
    assert val == pytest.approx(tref, rel=1e-9)
    assert val > 0.0
    # Monotonic discrimination: a higher pressure raises the dew point.
    assert ga.Tdew(vol, pref * 1.5) > val


@pytest.mark.physics_invariant
def test_tdew_monotonic_and_co2_sublimation_branch():
    """Dew point rises with pressure and the CO2 sublimation branch stays valid.

    The dew point increases with the imposed partial pressure (Clausius-Clapeyron
    run backwards), verified for CH4 and N2 near their reference pressures. Below
    their triple pressures both CO2 and water take their sublimation branches,
    which must still return finite, strictly positive temperatures; because the
    imposed pressure sits below the triple-point reference pressure the recovered
    dew point lies below the sublimation anchor (216.58 K for CO2, 273.15 K for
    water), discriminating the vaporization branch that would sit above it.
    """
    for vol in ('CH4', 'N2'):
        p_ref = ga.p_sat(vol, _TDEW_REF[vol])
        assert ga.Tdew(vol, 2.0 * p_ref) > ga.Tdew(vol, p_ref)

    # CO2 below its triple pressure (5.11e5 Pa) uses the sublimation reference.
    t_sub = ga.Tdew('CO2', 1.0e5)
    assert np.isfinite(t_sub)
    assert t_sub > 0.0
    assert t_sub < 216.58

    # Water below its triple pressure (611.657 Pa) uses the sublimation branch.
    t_sub_h2o = ga.Tdew('H2O', 100.0)
    assert np.isfinite(t_sub_h2o)
    assert t_sub_h2o > 0.0
    assert t_sub_h2o < 273.15


def test_tdew_unknown_species_returns_zero():
    """An unrecognized species yields a zero dew point sentinel.

    The match default returns 0.0 K, signalling "no dew point" for a species the
    routine does not know. This is the contract path; a known condensable must
    instead return a physical, strictly positive dew point, discriminating a
    regression that lost the species dispatch entirely.
    """
    assert ga.Tdew('NotAGas', 1.0e5) == pytest.approx(0.0, abs=1e-12)
    assert ga.Tdew('N2', 1.0e5) > 50.0


# ---------------------------------------------------------------------------
# get_beta
# ---------------------------------------------------------------------------


@pytest.mark.physics_invariant
def test_get_beta_water_lookup_differs_from_constant_formula():
    """The water lookup phase gradient differs from the constant-L formula.

    Inside the water triple-to-critical window with water_lookup enabled, get_beta
    reads the steam-table phase gradient d ln(psat)/d ln(T); disabled, it uses the
    constant-latent-heat form L/(R T). Both are strictly positive, and at 400 K the
    two must disagree, so a regression that ignored the lookup flag would be
    caught. Below the triple point the lookup window does not apply, so both modes
    must coincide there (the boundary limit of the lookup branch).
    """
    beta_lookup = ga.get_beta('H2O', 400.0, water_lookup=True)
    beta_const = ga.get_beta('H2O', 400.0, water_lookup=False)
    assert beta_lookup > 0.0
    assert beta_const > 0.0
    # Inside the window the two modes genuinely differ.
    assert abs(beta_lookup - beta_const) > 0.5
    # Below the triple point the lookup is inactive; both modes agree.
    beta_cold_lookup = ga.get_beta('H2O', 250.0, water_lookup=True)
    beta_cold_const = ga.get_beta('H2O', 250.0, water_lookup=False)
    assert beta_cold_lookup == pytest.approx(beta_cold_const, rel=1e-12)


# ---------------------------------------------------------------------------
# get_T_crit
# ---------------------------------------------------------------------------


@pytest.mark.physics_invariant
@pytest.mark.parametrize('vol', ['CH4', 'CO2', 'CO', 'N2', 'O2', 'O3', 'H2', 'He', 'NH3'])
def test_get_T_crit_matches_tabulated_critical_temperature(vol):
    """The critical temperature equals the tabulated per-volatile value.

    get_T_crit dispatches by species to the critical temperature stored in phys.
    Each returned value must equal that constant and be strictly positive (a
    physical critical temperature). An unknown species falls to the default 0.0 K
    sentinel, the contract path that the known-species pin discriminates against.
    """
    g = getattr(phys, vol)
    val = ga.get_T_crit(vol)
    assert val == pytest.approx(g.CriticalPointT, rel=1e-12)
    assert val > 0.0
    # Contract: an unrecognized species returns the zero sentinel.
    assert ga.get_T_crit('NotAGas') == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# L_heat
# ---------------------------------------------------------------------------


@pytest.mark.physics_invariant
@pytest.mark.parametrize('vol', ['CH4', 'CO2', 'CO', 'N2', 'O2', 'O3', 'H2', 'He', 'NH3'])
def test_L_heat_selects_phase_and_vanishes_above_critical(vol):
    """Latent heat follows the phase regime and drops to zero above critical.

    Below the triple point L_heat returns the sublimation latent heat, between the
    triple and critical points the vaporization value, and above the critical
    point exactly zero (no distinct condensed phase). Each regime is pinned to the
    molar conversion L * MolecularWeight * 1e-3 from phys, so a dropped 1e-3 or a
    wrong molecular weight fails. The super-critical zero is the boundary limit;
    the vaporization value is strictly positive for every condensable.
    """
    g = getattr(phys, vol)
    t_below = 0.9 * g.TriplePointT
    t_between = 0.5 * (g.TriplePointT + g.CriticalPointT)
    t_above = 1.1 * g.CriticalPointT

    assert ga.L_heat(vol, t_below) == pytest.approx(_expected_L_heat(vol, t_below), rel=1e-9)
    l_vap = ga.L_heat(vol, t_between)
    assert l_vap == pytest.approx(_expected_L_heat(vol, t_between), rel=1e-9)
    assert l_vap > 0.0
    # Super-critical: no condensed phase, latent heat vanishes.
    assert ga.L_heat(vol, t_above) == pytest.approx(0.0, abs=1e-30)


@pytest.mark.physics_invariant
def test_L_heat_water_lookup_differs_and_unknown_is_zero():
    """Water lookup latent heat differs from the constant value; unknown is zero.

    With water_lookup enabled in the vaporization regime, L_heat reads the
    temperature-dependent steam-table latent heat, which at 400 K differs from the
    constant phys value; both are strictly positive. An unknown species returns
    zero at every temperature (the default match with no tabulated data), the
    contract path.
    """
    l_lookup = ga.L_heat('H2O', 400.0, water_lookup=True)
    l_const = ga.L_heat('H2O', 400.0, water_lookup=False)
    assert l_lookup > 0.0
    assert l_const > 0.0
    assert abs(l_lookup - l_const) > 100.0
    # Unknown species: no latent heat in any regime.
    assert ga.L_heat('NotAGas', 300.0) == pytest.approx(0.0, abs=1e-30)


# ---------------------------------------------------------------------------
# slopeRay
# ---------------------------------------------------------------------------


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_slopeRay_dry_limit_and_moist_reduction():
    """The water/air slope reaches the dry adiabat lapse and flattens when moist.

    slopeRay returns d ln T / d ln p for a saturated water/air mixture. As the
    water loading vanishes (cold, high total pressure) it must approach the dry
    adiabat lapse R_air / cp_air, the analytical dry limit. When water is
    abundant (warm, low pressure) latent heating flattens the adiabat, so the
    moist slope is strictly smaller. Both slopes are strictly positive.
    """
    dry_limit = phys.air.R / phys.air.cp
    slope_dry = ga.slopeRay(math.log(1.0e5), math.log(200.0))
    slope_moist = ga.slopeRay(math.log(1.0e4), math.log(350.0))
    # Dry limit pins the analytical lapse to a few parts in 1e3.
    assert slope_dry == pytest.approx(dry_limit, rel=1e-2)
    assert slope_dry > 0.0
    assert slope_moist > 0.0
    # Latent heating flattens the adiabat: the moist slope is smaller.
    assert slope_moist < slope_dry


# ---------------------------------------------------------------------------
# dry_adiabat / dry_adiabat_pressure
# ---------------------------------------------------------------------------


@pytest.mark.physics_invariant
def test_dry_adiabat_poisson_profile():
    """The dry adiabat follows the Poisson relation and skips zero heat capacity.

    dry_adiabat returns T_surf (p / p_surf) ** (R / cp). At the surface pressure
    the ratio is unity, so temperature equals the surface value (the identity
    limit); aloft the temperature falls with pressure. A molar heat capacity of
    29 J/mol/K gives an exponent near 2/7 so the profile is well separated from
    the surface value, and a level with zero heat capacity is left at zero by the
    guard. A flipped exponent (cp/R rather than R/cp) would move the mid-level
    temperature by more than a factor of two, far past tolerance.
    """
    p_arr = np.array([1.0e5, 1.0e4, 1.0e3])
    cp = 29.0
    t_surf = 300.0
    t_dry = ga.dry_adiabat(t_surf, p_arr, cp)

    # Surface identity: at the maximum pressure the profile equals T_surf.
    assert t_dry[0] == pytest.approx(t_surf, rel=1e-12)
    # Poisson value aloft, pinned against the closed form.
    expected_mid = t_surf * (1.0e4 / 1.0e5) ** (phys.R_gas / cp)
    assert t_dry[1] == pytest.approx(expected_mid, rel=1e-9)
    # Monotone increasing with pressure.
    assert t_dry[2] < t_dry[1] < t_dry[0]
    # Wrong exponent (cp/R) would be far away, so the pin discriminates it.
    wrong_mid = t_surf * (1.0e4 / 1.0e5) ** (cp / phys.R_gas)
    assert abs(expected_mid - wrong_mid) > 1.0

    # Guard: a level with non-positive heat capacity is left at zero.
    t_guard = ga.dry_adiabat(t_surf, p_arr, [29.0, 0.0, 29.0])
    assert t_guard[1] == pytest.approx(0.0, abs=1e-30)


@pytest.mark.physics_invariant
def test_dry_adiabat_pressure_inverts_temperature_profile():
    """The pressure form inverts the Poisson relation and skips zero heat capacity.

    dry_adiabat_pressure returns P_surf (T / T_surf) ** (cp / R). At the maximum
    temperature the pressure equals P_surf (the identity limit); cooler levels sit
    at lower pressure. The mid-level value is pinned against the closed form, the
    profile decreases monotonically with temperature, and a level with zero heat
    capacity is left at zero by the guard. A flipped exponent (R/cp) would shift
    the mid-level pressure well past tolerance.
    """
    t_arr = np.array([300.0, 200.0, 150.0])
    cp = 29.0
    p_surf = 1.0e5
    p_dry = ga.dry_adiabat_pressure(p_surf, t_arr, cp)

    assert p_dry[0] == pytest.approx(p_surf, rel=1e-12)
    expected_mid = p_surf * (200.0 / 300.0) ** (cp / phys.R_gas)
    assert p_dry[1] == pytest.approx(expected_mid, rel=1e-9)
    assert p_dry[2] < p_dry[1] < p_dry[0]
    wrong_mid = p_surf * (200.0 / 300.0) ** (phys.R_gas / cp)
    assert abs(expected_mid - wrong_mid) > 1.0

    p_guard = ga.dry_adiabat_pressure(p_surf, t_arr, [29.0, 0.0, 29.0])
    assert p_guard[1] == pytest.approx(0.0, abs=1e-30)


# ---------------------------------------------------------------------------
# atm_z
# ---------------------------------------------------------------------------


@pytest.mark.physics_invariant
def test_atm_z_warns_when_height_step_descends(caplog):
    """A pressure inversion forces a descending height step and logs a warning.

    atm_z integrates hydrostatic height upward: dz = -dp / rho. A physical column
    has pressure falling with height, so dz is positive. Forcing the next level to
    a higher pressure than the current one drives dp positive and dz negative; the
    routine must both compute the negative step (height decreases across the
    inversion) and emit a warning naming dz. The height decrement is the numeric
    contract that the log line only annotates.
    """
    atm = _make_atm({'H2O': 0.1, 'CO2': 0.2, 'N2': 0.7})
    atm.p[0] = 1.0e5
    atm.p[1] = 2.0e5  # inverted: pressure rises with height
    atm.pl[0] = 1.0e5
    atm.pl[1] = 2.0e5
    atm.mu[0] = 0.029
    atm.tmp[0] = 300.0
    atm.xd[0] = 1.0  # non-zero gas fraction so density stays finite
    atm.xc[0] = 0.0
    atm.z[0] = 0.0

    with caplog.at_level('WARNING', logger='fwl.janus.utils.GeneralAdiabat'):
        ga.atm_z(atm, 0)

    # The height step is negative across the inversion.
    assert atm.z[1] < atm.z[0]
    assert atm.z[1] < 0.0
    # The routine surfaces the non-physical step in the log.
    assert any('dz' in rec.message for rec in caplog.records)


# ---------------------------------------------------------------------------
# condensation
# ---------------------------------------------------------------------------


@pytest.mark.physics_invariant
def test_condensation_renormalizes_inventory_and_classifies_species():
    """Surface condensation renormalizes the mixing ratios and sorts wet vs dry.

    When the surface mixing ratios do not sum to one, condensation renormalizes
    them so the total volatile inventory is conserved (sum returns to unity). It
    then classifies each species: water at 300 K is supersaturated and enters the
    wet list, while N2 and CO2 stay far below saturation and enter the dry list.
    The renormalized partial pressures are strictly positive. A regression that
    dropped the renormalization would leave the inflated sum in place.
    """
    atm = _make_atm({'H2O': 0.1, 'CO2': 0.2, 'N2': 0.7})
    atm.vol_list['N2'] = atm.vol_list['N2'] + 0.5  # break the unit-sum invariant
    atm.xd[0] = 1e-10  # avoid a divide-by-zero when called outside the driver
    assert sum(atm.vol_list.values()) > 1.0

    atm, wet_list, dry_list = ga.condensation(atm, 0, [], [], prs_reset=False)

    # Conservation: the renormalized mixing ratios sum back to unity.
    assert sum(atm.vol_list.values()) == pytest.approx(1.0, rel=1e-9)
    # Water is supersaturated at the surface; N2 and CO2 are dry.
    assert 'H2O' in wet_list
    assert 'N2' in dry_list
    assert 'CO2' in dry_list
    # Partial pressures stay strictly positive after renormalization.
    assert atm.p_vol['H2O'][0] > 0.0
    assert atm.p_vol['N2'][0] > 0.0


# ---------------------------------------------------------------------------
# general_adiabat
# ---------------------------------------------------------------------------


@pytest.mark.physics_invariant
def test_general_adiabat_profile_monotone_and_bounded():
    """The assembled adiabat cools upward with a bounded, conserved gas phase.

    Integrating the general moist adiabat for a water/CO2/N2 surface must return a
    column whose pressure increases with index (top to surface) and whose
    temperature decreases upward, all strictly positive. The dry and wet gas molar
    concentrations stay within [0, 1], and with full rainout their sum equals one
    at every level (the gas-phase closure). A regression that inverted the profile
    would warm the upper column instead of cooling it.
    """
    atm = _make_atm({'H2O': 0.3, 'CO2': 0.2, 'N2': 0.5}, req_levels=15)
    atm = ga.general_adiabat(atm)

    assert np.all(atm.tmp > 0.0)
    assert np.all(atm.p > 0.0)
    # Pressure ascends with index; temperature descends toward the top.
    assert np.all(np.diff(atm.p) > 0.0)
    assert atm.tmp[0] < atm.tmp[-1]
    assert np.all(np.diff(atm.tmp) >= -1e-9)
    # Molar concentrations are bounded and the gas phase closes to unity.
    assert np.all((atm.xd >= 0.0) & (atm.xd <= 1.0 + 1e-9))
    assert np.all((atm.xv >= 0.0) & (atm.xv <= 1.0 + 1e-9))
    np.testing.assert_allclose(atm.xd + atm.xv, 1.0, rtol=1e-6)


# ---------------------------------------------------------------------------
# plot_adiabats
# ---------------------------------------------------------------------------


@pytest.mark.physics_invariant
def test_plot_adiabats_renders_labelled_axes(tmp_path, monkeypatch):
    """The diagnostic figure carries labelled axes and one line set per species.

    plot_adiabats draws the saturation curves, partial pressures, dry and general
    adiabats, and the phase molar concentrations for every volatile above the 1e-10
    presence floor. The left panel must be labelled in temperature and pressure and
    hold two lines per plotted species plus the three reference curves; the right
    panel holds two lines per species plus the gas-phase curve. A trace species
    (1e-15) stays below the floor and adds no lines, exercising the presence guard.
    The saved PDF is non-empty.
    """
    atm = _make_atm({'H2O': 0.3, 'CO2': 0.2, 'N2': 0.5, 'CH4': 1e-15}, req_levels=15)
    atm = ga.general_adiabat(atm)
    n_plotted = sum(1 for v in atm.vol_list if atm.vol_list[v] > 1e-10)

    captured = {}
    real_subplots = ga.plt.subplots

    def _spy(*args, **kwargs):
        fig, axs = real_subplots(*args, **kwargs)
        captured['fig'] = fig
        captured['axs'] = axs
        return fig, axs

    monkeypatch.setattr(ga.plt, 'subplots', _spy)

    out_file = tmp_path / 'general_adiabat.pdf'
    ga.plot_adiabats(atm, filename=str(out_file))

    ax1, ax2 = captured['axs']
    assert 'Temperature' in ax1.get_xlabel()
    assert 'Pressure' in ax1.get_ylabel()
    # Two lines per plotted species plus three references (partial-pressure sum,
    # dry adiabat, general adiabat); the trace species adds none.
    assert len(ax1.get_lines()) == 2 * n_plotted + 3
    assert len(ax2.get_lines()) == 2 * n_plotted + 1
    assert n_plotted == 3
    assert out_file.exists()
    assert out_file.stat().st_size > 0

    ga.plt.close('all')
