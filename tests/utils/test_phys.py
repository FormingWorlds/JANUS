"""Tests for src/janus/utils/phys.py.

phys.py holds the thermodynamic constants and saturation-vapour-pressure
relations the JANUS adiabat and height integrators depend on. This file
exercises:

* The simplified Clausius-Clapeyron law ``satvps`` and its exact
  reference-pressure identity ``satvps(T0, T0, e0, ...) == e0``.
* The Smithsonian liquid-water saturation curve ``satvpw`` against the
  ~1 atm steam-point anchor, and the ice curve ``satvpi``.
* Positivity of saturation pressures and monotonic increase with
  temperature (the Clausius-Clapeyron sign).
* Planck spectral radiance ``B``: the closed-form value at a pinned
  (nu, T), positivity, dB/dT > 0, the overflow-clamped Wien tail, and the
  zero-return guard at non-positive frequency.
* Self-consistency of the CODATA constant set through the Stefan-Boltzmann
  identity and the Rstar = 1000 k N_A gas-constant chain.
* Agreement between the ``molar_mass`` dict and the per-gas ``MolecularWeight``
  objects, and that a numerically constructed ``satvps_function`` never
  dispatches to the water lookup tables.

Invariant families exercised: positivity/boundedness, monotonicity, and
pinned-value-with-discrimination-guard. See docs/How-to/test.md.
"""

import math

import numpy as np
import pytest

import janus.utils.phys as phys
import janus.utils.water_tables as wt

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_satvps_returns_reference_pressure_at_reference_temperature():
    """Clausius-Clapeyron law recovers e0 at T=T0 and its exponent at T!=T0.

    ``satvps`` is the simplified Clausius-Clapeyron relation
    e(T)=e0*exp(-(L/Rv)(1/T-1/T0)). At the reference temperature the exponent
    vanishes, so e(T0) must equal e0 to machine precision; this degenerate
    anchor is insensitive to the L/Rv coefficient. Thirty kelvin below T0 the
    exponent is far from zero, so the value there pins the coefficient itself: a
    doubled or halved L/Rv (a factor slip in the latent heat or the specific gas
    constant) moves the result outside tolerance. Cooling must lower the
    saturation pressure and warming must raise it (dP/dT>0), discriminating a
    dropped minus sign in the exponent.
    """
    # Water-like reference point: e0=611 Pa at the triple point T0=273.16 K,
    # molecular weight 18 (g/mol), latent heat 2.5e6 J/kg.
    T0, e0, mol_weight, latent = 273.16, 611.0, 18.0, 2.5e6
    at_ref = phys.satvps(T0, T0, e0, mol_weight, latent)
    assert at_ref == pytest.approx(e0, rel=1e-12)  # exponent is exactly zero here

    colder = phys.satvps(T0 - 30.0, T0, e0, mol_weight, latent)
    warmer = phys.satvps(T0 + 30.0, T0, e0, mol_weight, latent)
    # Off-reference pin: 30 K below T0 the exponent is -2.44, giving 53.0 Pa for
    # the water-like L/Rv coefficient. Unlike the T0 identity this value depends
    # on the coefficient, so it pins L/Rv rather than collapsing to e0.
    assert colder == pytest.approx(53.02, rel=1e-2)
    # Coefficient guard: a factor slip in L or Rv leaves this band. Doubling the
    # coefficient collapses the value to ~4.6 Pa (below), halving it lifts the
    # value to ~180 Pa (above); the correct 53 Pa is the only value inside.
    assert 10.0 < colder < 100.0
    # Sign-of-exponent guard: a dropped minus sign would invert these.
    assert colder < e0 < warmer
    # Positivity: saturation pressure is strictly positive.
    assert colder > 0.0


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_satvpw_matches_one_atm_at_steam_point():
    """Smithsonian liquid-water curve gives ~1 atm at the 373.16 K steam point.

    ``satvpw`` is the Smithsonian saturation-vapour-pressure formula over
    liquid water. At the steam-point temperature 373.16 K it returns
    1.0151e5 Pa, 0.19% above one standard atmosphere (101325 Pa), consistent
    with the Smithsonian fit; this ~1 atm scale is the physical anchor. The
    curve must rise monotonically with temperature (Clausius-Clapeyron).
    """
    p_boil = phys.satvpw(373.16)
    # Smithsonian fit at the steam point; pins the regression value.
    assert p_boil == pytest.approx(101514.5, rel=1e-4)
    # Reference-scale guard: within ~0.5% of one standard atmosphere, ruling
    # out a Pa/hPa unit slip (would give ~1013 or ~1.01e7) or a 10x error.
    assert 1.00e5 < p_boil < 1.02e5
    # Clausius-Clapeyron monotonicity: warmer water has higher saturation
    # pressure, colder (still above the triple point) has lower.
    assert phys.satvpw(363.16) < p_boil < phys.satvpw(383.16)
    assert phys.satvpw(363.16) > 0.0


@pytest.mark.physics_invariant
def test_ice_saturation_below_liquid_and_both_positive():
    """Saturation over ice is positive and below that over supercooled water.

    Below the triple point the saturation vapour pressure over ice
    (``satvpi``) is lower than over supercooled liquid water (``satvpw``);
    this difference drives the Wegener-Bergeron-Findeisen process. Both must
    remain strictly positive.
    """
    temp = 253.16  # -20 C, the lower validity edge of the liquid formula
    p_ice = phys.satvpi(temp)
    p_liq = phys.satvpw(temp)
    assert p_ice > 0.0
    assert p_liq > 0.0
    assert p_ice < p_liq  # ice holds less vapour than supercooled water
    # The gap is several percent, resolvable well above rounding; a swapped
    # ice/water formula would flip the inequality.
    assert (p_liq - p_ice) / p_liq > 0.01


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_planck_radiance_positive_and_increases_with_temperature():
    """Planck radiance B(nu,T) matches its closed form and rises with T.

    ``B`` is the Planck function of frequency,
    B = (2 h nu**3 / c**2) / (exp(h nu / k T) - 1). At nu=5e13 Hz (~6 micron)
    and T=300 K the closed form evaluates to 6.193e-13 W m^-2 sr^-1 Hz^-1;
    pinning this value fixes the nu**3 prefactor and the h*nu/kT exponent
    together. For any physical (nu>0, T>0) the radiance is strictly positive,
    and at fixed frequency a hotter body radiates more (dB/dT>0). The overflow
    guard clamps h*nu/(kT) at 500 so the high-frequency Wien tail stays finite
    rather than overflowing; that clamp is the error-contract path.
    """
    nu = 5.0e13  # ~6 micron, thermal-IR
    b_cool = phys.B(nu, 300.0)
    b_warm = phys.B(nu, 1500.0)
    # Closed-form Planck value at (5e13 Hz, 300 K); pins the nu**3 prefactor and
    # the h*nu/kT exponent. A dropped factor of 2 or a nu**2 / nu**4 prefactor
    # slip leaves the value outside the surrounding band.
    assert b_cool == pytest.approx(6.193e-13, rel=1e-3)
    assert 5.0e-13 < b_cool < 7.0e-13
    assert b_cool > 0.0
    assert b_warm > b_cool  # dB/dT > 0 at fixed frequency
    # Wien-tail ratio: warming 300 -> 1500 K lifts the radiance ~750x at this
    # frequency, far from the ~5x a wrong linear-in-T law would give.
    assert b_warm / b_cool > 100.0
    # Edge case: an extreme frequency trips the u=500 overflow clamp and must
    # still return a finite, positive, vanishingly small radiance.
    b_extreme = phys.B(1.0e16, 300.0)
    assert math.isfinite(b_extreme)
    assert 0.0 < b_extreme < b_cool


@pytest.mark.physics_invariant
def test_planck_radiance_vanishes_at_nonpositive_frequency():
    """Planck radiance returns exactly zero at zero and negative frequency.

    B(nu, T) carries a nu**3 prefactor, so as nu -> 0 the spectral radiance
    tends to zero. At nu = 0 the closed form is a 0/0 division; the routine's
    guard returns 0.0 rather than raising ZeroDivisionError, and extends the
    same zero to unphysical negative frequencies. This guard is the error
    contract; the vanishing small-frequency limit is the physical edge case.
    """
    assert phys.B(0.0, 300.0) == pytest.approx(0.0, abs=1e-30)
    assert phys.B(-1.0e13, 300.0) == pytest.approx(0.0, abs=1e-30)
    # Approaching zero frequency from above, the radiance falls toward zero: in
    # the small-nu limit B ~ nu**2, so a decade drop in nu cuts B by ~100.
    b_small = phys.B(1.0e9, 300.0)
    b_smaller = phys.B(1.0e8, 300.0)
    assert b_small > 0.0
    assert b_smaller > 0.0
    assert b_smaller < b_small
    # Discrimination: the small-nu radiance is orders above the exact zero the
    # guard returns, so the guard is a genuine special case, not rounding.
    assert phys.B(0.0, 300.0) < b_smaller


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_fundamental_constants_are_codata_self_consistent():
    """The CODATA constant set closes the Stefan-Boltzmann identity.

    The Planck constant h, speed of light c, Boltzmann constant k and
    Stefan-Boltzmann constant sigma are not independent: sigma is fixed by the
    other three through sigma = 2 pi^5 k^4 / (15 h^3 c^2). With the SI-defined
    CODATA values this identity closes to the precision of the tabulated sigma.
    The universal gas constant is likewise tied to k and Avogadro's number by
    Rstar = 1000 k N_A, and the molar gas constant by R_gas = Rstar/1000. A
    mixed-vintage constant set breaks the closure at the fifth significant
    figure, which the discrimination guard rejects.
    """
    sigma_from_hck = 2.0 * math.pi**5 * phys.k**4 / (15.0 * phys.h**3 * phys.c**2)
    # The tabulated sigma reproduces the radiation-constant identity.
    assert phys.sigma == pytest.approx(sigma_from_hck, rel=1e-6)
    # Discrimination guard: an older mixed-vintage sigma (5.67051e-8) sits far
    # outside the 1e-6 closure band, so this pins self-consistency, not scale.
    assert abs(5.67051e-8 - sigma_from_hck) / sigma_from_hck > 1e-5
    # Universal / molar gas-constant chain.
    assert phys.Rstar == pytest.approx(1000.0 * phys.k * phys.N_avogadro, rel=1e-12)
    assert phys.Rstar == pytest.approx(1000.0 * phys.R_gas, rel=1e-12)
    # Scale guard: molar gas constant is 8.314 J/mol/K, not 8314 or 0.008314.
    assert 8.0 < phys.R_gas < 8.5
    assert phys.sigma > 0.0


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_molar_masses_agree_between_lookup_table_and_gas_objects():
    """Per-species molar masses match between the dict and the gas objects.

    Molar mass has two representations: the ``molar_mass`` dict in kg/mol used
    by the height and adiabat integrators, and the per-gas objects'
    ``MolecularWeight`` in g/mol used by the specific-gas-constant path
    R = Rstar/MolecularWeight. For every shared species the two must describe
    the same molecule, molar_mass[X]*1000 == gas.MolecularWeight; a drift would
    give one species two different masses and two different gas constants. CO's
    heat-capacity ratio is pinned to the diatomic 7/5 value as an independent
    anchor, and an integer-rounded H2 mass is rejected as the edge case.
    """
    pairs = {
        'H2O': phys.H2O,
        'CO2': phys.CO2,
        'CH4': phys.CH4,
        'O2': phys.O2,
        'H2': phys.H2,
        'N2': phys.N2,
        'CO': phys.CO,
    }
    for formula, gasobj in pairs.items():
        table_g_per_mol = phys.molar_mass[formula] * 1000.0
        assert gasobj.MolecularWeight == pytest.approx(table_g_per_mol, rel=1e-4)
    # Water is the tightest anchor: both representations read 18.01528 g/mol.
    assert phys.H2O.MolecularWeight == pytest.approx(18.01528, rel=1e-7)
    # Discrimination: an integer-rounded H2 mass (2.0) misses the tabulated
    # 2.01588 by ~0.8%, far past the 1e-4 cross-table consistency band.
    assert phys.H2.MolecularWeight == pytest.approx(2.01588, rel=1e-5)
    assert abs(2.0 - phys.H2.MolecularWeight) / phys.H2.MolecularWeight > 1e-4
    # CO is diatomic: gamma = 7/5, not a monatomic 5/3 or an absurd value.
    assert phys.CO.gamma == pytest.approx(1.4, rel=1e-6)
    assert 1.3 < phys.CO.gamma < 1.5


@pytest.mark.physics_invariant
def test_satvps_function_numeric_construction_ignores_water_lookup():
    """A numerically built satvps_function never dispatches to the water tables.

    ``satvps_function`` is built either from a gas object or from bare numbers
    (T0, e0, MolecularWeight, LatentHeat). In the numeric case it carries no gas
    identity, so it must fall through to the closed-form Clausius-Clapeyron law
    on every call, whether or not water_lookup is requested; asking for the
    steam-table path on such an instance must neither raise nor change the
    result. This gas-less request is the error-contract / edge path.
    """
    T0, e0, mol_weight, latent = 273.16, 611.0, 18.0, 2.5e6
    sat = phys.satvps_function(T0, e0, mol_weight, latent)
    # A numeric construction attaches no gas identity.
    assert sat.gas is None
    baseline = sat(300.0)
    # The Clausius-Clapeyron value is recovered directly.
    assert baseline == pytest.approx(phys.satvps(300.0, T0, e0, mol_weight, latent), rel=1e-12)
    # Error-contract path: requesting the water lookup on a gas-less instance
    # must neither raise nor alter the returned pressure.
    with_lookup = sat(300.0, water_lookup=True)
    assert with_lookup == pytest.approx(baseline, rel=1e-12)
    # Physical checks: positive, and 27 K above the reference it exceeds e0
    # (Clausius-Clapeyron sign).
    assert baseline > 0.0
    assert baseline > e0


@pytest.mark.physics_invariant
def test_planck_temperature_derivative_positive_and_matches_finite_difference():
    """dB(nu,T) reproduces the analytic dB/dT of the Planck function and stays positive.

    ``dB`` is the closed-form temperature derivative of the Planck radiance
    B(nu,T). At fixed frequency a black body brightens as it warms, so dB/dT is
    strictly positive; feeding it back against a centred finite difference of the
    independently coded ``B`` pins the analytic derivative formula. Two
    temperatures 300 K and 1500 K are used so the strong T dependence of the Wien
    tail is resolved and a dropped prefactor cannot hide inside a single point.
    """
    nu = 5.0e13  # ~6 micron, thermal-IR
    d_cool = phys.dB(nu, 300.0)
    d_warm = phys.dB(nu, 1500.0)
    # Positivity: warming a black body raises its radiance at every frequency.
    assert d_cool > 0.0
    assert d_warm > 0.0
    # The analytic derivative must match a centred finite difference of B to the
    # order of the step; a mis-derived prefactor would break this closure.
    for temp in (300.0, 1500.0):
        step = 1.0e-3 * temp
        fd = (phys.B(nu, temp + step) - phys.B(nu, temp - step)) / (2.0 * step)
        assert phys.dB(nu, temp) == pytest.approx(fd, rel=1e-4)
    # Discrimination: the derivative itself rises steeply with T at this
    # frequency (a factor of ~30 from 300 to 1500 K), far from the flat response
    # a wrong exponent would give.
    assert d_warm / d_cool > 10.0


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_heymsfield_saturation_matches_one_atm_and_rises_with_temperature():
    """The Heymsfield water-saturation fit gives ~1 atm at the steam point.

    ``satvpw_Heymsfield`` is an alternate saturation-vapour-pressure formula over
    liquid water. At the 373.16 K steam point it returns 1.0132e5 Pa, one standard
    atmosphere to four figures, and it must track the Smithsonian ``satvpw`` curve
    to better than one percent across the liquid range. The curve rises
    monotonically with temperature (Clausius-Clapeyron), and stays strictly
    positive; the sub-triple-point limit input is the physical edge.
    """
    p_boil = phys.satvpw_Heymsfield(373.16)
    # Published anchor: saturation over liquid water at the steam point is one
    # standard atmosphere (101325 Pa); pins the absolute scale.
    assert p_boil == pytest.approx(101325.0, rel=1e-3)
    # Cross-check against the Smithsonian formula: the two water fits agree to
    # better than 1 percent, a genuine cross-implementation constraint.
    assert p_boil == pytest.approx(phys.satvpw(373.16), rel=1e-2)
    # Scale guard: within 0.5 percent of one atmosphere, ruling out a Pa/hPa slip
    # (would give ~1013) or the un-converted dyn/cm^2 value (~1.01e6).
    assert 1.00e5 < p_boil < 1.02e5
    assert p_boil > 0.0
    # Clausius-Clapeyron monotonicity: a cooler parcel holds less vapour.
    assert phys.satvpw_Heymsfield(300.0) < p_boil


@pytest.mark.physics_invariant
def test_satvpg_blends_ice_and_water_across_freezing():
    """satvpg switches ice-to-water saturation and blends across the freezing band.

    ``satvpg`` is the GFDL saturation vapour pressure: it uses the ice curve
    below -20 C, the pure liquid curve above 0 C, and a linear blend of the two
    in between. This test exercises all three branches. Above freezing the value
    must equal ``satvpw`` and well below it must equal ``satvpi``; inside the
    blend band the value must lie strictly between the ice and water curves, which
    a dropped blend weight would violate. All values stay positive and the curve
    rises with temperature.
    """
    # Warm branch (T-273.16 > 0): pure liquid-water saturation.
    assert phys.satvpg(300.0) == pytest.approx(phys.satvpw(300.0), rel=1e-12)
    # Cold branch (T-273.16 < -20): pure ice saturation.
    assert phys.satvpg(250.0) == pytest.approx(phys.satvpi(250.0), rel=1e-12)
    # Blend branch (-20 <= T-273.16 <= 0): strictly between ice and water. This
    # is the discriminating case; a collapsed blend would snap to one endpoint.
    blend = phys.satvpg(260.0)
    assert phys.satvpi(260.0) < blend < phys.satvpw(260.0)
    # Positivity and monotonic rise from the ice regime to the warm regime.
    assert phys.satvpg(250.0) > 0.0
    assert phys.satvpg(250.0) < blend < phys.satvpg(300.0)


@pytest.mark.physics_invariant
def test_satvps_function_ice_and_liquid_flags_select_latent_heat():
    """The ice and liquid flags pick the sublimation and vaporization latent heats.

    Built from a gas object, ``satvps_function`` accepts an explicit phase flag.
    The 'ice' flag must load the latent heat of sublimation and 'liquid' the
    latent heat of vaporization. Because sublimation is the larger latent heat,
    the ice branch gives the steeper Clausius-Clapeyron curve and therefore a
    lower saturation pressure below the triple point; that inequality
    discriminates a swapped flag. Both saturation pressures stay positive at the
    sub-triple-point edge.
    """
    sat_ice = phys.satvps_function(phys.H2O, 'ice')
    sat_liq = phys.satvps_function(phys.H2O, 'liquid')
    # The flags select the two distinct latent heats of the water object.
    assert sat_ice.L == pytest.approx(phys.H2O.L_sublimation, rel=1e-12)
    assert sat_liq.L == pytest.approx(phys.H2O.L_vaporization, rel=1e-12)
    # Sublimation exceeds vaporization, so the ice curve is the steeper one.
    assert phys.H2O.L_sublimation > phys.H2O.L_vaporization
    # Below the triple point the steeper ice curve sits below the liquid curve; a
    # swapped flag would invert this. Pin the ice value against the closed form.
    cold = 250.0
    assert sat_ice(cold) == pytest.approx(
        phys.satvps(
            cold,
            phys.H2O.TriplePointT,
            phys.H2O.TriplePointP,
            phys.H2O.MolecularWeight,
            phys.H2O.L_sublimation,
        ),
        rel=1e-12,
    )
    assert sat_ice(cold) < sat_liq(cold)
    assert sat_ice(cold) > 0.0


@pytest.mark.physics_invariant
def test_satvps_function_water_lookup_dispatches_by_triple_point_and_flag():
    """Water lookup uses steam tables above the triple point, Clausius-Clapeyron below.

    For a water gas object with water_lookup requested, ``satvps_function``
    returns the tabulated IAPWS saturation pressure above the triple point, and
    below it falls back to the closed-form Clausius-Clapeyron law with a latent
    heat chosen by the phase flag: the sublimation heat for the default 'switch'
    instance and the tabulated vaporization heat for an explicit 'liquid'
    instance. The two sub-triple-point paths must give different pressures (the
    error-contract that a dropped flag branch would collapse), and the table
    branch must reproduce the water_tables lookup exactly.
    """
    sat_liq = phys.satvps_function(phys.H2O, 'liquid')
    warm = 300.0  # above the 273.16 K triple point -> steam-table branch
    assert sat_liq(warm, water_lookup=True) == pytest.approx(wt.lookup('psat', warm), rel=1e-12)
    # The tabulated value differs from the bare closed form, so the lookup branch
    # is a genuine special case and not an accidental match.
    closed = phys.satvps(
        warm,
        phys.H2O.TriplePointT,
        phys.H2O.TriplePointP,
        phys.H2O.MolecularWeight,
        phys.H2O.L_vaporization,
    )
    assert abs(sat_liq(warm, water_lookup=True) - closed) / closed > 1e-3

    cold = 250.0  # below the triple point -> Clausius-Clapeyron fallback
    liq_cold = sat_liq(cold, water_lookup=True)
    assert liq_cold == pytest.approx(
        phys.satvps(
            cold,
            phys.H2O.TriplePointT,
            phys.H2O.TriplePointP,
            phys.H2O.MolecularWeight,
            wt.L_vap[0],
        ),
        rel=1e-12,
    )
    sat_switch = phys.satvps_function(phys.H2O)  # default flag becomes 'switch'
    assert sat_switch.iceFlag == 'switch'
    switch_cold = sat_switch(cold, water_lookup=True)
    assert switch_cold == pytest.approx(
        phys.satvps(
            cold,
            phys.H2O.TriplePointT,
            phys.H2O.TriplePointP,
            phys.H2O.MolecularWeight,
            phys.H2O.L_sublimation,
        ),
        rel=1e-12,
    )
    # The switch instance uses the larger sublimation heat, so its sub-triple
    # saturation is strictly below the liquid instance; both stay positive.
    assert switch_cold < liq_cold
    assert switch_cold > 0.0


@pytest.mark.physics_invariant
def test_moist_adiabat_water_air_profile_cools_upward_and_stays_physical():
    """The single-condensable moist adiabat cools upward and stays thermodynamically valid.

    ``MoistAdiabat`` integrates the moist adiabat for a condensable (water) mixed
    with a non-condensing background (air). Driven from an Earth-like 300 K, 1 bar
    surface up to a 1000 Pa top, the returned column must start at the surface
    temperature, cool monotonically with decreasing pressure, keep pressure and
    temperature strictly positive and finite at every level, and keep the
    condensable molar concentration inside [0, 1]. The optional pressure-grid
    interpolation must return the requested grid and preserve the downward cooling
    trend; this is the alternate output-path edge case.
    """
    m = phys.MoistAdiabat(phys.H2O, phys.air)
    press, temp, molar_con, mass_con = m(1.0e5, 300.0, 1000.0)
    # Boundary condition: the profile starts at the prescribed surface temperature.
    assert temp[0] == pytest.approx(300.0, rel=1e-9)
    # Pressure decreases and temperature cools monotonically upward.
    assert np.all(np.diff(press) < 0.0)
    assert np.all(np.diff(temp) <= 0.0)
    # Positivity and finiteness across the whole column.
    assert np.all(temp > 0.0)
    assert np.all(press > 0.0)
    assert np.all(np.isfinite(temp)) and np.all(np.isfinite(press))
    # The condensable molar concentration is a physical fraction.
    assert molar_con.min() >= 0.0
    assert molar_con.max() <= 1.0
    assert mass_con.max() <= 1.0
    # The top of the integrated column is genuinely colder than the surface, a
    # discriminating drop far larger than any rounding tolerance.
    assert temp[-1] < temp[0] - 50.0

    # Alternate output path: interpolate onto a requested descending grid.
    grid = [5.0e4, 1.0e4, 5.0e3]
    p_grid, t_grid, mc_grid, q_grid = m(1.0e5, 300.0, 1000.0, grid)
    np.testing.assert_allclose(p_grid, grid, rtol=1e-12)
    assert np.all(np.diff(t_grid) < 0.0)  # still cooling upward on the coarse grid
    assert np.all(np.isfinite(t_grid))
