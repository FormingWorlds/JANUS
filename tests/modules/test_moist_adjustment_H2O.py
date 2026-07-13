"""Tests for src/janus/modules/moist_adjustment_H2O.py.

moist_adjustment_H2O.py provides ``moist_adj``, a constant-latent-heat moist
convective adjustment that relaxes supersaturated levels toward the local H2O
dew point. This file exercises:

* The two degenerate limits that produce no adjustment: a water-free column
  (early return) and a column everywhere above the dew point.
* A supersaturated (cold) column, whose relaxation wiring (Tdew-T)/tau and 1/tau
  scaling toward the local dew point are verified.
* The dew-point relaxation target itself, anchored to the water triple point and
  the constant-latent-heat Clausius-Clapeyron inversion.

Invariant families exercised: conservation-style limit behaviour,
positivity/boundedness, monotonicity, and pinned-value-with-discrimination-guard
against the water triple point. See docs/How-to/test.md.
"""

import math
from types import SimpleNamespace

import numpy as np
import pytest

import janus.utils.GeneralAdiabat as ga
import janus.utils.phys as phys
from janus.modules.moist_adjustment_H2O import moist_adj

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.physics_invariant
def test_moist_adjustment_zero_without_water_or_when_subsaturated():
    """moist_adj returns no adjustment absent water or above the dew point.

    Two degenerate limits must both yield an all-zero tendency: a column with
    no H2O in vol_list (the early return), and a warm column already above the
    local dew point everywhere (no level triggers adjustment). The warm-column
    case discriminates a broken "always adjust" path that would push levels
    down to the dew point and produce a spurious negative tendency.
    """
    nlev = 10
    press = np.logspace(5, 3, nlev)
    # No water at all -> early return of zeros.
    atm_dry = SimpleNamespace(
        p=press.copy(),
        tmp=np.full(nlev, 300.0),
        vol_list={'CO2': 1.0},
        p_vol={'CO2': press.copy()},
    )
    dt_dry = moist_adj(atm_dry, conv_timescale=1.0)
    assert dt_dry.shape == (nlev,)
    np.testing.assert_array_equal(dt_dry, np.zeros(nlev))

    # Water present but the column is far too warm to condense anywhere.
    pp_h2o = press.copy()
    atm_warm = SimpleNamespace(
        p=press.copy(),
        tmp=np.full(nlev, 2000.0),
        vol_list={'H2O': 1.0},
        p_vol={'H2O': pp_h2o},
    )
    dt_warm = moist_adj(atm_warm, conv_timescale=1.0)
    # Dew point at every level is far below the 2000 K column, so nothing
    # condenses; a broken "always adjust" path would yield a negative tendency.
    max_tdew = max(ga.Tdew('H2O', float(pp)) for pp in pp_h2o)
    assert max_tdew < 2000.0
    np.testing.assert_array_equal(dt_warm, np.zeros(nlev))


@pytest.mark.physics_invariant
def test_supersaturated_column_relaxes_toward_local_dew_point():
    """A supersaturated cold column relaxes toward the local H2O dew point.

    When a level sits below the local H2O dew point, moist_adj raises its
    temperature to ga.Tdew('H2O', pp_h2o), so the tendency is (Tdew-T)/tau. This
    test verifies the relaxation wiring: the tendency on the adjusted interior
    levels equals (Tdew-T)/tau, latent heating warms the cold column (dT>0), and
    the tendency scales as 1/tau, discriminating a wrong timescale factor. The
    dew-point target ga.Tdew is anchored to the water triple point separately in
    test_dew_point_target_matches_clausius_clapeyron_reference.
    """
    nlev = 6
    press = np.logspace(5, 3, nlev)
    pp_h2o = press.copy()
    t_cold = 150.0  # below the dew point at every level's partial pressure
    tau = 100.0
    atm = SimpleNamespace(
        p=press.copy(),
        tmp=np.full(nlev, t_cold),
        vol_list={'H2O': 1.0},
        p_vol={'H2O': pp_h2o},
    )
    dt_conv = moist_adj(atm, conv_timescale=tau)

    # Dew-point relaxation target. The downward pass adjusts levels 0..nlev-2;
    # the top level is outside that range, so restrict the pin to the interior
    # levels where the physics applies.
    tdew = np.array([ga.Tdew('H2O', float(pp)) for pp in pp_h2o])
    expected_interior = (tdew[:-1] - t_cold) / tau
    np.testing.assert_allclose(dt_conv[:-1], expected_interior, rtol=1e-9)
    # Latent heat release warms a cold column everywhere it condenses.
    assert np.all(dt_conv[:-1] > 0.0)

    # 1/tau scaling: doubling the convective timescale halves the tendency.
    atm_slow = SimpleNamespace(
        p=press.copy(),
        tmp=np.full(nlev, t_cold),
        vol_list={'H2O': 1.0},
        p_vol={'H2O': pp_h2o},
    )
    dt_slow = moist_adj(atm_slow, conv_timescale=2.0 * tau)
    np.testing.assert_allclose(dt_slow[:-1], dt_conv[:-1] / 2.0, rtol=1e-9)


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_dew_point_target_matches_clausius_clapeyron_reference():
    """The moist-adjustment dew-point target reproduces the water triple point.

    moist_adj relaxes each supersaturated level to ga.Tdew('H2O', pp), so the
    module's physics is anchored by the dew-point curve. At the triple-point
    pressure 611.657 Pa the dew point must equal the water triple-point
    temperature 273.15 K (the published thermodynamic fixed point); this pins
    the reference. Because the exponent term vanishes there, that single point
    is insensitive to the latent heat, so an off-reference round-trip is added:
    ga.Tdew is the closed-form inverse of the constant-latent-heat
    Clausius-Clapeyron relation p(T)=pref*exp((L/R)(1/Tref-1/T)); feeding a
    hand-computed p(T) back through ga.Tdew must recover T, and a doubled latent
    heat must shift the recovered temperature well past tolerance.
    """
    # Published anchor: water triple point, 273.16 K at 611.657 Pa. The routine
    # uses 273.15 K as its p<=p_triple reference temperature; pin within 0.05 K.
    assert ga.Tdew('H2O', 611.657) == pytest.approx(273.15, abs=0.05)

    # Off-reference round-trip that depends on the latent heat. The p>p_triple
    # branch of ga.Tdew anchors on the 1 atm reference pressure at Tref=373.15 K,
    # so ga.Tdew is the exact inverse of the constant-L Clausius-Clapeyron curve
    # p(T)=pref*exp((L/R)(1/Tref-1/T)) with pref=1e5; the round-trip then closes
    # to machine precision. t_target=340 K keeps both the forward and the
    # doubled-L pressures above the triple point, inside this branch.
    Tref = 373.15
    pref = 1e5  # 1 atm reference pressure at Tref, the p>p_triple branch anchor
    latent_molar = phys.water.L_vaporization * phys.water.MolecularWeight * 1e-3
    coeff = latent_molar / phys.R_gas
    t_target = 340.0
    p_forward = pref * math.exp(coeff * (1.0 / Tref - 1.0 / t_target))
    assert ga.Tdew('H2O', p_forward) == pytest.approx(t_target, rel=1e-9)
    # Latent-heat guard: feeding the same pressure but with a doubled coefficient
    # recovers a temperature tens of kelvin away, far past the round-trip
    # tolerance, so the anchor actually pins L rather than being degenerate.
    p_wrong_l = pref * math.exp(2.0 * coeff * (1.0 / Tref - 1.0 / t_target))
    assert abs(ga.Tdew('H2O', p_wrong_l) - t_target) > 1.0

    # Tie the published anchor to the module output: a cold column with one
    # interior level at the triple-point pressure relaxes that level's target to
    # 273.15 K. press descends so index 1 (611.657 Pa) is an interior level the
    # downward pass adjusts; the top node (index -1) is left untouched.
    pp_triple = 611.657
    press = np.array([1.0e5, pp_triple, 1.0e3])
    t_cold, tau = 200.0, 50.0
    atm = SimpleNamespace(
        p=press.copy(),
        tmp=np.full(3, t_cold),
        vol_list={'H2O': 1.0},
        p_vol={'H2O': press.copy()},
    )
    dt_conv = moist_adj(atm, conv_timescale=tau)
    target_triple = t_cold + dt_conv[1] * tau  # reconstructed relaxation target
    assert target_triple == pytest.approx(273.15, abs=0.1)


@pytest.mark.physics_invariant
def test_moist_adjustment_skips_dry_levels_and_completes_single_pass():
    """moist_adj skips near-vacuum H2O levels and runs a full non-breaking pass.

    A level whose H2O partial pressure is below 1e-10 Pa carries no condensable,
    so both the downward and the upward sweep must skip it rather than call the
    dew-point inversion on a vanishing pressure. With the number of passes set to
    one and a genuinely supersaturated interior, the sweep performs an adjustment
    and therefore does not take the early-exit break, exercising the full single
    pass to completion. The dry top level (index 0) must be left exactly at its
    original temperature, the discriminating check that the skip fired.
    """
    press = np.logspace(5, 3, 4)
    # Index 0 is dry (partial pressure below the 1e-10 Pa floor); the interior
    # levels hold water and sit below their local dew point (cold column).
    pp_h2o = np.array([0.0, 100.0, 200.0, 300.0])
    t_cold = 120.0
    tau = 1000.0
    atm = SimpleNamespace(
        p=press.copy(),
        tmp=np.full(4, t_cold),
        vol_list={'H2O': 1.0},
        p_vol={'H2O': pp_h2o.copy()},
    )
    dt_conv = moist_adj(atm, conv_timescale=tau, nb_convsteps=1)
    assert dt_conv.shape == (4,)
    # The dry level carries no water, so it receives exactly zero tendency.
    assert dt_conv[0] == pytest.approx(0.0, abs=1e-30)
    # The wet interior levels the downward sweep reaches (indices 1 and 2) are
    # relaxed toward their local dew point, so their tendency is strictly warming.
    tdew_interior = np.array([ga.Tdew('H2O', float(pp_h2o[i])) for i in (1, 2)])
    expected = (tdew_interior - t_cold) / tau
    np.testing.assert_allclose(dt_conv[1:3], expected, rtol=1e-9)
    assert np.all(dt_conv[1:3] > 0.0)
