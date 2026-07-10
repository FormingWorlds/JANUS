"""Tests for src/janus/modules/set_stratosphere.py.

set_stratosphere.py rewrites the isothermal stratosphere above the tropopause:
it pins temperature to the tropopause value, resets the volatile mixing ratios
to their tropopause partition, and rebuilds the layer heat capacity and mean
molar mass. The heat capacity is accumulated over every volatile and then
normalised once by the gas-phase molar concentration; the mean molar mass is
accumulated over the volatile partial pressures and normalised once by the total
pressure. Both normalisations apply to the whole-column sum, not per volatile.
This file exercises:

* Single whole-column normalisation of cp and mu: the result matches the
  gas-weighted sum divided once, and a per-volatile (repeated) normalisation
  would drive mu far below any physical molar mass.
* The mean molar mass landing in the CO2/N2 physical band.
* The below-tropopause guard return and the minimum-temperature clip.

Invariant families exercised: conservation/aggregation symmetry (whole-column
normalisation), positivity/boundedness (physical molar-mass band), and
pinned-value-with-discrimination-guard. See docs/How-to/test.md.
"""

from types import SimpleNamespace

import numpy as np
import pytest

import janus.utils.GeneralAdiabat as ga
import janus.utils.phys as phys
from janus.modules.set_stratosphere import set_stratosphere

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

VOLS = ['CO2', 'N2']


def _build_atm(trpp_idx, trpp_tmp=150.0, min_t=10.0):
    """Assemble a minimal two-volatile atmosphere column for set_stratosphere.

    Index 0 is the top of the column (lowest pressure); the tropopause sits at
    ``trpp_idx``. Mixing ratios are gas-dominated with a condensate fraction so
    the stratosphere rebuild has non-trivial cp and mu to accumulate.
    """
    nlev = 6
    press = np.array([1.0e2, 3.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4])
    pressl = np.array([5.0e1, 1.5e2, 5.0e2, 1.5e3, 5.0e3, 1.5e4, 4.5e4])
    x_gas = {'CO2': np.full(nlev, 0.30), 'N2': np.full(nlev, 0.50)}
    x_cond = {'CO2': np.full(nlev, 0.15), 'N2': np.full(nlev, 0.05)}
    return SimpleNamespace(
        trppidx=trpp_idx,
        trppP=press[trpp_idx] if trpp_idx >= 0 else press[0],
        trppT=trpp_tmp,
        p=press.copy(),
        pl=pressl.copy(),
        tmp=np.full(nlev, 260.0),
        tmpl=np.full(nlev + 1, 260.0),
        minT=min_t,
        vol_list={v: 1.0 for v in VOLS},
        x_cond={v: x_cond[v].copy() for v in VOLS},
        x_gas={v: x_gas[v].copy() for v in VOLS},
        p_vol={v: x_gas[v].copy() * press for v in VOLS},
        xc=np.full(nlev, 0.20),
        xv=np.full(nlev, 0.80),
        xd=np.zeros(nlev),
        cp=np.full(nlev, 999.0),  # sentinel: overwritten only if the loop runs
        mu=np.full(nlev, 999.0),
        water_lookup=False,
        alpha_cloud=1.0,
    )


@pytest.mark.physics_invariant
def test_stratosphere_normalises_cp_and_mu_once_over_the_whole_column():
    """cp and mu are normalised once over all volatiles, not once per volatile.

    After the rebuild each stratosphere level must satisfy cp = (sum over
    volatiles of gas and retained-condensate heat capacity) / (xd + xv) and
    mu = (sum over volatiles of molar_mass * partial pressure) / p_total, with a
    single division. Recomputing from the final arrays reproduces the routine to
    machine precision. A per-volatile normalisation would divide the running sum
    once per species: mu would pick up an extra factor of the total pressure and
    collapse to order 1e-4 kg/mol, far below any molecular weight, and cp on a
    layer whose gas concentration differs from unity would shift by more than a
    J/mol/K.
    """
    trpp_idx = 3
    atm = _build_atm(trpp_idx)
    set_stratosphere(atm)

    worst_cp_gap = 0.0
    for idx in range(trpp_idx):
        temp = atm.tmp[idx]
        # Stratosphere temperature is pinned to the tropopause value.
        assert temp == pytest.approx(150.0, rel=1e-12)
        conc = atm.xd[idx] + atm.xv[idx]
        cp_sum = sum(
            atm.x_gas[v][idx] * ga.cpv(v, temp)
            + atm.x_cond[v][idx] * ga.cp_cond(v, temp) * atm.alpha_cloud
            for v in VOLS
        )
        mu_sum = sum(phys.molar_mass[v] * atm.p_vol[v][idx] for v in VOLS)
        cp_single = cp_sum / conc
        mu_single = mu_sum / atm.p[idx]

        # Whole-column single normalisation reproduces the routine exactly.
        assert atm.cp[idx] == pytest.approx(cp_single, rel=1e-12)
        assert atm.mu[idx] == pytest.approx(mu_single, rel=1e-12)
        # Physical band: mean molar mass sits between the light and heavy end of
        # a CO2 (0.044) / N2 (0.028) mix, well above the 1e-4 a repeated
        # normalisation would produce.
        assert 0.02 < atm.mu[idx] < 0.05
        assert atm.mu[idx] / atm.p[idx] < 1.0e-3  # the per-volatile-normalised value
        # Track the layer where the gas concentration is furthest from unity, so
        # the doubled cp normalisation is resolvable.
        if abs(conc - 1.0) > abs(worst_cp_gap):
            worst_cp_gap = conc - 1.0
            worst = (idx, atm.cp[idx], cp_single / conc)

    # On the layer with the largest concentration offset, a second cp division
    # differs from the routine by more than a J/mol/K (discrimination guard).
    assert abs(worst_cp_gap) > 0.1
    idx_w, cp_routine, cp_double = worst
    assert abs(cp_routine - cp_double) > 1.0


@pytest.mark.physics_invariant
def test_stratosphere_guard_returns_below_tropopause_and_clips_to_minimum():
    """A negative tropopause index is a no-op, and the minT clip raises the floor.

    When ``trppidx`` is negative the column has no resolved tropopause, so the
    routine returns immediately and leaves cp and mu at their sentinels (the
    guard return / error contract). With a valid tropopause the isothermal value
    is written first and then clipped up to ``minT``: a tropopause temperature
    below the floor must be replaced by the floor, both on standard and on
    staggered nodes.
    """
    atm_guard = _build_atm(-1)
    set_stratosphere(atm_guard)
    # Guard return leaves the accumulators untouched at their sentinel.
    np.testing.assert_allclose(atm_guard.cp, 999.0, rtol=0, atol=0)
    np.testing.assert_allclose(atm_guard.mu, 999.0, rtol=0, atol=0)

    # Empty adjustment loop (trpp_idx=0) isolates the temperature clip: the
    # tropopause value 150 K is raised to the 180 K floor at the top node only.
    atm_clip = _build_atm(0, trpp_tmp=150.0, min_t=180.0)
    set_stratosphere(atm_clip)
    assert atm_clip.tmp[0] == pytest.approx(180.0, rel=1e-12)
    assert atm_clip.tmpl[0] == pytest.approx(180.0, rel=1e-12)
    # Levels below the tropopause pressure are untouched by the isothermal write
    # and already sit above the floor.
    assert atm_clip.tmp[-1] == pytest.approx(260.0, rel=1e-12)
    assert np.all(atm_clip.tmp >= 180.0)
