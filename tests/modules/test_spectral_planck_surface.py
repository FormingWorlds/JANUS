"""Tests for src/janus/modules/spectral_planck_surface.py.

spectral_planck_surface.py provides ``surf_Planck_nu``, the per-band surface
Planck emission in W/m^2 used as the lower radiative boundary condition. This
file exercises:

* The band-integrated flux against the grey Stefan-Boltzmann limit
  (1-albedo)*sigma*Ts^4 and its Ts^4 scaling.
* The (1-albedo_s) prefactor: emission vanishes for a perfect reflector and
  scales linearly with (1-albedo), and every per-band term is non-negative.

Invariant families exercised: positivity/boundedness, monotonicity, and
pinned-value-with-discrimination-guard. See docs/How-to/test.md.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from janus.modules.spectral_planck_surface import surf_Planck_nu

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# CODATA Stefan-Boltzmann constant, W/m^2/K^4.
SIGMA_SB = 5.670374419e-8


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_band_integrated_surface_flux_matches_stefan_boltzmann():
    """Summing surf_Planck_nu over a fine grid recovers (1-a)*sigma*Ts^4.

    ``surf_Planck_nu`` returns the per-band surface Planck emission. Summed
    over a wavenumber grid dense enough to resolve the Planck peak it must
    equal the grey Stefan-Boltzmann surface flux (1-albedo)*sigma*Ts^4
    (analytical black/grey-body limit). Halving Ts must drop the integrated
    flux by 2^4=16, which discriminates a wrong Planck exponent.
    """
    ts, albedo = 1500.0, 0.0
    wavenum = np.linspace(1.0, 30000.0, 60000)  # cm^-1, fine enough to converge
    widths = np.gradient(wavenum)
    atm = SimpleNamespace(band_centres=wavenum, band_widths=widths, ts=ts, albedo_s=albedo)
    total = float(np.sum(surf_Planck_nu(atm)))
    stefan_boltzmann = (1.0 - albedo) * SIGMA_SB * ts**4
    # Grid converges to ~4e-6; 2e-3 leaves generous headroom.
    assert total == pytest.approx(stefan_boltzmann, rel=2e-3)

    # T^4 scaling guard: halving Ts drops the integrated flux by 2^4.
    atm_cool = SimpleNamespace(
        band_centres=wavenum, band_widths=widths, ts=ts / 2.0, albedo_s=albedo
    )
    total_cool = float(np.sum(surf_Planck_nu(atm_cool)))
    assert total / total_cool == pytest.approx(16.0, rel=1e-2)


@pytest.mark.physics_invariant
def test_surface_flux_scales_with_albedo_and_vanishes_at_unit_albedo():
    """Surface emission scales as (1-albedo) and is zero for a perfect mirror.

    The (1-albedo_s) prefactor means a perfectly reflecting surface
    (albedo=1) emits nothing, emission is strictly positive for albedo<1, and
    an albedo of 0.5 halves the black-surface emission. Every per-band
    contribution must be non-negative.
    """
    ts = 1000.0
    wavenum = np.linspace(1.0, 20000.0, 4000)
    widths = np.gradient(wavenum)

    def flux(albedo):
        atm = SimpleNamespace(band_centres=wavenum, band_widths=widths, ts=ts, albedo_s=albedo)
        return surf_Planck_nu(atm)

    b_black = flux(0.0)
    b_grey = flux(0.5)
    b_mirror = flux(1.0)
    assert np.all(b_black >= 0.0)  # per-band non-negativity
    assert np.sum(b_black) > 0.0
    # Perfect reflector emits nothing (error-contract limit).
    assert np.sum(b_mirror) == pytest.approx(0.0, abs=1e-9)
    # (1-albedo) linearity: albedo 0.5 halves the black-surface emission.
    assert np.sum(b_grey) == pytest.approx(0.5 * np.sum(b_black), rel=1e-12)
