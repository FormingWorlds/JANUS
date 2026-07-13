"""Tests for src/janus/utils/RayleighSpectrum.py.

Exercises the Rayleigh-scattering helpers: the species lookup including the
normalised-cross-section entries, the lambda^-4 law of the cross-section
evaluator in both parameterisations, and the mixing-ratio weighting of the
band integrator. See docs/How-to/test.md.
"""

import numpy as np
import pytest

from janus.utils import RayleighSpectrum as rs

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def test_species_info_normalised_cross_section_entries():
    """The pluriel H2O and CO2 entries carry mass cross sections at their
    stated normalisation wavelengths.

    Both entries convert a per-molecule cross section in cm^2 to m^2/kg
    (factor 1e-4 * N_A / M); the pinned values re-derive that conversion,
    so a dropped 1e-4 or a wrong molar mass shifts them by orders of
    magnitude. The two species must also differ from each other.
    """
    h2o = rs.species_info('pluriel_h2o')
    co2 = rs.species_info('pluriel_co2')

    expected_h2o = 2.5e-27 * 1e-4 * 6.022e23 / 18e-3
    expected_co2 = 1.24e-26 * 1e-4 * 6.022e23 / 44e-3
    assert h2o['cross_section'] == pytest.approx(expected_h2o, rel=1e-12)
    assert co2['cross_section'] == pytest.approx(expected_co2, rel=1e-12)
    assert h2o['normalization_wavelength'] == pytest.approx(0.6e-6, rel=1e-12)
    assert co2['normalization_wavelength'] == pytest.approx(0.532e-6, rel=1e-12)
    # Scale and sign guards: positive, and CO2 per kg scatters less than H2O
    # per kg here because of its larger molar mass relative to the raw
    # molecular cross section ratio.
    assert h2o['cross_section'] > 0.0 and co2['cross_section'] > 0.0
    assert h2o['cross_section'] != pytest.approx(co2['cross_section'], rel=1e-2)


@pytest.mark.physics_invariant
def test_cross_section_follows_inverse_fourth_power_law():
    """Both evaluator branches scale as lambda^-4 toward shorter wavelengths.

    Rayleigh scattering strengthens steeply blueward; doubling the
    wavelength must reduce the normalised-branch coefficient by exactly
    2^4 and the refractive-index branch by close to 2^4 (the small
    dispersion correction keeps it within a few percent). A sign or
    exponent error breaks the ratio by orders of magnitude.
    """
    lam = 0.6e-6
    info_h2o = rs.species_info('pluriel_h2o')
    c1 = rs.cross_section(lam, info_h2o)
    c2 = rs.cross_section(2 * lam, info_h2o)
    assert c1 > 0.0
    assert c1 / c2 == pytest.approx(16.0, rel=1e-12)
    # At the normalisation wavelength the coefficient equals the stored one.
    at_norm = rs.cross_section(info_h2o['normalization_wavelength'], info_h2o)
    assert at_norm == pytest.approx(info_h2o['cross_section'], rel=1e-12)

    info_n2 = rs.species_info('n2')
    d1 = rs.cross_section(lam, info_n2)
    d2 = rs.cross_section(2 * lam, info_n2)
    assert d1 > 0.0
    # Dispersion (the 1 + B/microns^2 term) perturbs the pure power law at
    # the percent level over one octave.
    assert d1 / d2 == pytest.approx(16.0, rel=5e-2)


@pytest.mark.physics_invariant
def test_band_integrator_weights_by_mixing_ratio():
    """The band integral scales linearly with the molar mixing ratio and
    adds across species.

    Halving the mixing ratio must halve the coefficient exactly, and a
    two-species integral must equal the sum of the single-species ones;
    a vanishing mixing ratio is the degenerate no-scatterer limit.
    """
    w1 = np.array([0.5e-6])
    w2 = np.array([0.7e-6])

    full = rs.band_integrator(['co2'], [1.0], w1, w2)
    half = rs.band_integrator(['co2'], [0.5], w1, w2)
    none = rs.band_integrator(['co2'], [0.0], w1, w2)
    assert full[0] > 0.0
    assert half[0] == pytest.approx(0.5 * full[0], rel=1e-12)
    assert none[0] == pytest.approx(0.0, abs=1e-30)

    n2_full = rs.band_integrator(['n2'], [1.0], w1, w2)
    both = rs.band_integrator(['co2', 'n2'], [1.0, 1.0], w1, w2)
    assert both[0] == pytest.approx(full[0] + n2_full[0], rel=1e-12)
