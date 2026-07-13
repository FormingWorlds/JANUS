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
    # Scale and sign guards: positive, and CO2 per kg scatters about twice
    # as strongly as H2O here: its 5.0x larger molecular cross section
    # outweighs its 2.4x larger molar mass.
    assert h2o['cross_section'] > 0.0 and co2['cross_section'] > 0.0
    assert co2['cross_section'] == pytest.approx(2.03 * h2o['cross_section'], rel=1e-2)


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


SPECTRAL_TOY = """*BLOCK: TYPE =    1: SUBTYPE =    0: VERSION =    0
Specification of spectral intervals
Number of spectral intervals =     2
Band        Lower limit         Upper limit
    1        4.000000000E-07     5.000000000E-07
    2        5.000000000E-07     1.000000000E-06
*END
*BLOCK: TYPE =    4: SUBTYPE =    0: VERSION =    0
Gas absorber data placeholder
*END
"""


def _read_block3_coeffs(path):
    """Parse the (band, coefficient) rows out of a rewritten BLOCK 3."""
    lines = path.read_text().splitlines()
    start = next(i for i, ln in enumerate(lines) if '*BLOCK: TYPE =    3' in ln)
    rows = []
    for ln in lines[start:]:
        parts = ln.split()
        if len(parts) == 2 and parts[0].isdigit():
            rows.append((int(parts[0]), float(parts[1])))
        if '*END' in ln and rows:
            break
    return dict(rows)


@pytest.mark.physics_invariant
def test_rayleigh_coeff_adder_inserts_band_averaged_block3(tmp_path):
    """The spectral-file rewrite inserts a BLOCK 3 with band-averaged
    lambda^-4 coefficients between the band block and the gas block.

    The blue band (400-500 nm) must carry a larger average coefficient than
    the red band (500-1000 nm), the band-1 value must match an independent
    quadrature of the same cross section over the band, the surrounding
    blocks must survive the rewrite in order, and the temporary
    wavelength-band scratch file must be cleaned up.
    """
    quad = pytest.importorskip('scipy.integrate').quad

    sf = tmp_path / 'toy.sf'
    sf.write_text(SPECTRAL_TOY)
    wldummy = tmp_path / 'wl_bands.txt'

    rs.rayleigh_coeff_adder(
        species_list=['co2'],
        mixing_ratio_list=[1.0],
        spectral_file_path=str(sf),
        wavelength_dummy_file_path=str(wldummy),
    )

    coeffs = _read_block3_coeffs(sf)
    assert set(coeffs) == {1, 2}
    assert coeffs[1] > 0.0 and coeffs[2] > 0.0
    # lambda^-4: the shorter-wavelength band average must dominate.
    assert coeffs[1] > coeffs[2]
    # Independent band-1 average: quadrature of the same cross section over
    # 400-500 nm divided by the 100 nm width; pins the averaging and the
    # numeric formatting of the rewritten block.
    info = rs.species_info('co2')
    expected_b1 = quad(rs.cross_section, 4.0e-7, 5.0e-7, args=(info,))[0] / 1.0e-7
    assert coeffs[1] == pytest.approx(expected_b1, rel=1e-6)

    # Structural integrity: band block before BLOCK 3, gas block after it,
    # and the scratch file removed.
    text = sf.read_text()
    assert text.index('TYPE =    1') < text.index('TYPE =    3') < text.index('TYPE =    4')
    assert not wldummy.exists()
    assert not (tmp_path / 'toy.sf_temp_spectral_file').exists()


def test_rayleigh_coeff_adder_replaces_existing_block_and_weights_mixing(tmp_path):
    """An existing BLOCK 3 is replaced, and coefficients scale with the
    molar mixing ratio.

    The stale sentinel coefficient must vanish from the rewritten file, and
    halving the mixing ratio must exactly halve every band coefficient; a
    rewrite that appended instead of replacing, or ignored the weighting,
    fails both pins.
    """
    stale = (
        '*BLOCK: TYPE =    3: SUBTYPE =    0: VERSION =    0\n'
        'Rayleigh mass scatering coefficients at STP: unit m**2/kg\n'
        'Band        Rayleigh coefficient\n'
        '                 m2/kg\n'
        '    1        9.876543210E-01\n'
        '    2        9.876543210E-01\n'
        '*END\n'
    )
    toy = SPECTRAL_TOY.replace('*BLOCK: TYPE =    4', stale + '*BLOCK: TYPE =    4', 1)

    sf_full = tmp_path / 'full.sf'
    sf_full.write_text(toy)
    rs.rayleigh_coeff_adder(
        species_list=['co2'],
        mixing_ratio_list=[1.0],
        spectral_file_path=str(sf_full),
        wavelength_dummy_file_path=str(tmp_path / 'wl1.txt'),
    )
    text = sf_full.read_text()
    assert '9.876543210E-01' not in text
    assert text.count('*BLOCK: TYPE =    3') == 1
    full = _read_block3_coeffs(sf_full)

    sf_half = tmp_path / 'half.sf'
    sf_half.write_text(SPECTRAL_TOY)
    rs.rayleigh_coeff_adder(
        species_list=['co2'],
        mixing_ratio_list=[0.5],
        spectral_file_path=str(sf_half),
        wavelength_dummy_file_path=str(tmp_path / 'wl2.txt'),
    )
    half = _read_block3_coeffs(sf_half)
    for band in (1, 2):
        assert half[band] == pytest.approx(0.5 * full[band], rel=1e-9)
