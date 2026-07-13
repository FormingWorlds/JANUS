"""Tests for src/janus/modules/plot_emission_spectrum.py.

Exercises the top-of-atmosphere emission-spectrum figure: the per-band flux
normalisation and unit conversion of the emission curve, the Wien-peak
placement and albedo scaling of the surface Planck overlay, the optional
stellar-flux and band-edge overlays, and the degenerate zero-emission column.
See docs/How-to/test.md.
"""

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from janus.modules.plot_emission_spectrum import plot_emission  # noqa: E402

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

NBANDS = 4
NLEV = 3


class _SpectralColumn:
    """Minimal stand-in carrying exactly the attributes plot_emission reads."""

    def __init__(self, ts=300.0, albedo_s=0.1):
        self.nbands = NBANDS
        # Band grid in nm spanning thermal wavelengths.
        self.band_edges = np.array([2.0e3, 6.0e3, 1.2e4, 2.5e4, 5.0e4])
        self.band_centres = 0.5 * (self.band_edges[:-1] + self.band_edges[1:])
        self.band_widths = np.diff(self.band_edges)
        # Plausible spectral fluxes in W m-2 per band, columns are levels.
        base = np.linspace(4.0, 1.0, NBANDS)
        self.LW_spectral_flux_up = np.outer(base, np.linspace(1.0, 0.7, NLEV))
        self.SW_spectral_flux_up = 0.05 * self.LW_spectral_flux_up
        self.SW_spectral_flux_down = 0.5 * self.LW_spectral_flux_up.copy()
        self.ts = ts
        self.albedo_s = albedo_s


@pytest.fixture()
def captured_ax(monkeypatch):
    """Record the axes plot_emission draws on, surviving its internal close."""
    captured = {}
    real_subplots = plt.subplots

    def recording_subplots(*args, **kwargs):
        fig, ax = real_subplots(*args, **kwargs)
        captured['fig'] = fig
        captured['ax'] = ax
        return fig, ax

    monkeypatch.setattr(plt, 'subplots', recording_subplots)
    yield captured
    plt.close('all')


def _expected_emission(atm, level_idx):
    """Per-band upward flux converted to erg s-1 cm-2 nm-1, as the plot does."""
    y = atm.LW_spectral_flux_up[:, level_idx] + atm.SW_spectral_flux_up[:, level_idx]
    return y / atm.band_widths * 1e3


@pytest.mark.physics_invariant
def test_plot_emission_curve_normalisation_and_planck_peak(tmp_path, captured_ax):
    """The emission curve carries the band-width normalisation and the Planck
    overlay peaks at the Wien wavelength.

    The per-band flux is divided by the band width and scaled by 1e3; a
    dropped width division or unit factor rescales the whole curve, so the
    line data is pinned against an independent recomputation. The surface
    overlay for ts = 300 K must peak near the Wien wavelength
    lambda_max = 2.898e6 nm K / 300 K, discriminating a misplaced exponent
    in the Planck function.
    """
    atm = _SpectralColumn(ts=300.0, albedo_s=0.1)
    out = tmp_path / 'emission.png'
    plot_emission(atm, filename=str(out), incoming_solar=False, planck_surface=True)

    ax = captured_ax['ax']
    # Two lines: surface Planck overlay first, then the emission spectrum.
    assert len(ax.lines) == 2
    planck_line, emission_line = ax.lines
    np.testing.assert_allclose(
        emission_line.get_ydata(), _expected_emission(atm, 0), rtol=1e-12
    )
    # Planck overlay: positive, finite, Wien peak within the sampled range.
    yp = planck_line.get_ydata()
    xp = planck_line.get_xdata()
    assert np.all(yp > 0.0)
    assert np.all(np.isfinite(yp))
    wien_nm = 2.898e6 / atm.ts
    peak_nm = xp[np.argmax(yp)]
    # The sampled grid spans the peak; factor 1.5 resolves an exponent slip,
    # which would push the apparent peak to a grid boundary.
    assert wien_nm / 1.5 < peak_nm < wien_nm * 1.5
    assert out.stat().st_size > 0


@pytest.mark.physics_invariant
def test_plot_emission_planck_overlay_scales_with_albedo(tmp_path, captured_ax):
    """The surface overlay scales by (1 - albedo).

    Doubling the reflectance from 0 to 0.5 must halve the overlay
    everywhere; an overlay computed without the albedo factor fails the
    ratio pin.
    """
    out = tmp_path / 'emission_albedo.png'
    atm_dark = _SpectralColumn(albedo_s=0.0)
    plot_emission(atm_dark, filename=str(out), incoming_solar=False)
    y_dark = captured_ax['ax'].lines[0].get_ydata().copy()

    atm_half = _SpectralColumn(albedo_s=0.5)
    plot_emission(atm_half, filename=str(out), incoming_solar=False)
    y_half = captured_ax['ax'].lines[0].get_ydata().copy()

    assert y_dark.shape == y_half.shape
    np.testing.assert_allclose(y_half, 0.5 * y_dark, rtol=1e-12)
    # Sign and scale guard: both overlays strictly positive.
    assert np.all(y_dark > 0.0) and np.all(y_half > 0.0)


def test_plot_emission_optional_overlays_and_axes(tmp_path, captured_ax):
    """Stellar-flux and band-edge overlays add the expected artists.

    With the stellar overlay on and band edges shown, the axes carry the
    solar line plus one axvline per band edge, log-log scaling, and
    wavelength/flux labels.
    """
    atm = _SpectralColumn()
    # The stellar overlay is the level-0 downward shortwave flux with the
    # same band-width normalisation; snapshot it before the call because the
    # overlay is computed on a view of the stored array.
    expected_solar = atm.SW_spectral_flux_down[:, 0].copy() / atm.band_widths * 1e3
    out = tmp_path / 'emission_full.png'
    plot_emission(
        atm,
        filename=str(out),
        incoming_solar=True,
        planck_surface=False,
        show_bands=True,
    )

    ax = captured_ax['ax']
    # solar line + 2*nbands band-edge vlines + emission line.
    assert len(ax.lines) == 1 + 2 * NBANDS + 1
    solar_line = ax.lines[0]
    np.testing.assert_allclose(solar_line.get_ydata(), expected_solar, rtol=1e-12)
    assert ax.get_xscale() == 'log'
    assert ax.get_yscale() == 'log'
    assert 'Wavelength' in ax.get_xlabel()
    assert 'nm' in ax.get_xlabel()
    assert 'Flux density' in ax.get_ylabel()
    assert out.stat().st_size > 0


def test_plot_emission_zero_flux_column_and_level_contract(tmp_path, captured_ax):
    """A zero-emission column renders finite output and a bad level index
    raises.

    Zero upward flux is the limit of a perfectly absorbing column; the
    figure must still be produced with an all-zero emission line. Indexing a
    level outside the stored profile is the documented failure mode of the
    level_idx argument and must raise IndexError rather than plot garbage.
    """
    atm = _SpectralColumn()
    atm.LW_spectral_flux_up = np.zeros((NBANDS, NLEV))
    atm.SW_spectral_flux_up = np.zeros((NBANDS, NLEV))
    out = tmp_path / 'emission_zero.png'
    plot_emission(atm, filename=str(out), incoming_solar=False, planck_surface=False)

    ax = captured_ax['ax']
    np.testing.assert_allclose(ax.lines[0].get_ydata(), np.zeros(NBANDS), atol=0.0)
    assert out.stat().st_size > 0

    with pytest.raises(IndexError):
        plot_emission(
            _SpectralColumn(),
            filename=str(tmp_path / 'bad.png'),
            level_idx=NLEV + 5,
            incoming_solar=False,
            planck_surface=False,
        )
