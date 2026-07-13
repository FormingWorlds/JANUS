"""Tests for src/janus/modules/plot_flux_balance.py.

Exercises the flux-balance figure builder: line content against the input
flux arrays (including the sign convention that downward fluxes are drawn
negated), axis scales and orientation for a pressure column, legend labels,
and the degenerate all-zero-flux column. See docs/How-to/test.md.
"""

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from janus.modules.plot_flux_balance import plot_fluxes  # noqa: E402

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


class _FluxColumn:
    """Minimal stand-in carrying exactly the attributes plot_fluxes reads."""

    def __init__(self, nlev=6):
        # Pressure levels in Pa, monotone decreasing upward like a real column.
        self.pl = np.logspace(5, 2, nlev)
        # Physically plausible flux magnitudes in W m-2: shortwave decays with
        # depth, longwave grows toward the surface, net is their difference.
        self.SW_flux_up = np.linspace(5.0, 30.0, nlev)
        self.LW_flux_up = np.linspace(240.0, 300.0, nlev)
        self.flux_up_total = self.SW_flux_up + self.LW_flux_up
        self.SW_flux_down = np.linspace(340.0, 100.0, nlev)
        self.LW_flux_down = np.linspace(0.5, 250.0, nlev)
        self.flux_down_total = self.SW_flux_down + self.LW_flux_down
        self.net_flux = self.flux_up_total - self.flux_down_total


@pytest.fixture()
def captured_ax(monkeypatch):
    """Record the axes plot_fluxes draws on, surviving its internal close."""
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


@pytest.mark.physics_invariant
def test_plot_fluxes_draws_signed_fluxes_against_pressure(tmp_path, captured_ax):
    """Upward fluxes plot as given and downward fluxes are negated.

    The sign convention of the figure is that the horizontal axis is the
    upward-directed flux, so a dropped -1 on the downward branches would
    flip half the figure; the line data pins the convention.
    """
    atm = _FluxColumn()
    out = tmp_path / 'fluxes.png'
    plot_fluxes(atm, filename=str(out))

    ax = captured_ax['ax']
    # One zero axvline plus seven flux curves.
    assert len(ax.lines) == 8
    # Line 1 is the total upward flux, drawn unmodified.
    np.testing.assert_allclose(ax.lines[1].get_xdata(), atm.flux_up_total, rtol=1e-12)
    # Line 4 is the total downward flux, drawn negated (sign convention).
    np.testing.assert_allclose(ax.lines[4].get_xdata(), -atm.flux_down_total, rtol=1e-12)
    # Downward fluxes are positive inputs, so the drawn curve must be negative.
    assert np.all(ax.lines[4].get_xdata() < 0.0)
    # The vertical coordinate is pressure in bar (Pa * 1e-5).
    np.testing.assert_allclose(ax.lines[1].get_ydata(), atm.pl * 1e-5, rtol=1e-12)


def test_plot_fluxes_axes_scales_labels_and_output(tmp_path, captured_ax):
    """The figure uses symlog/log axes, labels both axes, and writes the file.

    A pressure column is plotted with the surface at the bottom, so the
    y-limits must be inverted (larger pressure below smaller pressure).
    """
    atm = _FluxColumn()
    out = tmp_path / 'fluxes.png'
    plot_fluxes(atm, filename=str(out))

    ax = captured_ax['ax']
    assert ax.get_xscale() == 'symlog'
    assert ax.get_yscale() == 'log'
    assert 'flux' in ax.get_xlabel().lower()
    assert 'W m-2' in ax.get_xlabel()
    assert 'Pressure' in ax.get_ylabel()
    # Inverted vertical axis: top limit smaller than bottom limit.
    bottom, top = ax.get_ylim()
    assert top < bottom
    legend_texts = [t.get_text() for t in ax.get_legend().get_texts()]
    assert 'UP' in legend_texts and 'DN' in legend_texts and 'NET' in legend_texts
    assert out.stat().st_size > 0


def test_plot_fluxes_zero_flux_column_is_finite(tmp_path, captured_ax):
    """An all-zero flux column (radiative equilibrium limit) still renders.

    The symlog x-axis exists precisely so zero-crossing fluxes plot without
    a log-domain error; the degenerate all-zero column is its edge case.
    """
    atm = _FluxColumn()
    zeros = np.zeros_like(atm.pl)
    for name in (
        'SW_flux_up',
        'LW_flux_up',
        'flux_up_total',
        'SW_flux_down',
        'LW_flux_down',
        'flux_down_total',
        'net_flux',
    ):
        setattr(atm, name, zeros.copy())
    out = tmp_path / 'fluxes_zero.png'
    plot_fluxes(atm, filename=str(out))

    ax = captured_ax['ax']
    np.testing.assert_allclose(ax.lines[7].get_xdata(), zeros, atol=0.0)
    assert np.all(np.isfinite(ax.lines[7].get_xdata()))
    assert out.stat().st_size > 0
