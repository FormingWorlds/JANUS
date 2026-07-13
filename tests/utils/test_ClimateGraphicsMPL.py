"""Tests for src/janus/utils/ClimateGraphicsMPL.py.

Exercises the legacy Curve-plotting driver: multi-curve line plots with
labels, scatter markers, axis switching, log scaling, and axis reversal;
the contour plotter with and without explicit coordinate axes; and the
small resource/plotObj helper classes. See docs/How-to/test.md.
"""

import matplotlib

matplotlib.use('Agg')

import numpy as np  # noqa: E402
import pytest  # noqa: E402

from janus.utils import ClimateGraphicsMPL as cgm  # noqa: E402
from janus.utils.ClimateUtilities import Curve  # noqa: E402

# Importing the module flips the global backend request to TkAgg and turns on
# interactive mode; force the headless backend back on before any figure is
# created so the whole test session stays window-free.
matplotlib.use('Agg', force=True)

import matplotlib.pyplot as plt  # noqa: E402

plt.ioff()

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.fixture(autouse=True)
def _close_figures():
    """Close every figure the driver leaves open."""
    yield
    plt.close('all')


def _curve_two_columns():
    """A Curve with an X column and two Y columns of distinct amplitude."""
    c = Curve()
    x = np.linspace(0.0, 1.0, 5)
    c.addCurve(x, 'x')
    c.addCurve(2.0 * x, 'ya', label='double')
    c.addCurve(x**2, 'yb')  # no label: the legend must fall back to the id
    c.Xlabel = 'time [s]'
    c.Ylabel = 'signal [m]'
    c.PlotTitle = 'ramp'
    return c


def test_plot_draws_all_columns_with_labels_and_legend():
    """One line per Y column, data unchanged, legend from label or id.

    The two Y columns have distinct amplitudes so a column mix-up moves the
    drawn data far outside tolerance; the unlabeled column must appear in
    the legend under its id.
    """
    c = _curve_two_columns()
    obj = cgm.plot(c)

    ax = obj.plot.axes[0]
    assert len(ax.get_lines()) == 2
    np.testing.assert_allclose(ax.get_lines()[0].get_ydata(), 2.0 * c.X(), rtol=1e-12)
    np.testing.assert_allclose(ax.get_lines()[1].get_ydata(), c.X() ** 2, rtol=1e-12)
    assert ax.get_xlabel() == 'time [s]'
    assert ax.get_ylabel() == 'signal [m]'
    assert ax.get_title() == 'ramp'
    legend_texts = [t.get_text() for t in ax.get_legend().get_texts()]
    assert 'double' in legend_texts
    assert 'yb' in legend_texts


def test_plot_switchxy_swaps_axes_and_data():
    """switchXY exchanges the plotted coordinates and the axis labels.

    With switchXY the Y column supplies the horizontal coordinate, so the
    first line's x-data must equal the column values and the axis labels
    must swap; asserting both discriminates a label-only swap from a real
    coordinate swap.
    """
    c = _curve_two_columns()
    c.switchXY = 1
    obj = cgm.plot(c)

    ax = obj.plot.axes[0]
    np.testing.assert_allclose(ax.get_lines()[0].get_xdata(), 2.0 * c.X(), rtol=1e-12)
    np.testing.assert_allclose(ax.get_lines()[0].get_ydata(), c.X(), rtol=1e-12)
    assert ax.get_xlabel() == 'signal [m]'
    assert ax.get_ylabel() == 'time [s]'


def test_plot_log_axes_reversal_and_scatter_branch():
    """Log axes apply to both directions, reversal inverts limits, scatter
    columns draw markers without a connecting line.

    Reversal happens after the data are drawn, so the x-limits must come out
    decreasing; the scatter branch must set a marker and drop the line style.
    """
    c = Curve()
    x = np.logspace(0, 2, 4)
    c.addCurve(x, 'x')
    c.addCurve(10.0 * x, 'y', label='prop')
    c.scatter['y'] = True
    c.XlogAxis = 1
    c.YlogAxis = 1
    c.reverseX = 1
    c.reverseY = 1
    obj = cgm.plot(c)

    ax = obj.plot.axes[0]
    assert ax.get_xscale() == 'log'
    assert ax.get_yscale() == 'log'
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    assert x0 > x1
    assert y0 > y1
    line = ax.get_lines()[0]
    assert line.get_marker() != 'None'
    assert line.get_linestyle() in ('None', '')


def test_plot_multiple_curves_of_different_resolution():
    """Two Curve objects of different lengths plot onto one set of axes.

    The driver's stated purpose is overlaying data of varying resolution;
    each input column must keep its own length and values.
    """
    c1 = Curve()
    x1 = np.linspace(0.0, 1.0, 4)
    c1.addCurve(x1, 'x')
    c1.addCurve(np.sin(x1), 'y1', label='coarse')

    c2 = Curve()
    x2 = np.linspace(0.0, 1.0, 9)
    c2.addCurve(x2, 'x')
    c2.addCurve(np.cos(x2), 'y2', label='fine')

    obj = cgm.plot(c1, c2)
    ax = obj.plot.axes[0]
    lines = ax.get_lines()
    assert len(lines) == 2
    assert lines[0].get_xdata().size == 4
    assert lines[1].get_xdata().size == 9
    np.testing.assert_allclose(lines[1].get_ydata(), np.cos(x2), rtol=1e-12)


def test_contour_with_and_without_axes_and_helper_classes(capsys):
    """contour() accepts explicit or implicit coordinates; helpers behave.

    Without keyword axes the array indices are the coordinates; with x/y
    supplied the plot spans those ranges and adds a colorbar axes. The
    resource defaults are all False (nothing reversed, nothing logarithmic)
    and the plotObj save/delete stubs only print guidance.
    """
    a = np.outer(np.linspace(0.0, 1.0, 3), np.linspace(1.0, 2.0, 4))

    obj_plain = cgm.contour(a)
    assert len(obj_plain.plot.axes) == 2  # contour axes plus colorbar

    xs = np.linspace(10.0, 40.0, 4)
    ys = np.linspace(-5.0, 5.0, 3)
    obj_xy = cgm.contour(a, x=xs, y=ys)
    ax = obj_xy.plot.axes[0]
    x0, x1 = ax.get_xlim()
    assert x0 == pytest.approx(10.0, abs=1e-9)
    assert x1 == pytest.approx(40.0, abs=1e-9)

    r = cgm.resource()
    assert r.trYReverse is False
    assert r.trXLog is False and r.trYLog is False

    d = cgm.Dummy()
    d.anything = 3
    assert d.anything == 3

    p = cgm.plotObj(None, obj_plain.plot)
    p.delete()
    p.save('unused')
    printed = capsys.readouterr().out
    assert 'MatPlotLib' in printed
    assert p.plot is obj_plain.plot
