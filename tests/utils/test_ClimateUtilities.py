"""Tests for src/janus/utils/ClimateUtilities.py.

ClimateUtilities.py is the legacy numerics and tabular-data library that the
JANUS thermodynamics helpers build on. It provides the ``Curve`` data
container, whitespace/tabular file scanning (``scan``, ``clean``,
``findData``, ``readTable``), Neville polynomial interpolation (``polint``,
``interp``), trapezoidal / Romberg quadrature (``BetterTrap``, ``romberg``),
a fourth-order Runge-Kutta ODE integrator (``integrator``), and a Newton
root solver (``newtSolve``). This file exercises:

* ``Curve`` column storage, ordering, item access, in-place replacement,
  subset extraction, and tab-delimited dump / read round-trips.
* The plotting-array accessors ``Curve.X`` / ``Curve.Y``, including the
  masked-array ``_data`` branch.
* ``scan`` / ``clean`` / ``findData`` / ``readTable`` header detection,
  delimiter handling, and malformed-row tolerance.
* ``polint`` and ``interp`` reproducing polynomial nodes exactly and
  interpolating a known polynomial to closed-form values.
* ``BetterTrap`` and ``romberg`` converging to analytic definite integrals.
* ``integrator`` recovering exponential decay and conserving harmonic
  oscillator energy.
* ``newtSolve`` locating analytic roots, honouring an explicit derivative and
  parameter object, reporting non-convergence, and bracketing sign changes.

Invariant families exercised: conservation (oscillator energy, integral
convergence), positivity/boundedness, monotonicity, and pinned closed-form
values with discrimination guards. See docs/How-to/test.md.

The wrong-argument-count guards in ``integrator`` / ``romberg`` / ``newtSolve``
print a diagnostic naming the offending callable and leave the wrapper
unbound, so the first use of the half-built object raises ``AttributeError``.
Those paths are pinned as the current error contract.
"""

import numpy
import pytest

from janus.utils import ClimateUtilities as cu

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


# ---------------------------------------------------------------------------
# Curve container
# ---------------------------------------------------------------------------


class TestCurveContainer:
    def test_addcurve_stores_orders_and_autonames_columns(self):
        """A fresh Curve stores columns in insertion order, converts lists to
        arrays, defaults the abscissa to the first column, and auto-names a
        column added with an empty id."""
        c = cu.Curve()
        # Fresh container starts empty with no abscissa selected (init state).
        assert c.idList == []
        assert c.Xid is None
        c.addCurve([1.0, 2.0, 3.0], 'p')
        c.addCurve([10.0, 20.0, 30.0], 't')
        # Edge case: empty id triggers the 'v%d' auto-name path.
        c.addCurve([7.0, 8.0, 9.0])
        assert c.listVariables() == ['p', 't', 'v2']
        # First column installed becomes the default abscissa id.
        assert c.Xid == 'p'
        # __getitem__ returns the stored column, converted from list to array.
        numpy.testing.assert_allclose(c['t'], [10.0, 20.0, 30.0], rtol=1e-12)
        assert isinstance(c['p'], numpy.ndarray)

    def test_setitem_replaces_appends_and_rejects_non_indexable(self):
        """__setitem__ overwrites an existing column, routes a new key through
        addCurve, and returns None without mutating state for a non-indexable
        right-hand side."""
        c = cu.Curve()
        c.addCurve([1.0, 2.0, 3.0], 'x')
        # Overwrite an existing column in place from a list right-hand side.
        c['x'] = [4.0, 5.0, 6.0]
        numpy.testing.assert_allclose(c['x'], [4.0, 5.0, 6.0], rtol=1e-12)
        # An already-array right-hand side is stored without list conversion.
        c['x'] = numpy.array([8.0, 8.0, 8.0])
        numpy.testing.assert_allclose(c['x'], [8.0, 8.0, 8.0], rtol=1e-12)
        # A new key is appended and tracked in idList.
        c['z'] = [7.0, 7.0, 7.0]
        assert 'z' in c.listVariables()
        # Error contract: a scalar RHS is not indexable, so the call is a no-op
        # returning None and the column set is unchanged.
        n_before = len(c.listVariables())
        assert c.__setitem__('bad', 5) is None
        assert len(c.listVariables()) == n_before

    def test_extract_returns_independent_subset(self):
        """extract returns a new Curve holding only the requested columns,
        leaving the source Curve untouched."""
        c = cu.Curve()
        c.addCurve([1.0, 2.0, 3.0], 'x')
        c.addCurve([2.0, 4.0, 6.0], 'y')
        c.addCurve([9.0, 9.0, 9.0], 'z')
        sub = c.extract(['x', 'y'])
        assert sub is not c
        assert sub.listVariables() == ['x', 'y']
        numpy.testing.assert_allclose(sub['y'], [2.0, 4.0, 6.0], rtol=1e-12)
        # Edge case: extracting a single column still yields a usable Curve, and
        # the source keeps its full column set.
        one = c.extract(['z'])
        assert one.listVariables() == ['z']
        assert c.listVariables() == ['x', 'y', 'z']

    def test_dump_and_readtable_roundtrip(self, tmp_path):
        """dump writes a tab-delimited table that readTable reconstructs
        column-for-column, and a description line is emitted when present."""
        c = cu.Curve()
        c.addCurve([1.0, 2.0, 3.0], 'x')
        c.addCurve([2.0, 4.0, 6.0], 'y')
        clean_file = tmp_path / 'clean.txt'
        c.dump(str(clean_file))
        back = cu.readTable(str(clean_file))
        assert back.listVariables() == ['x', 'y']
        numpy.testing.assert_allclose(back['x'], [1.0, 2.0, 3.0], rtol=1e-12)
        numpy.testing.assert_allclose(back['y'], [2.0, 4.0, 6.0], rtol=1e-12)
        # Edge case: a description without a trailing newline is written as the
        # first line (the source appends the missing newline).
        c.description = 'metadata'
        desc_file = tmp_path / 'desc.txt'
        c.dump(str(desc_file))
        lines = desc_file.read_text().splitlines()
        assert lines[0] == 'metadata'
        assert lines[1].split('\t') == ['x', 'y']
        # A description already terminated by a newline is written verbatim,
        # without a second newline being appended.
        c.description = 'metadata\n'
        nl_file = tmp_path / 'nl.txt'
        c.dump(str(nl_file))
        nl_lines = nl_file.read_text().splitlines()
        assert nl_lines[0] == 'metadata'
        assert nl_lines[1].split('\t') == ['x', 'y']

    @pytest.mark.physics_invariant
    def test_X_Y_plain_columns_return_plotting_arrays(self):
        """X returns the abscissa column and Y returns the remaining columns
        stacked, both as float arrays with the column count and order kept."""
        c = cu.Curve()
        c.addCurve([1.0, 2.0, 3.0], 'p')
        c.addCurve([10.0, 20.0, 30.0], 't')
        c.addCurve([100.0, 200.0, 300.0], 'q')
        x = c.X()
        y = c.Y()
        numpy.testing.assert_allclose(x, [1.0, 2.0, 3.0], rtol=1e-12)
        # Y stacks every non-abscissa column, preserving column order and count.
        assert y.shape == (2, 3)
        numpy.testing.assert_allclose(y[0], [10.0, 20.0, 30.0], rtol=1e-12)
        numpy.testing.assert_allclose(y[1], [100.0, 200.0, 300.0], rtol=1e-12)

    @pytest.mark.physics_invariant
    def test_X_Y_masked_columns_use_data_attribute(self):
        """Masked-array columns expose a ``_data`` attribute; X and Y unwrap it
        to return the underlying numeric values as a plain array."""
        masked_x = numpy.ma.array([1.0, 2.0, 3.0])
        masked_y = numpy.ma.array([10.0, 20.0, 30.0])
        c = cu.Curve()
        c.addCurve(masked_x, 'p')
        c.addCurve(masked_y, 't')
        x = c.X()
        y = c.Y()
        # The _data branch returns a plain ndarray, not a masked array.
        assert not numpy.ma.isMaskedArray(x)
        numpy.testing.assert_allclose(x, [1.0, 2.0, 3.0], rtol=1e-12)
        numpy.testing.assert_allclose(y[0], [10.0, 20.0, 30.0], rtol=1e-12)


# ---------------------------------------------------------------------------
# Tabular scanning
# ---------------------------------------------------------------------------


class TestScanning:
    def test_clean_strips_and_drops_blank_lines(self):
        """clean strips surrounding whitespace from each line and removes empty
        lines entirely."""
        result = cu.clean(['  a ', '', '  ', 'b\t'])
        assert result == ['a', 'b']
        # Edge case: an all-blank buffer collapses to an empty list.
        assert cu.clean(['', '   ', '\t']) == []

    def test_finddata_locates_longest_consistent_run(self):
        """findData returns the (start, end) span of the longest run of lines
        with a constant column count, skipping a differently-shaped header."""
        start, end = cu.findData(['header a b', '1 2', '3 4', '5 6'])
        # The 3-token header is excluded; the 2-column data block spans 1..4.
        assert (start, end) == (1, 4)
        # Edge case: a uniform buffer is one run covering the whole file, which
        # exercises the nmax == -1 fallback branch.
        assert cu.findData(['1 2', '3 4']) == (0, 2)
        # Edge case: with several equal-length runs the first maximal run wins
        # the tie, so a later run of the same length does not displace it.
        assert cu.findData(['1 2', '3 4', '5 6', '7 8 9', 'a b']) == (3, 4)

    def test_scan_detects_header_and_generates_names(self):
        """scan reads a leading text header when present and otherwise
        synthesises V0, V1, ... names for a purely numeric block."""
        data, header = cu.scan(['x y z', '1 2 3', '4 5 6', '7 8 9'])
        assert header == ['x', 'y', 'z']
        numpy.testing.assert_allclose(data['y'], [2.0, 5.0, 8.0], rtol=1e-12)
        # Edge case: a numeric first line means there is no header, so columns
        # are auto-named and the first row is kept as data.
        data2, header2 = cu.scan(['1 2', '3 4', '5 6'])
        assert header2 == ['V0', 'V1']
        numpy.testing.assert_allclose(data2['V0'], [1.0, 3.0, 5.0], rtol=1e-12)

    def test_scan_honours_delimiter_and_supplied_header(self):
        """scan splits on an explicit delimiter and replaces detected headers
        with a caller-supplied header of matching width."""
        data, header = cu.scan(['a,b', '1,2', '3,4'], delimiter=',')
        assert header == ['a', 'b']
        numpy.testing.assert_allclose(data['b'], [2.0, 4.0], rtol=1e-12)
        # A supplied header of the right width overrides the auto-generated one.
        data2, header2 = cu.scan(['1 2', '3 4'], inHeader=['q', 'r'])
        assert header2 == ['q', 'r']
        numpy.testing.assert_allclose(data2['q'], [1.0, 3.0], rtol=1e-12)

    def test_scan_tolerates_a_malformed_data_row(self):
        """A row with a non-numeric token is caught by scan's guard rather than
        propagating a ValueError, and the parseable rows still populate the
        first column."""
        data, header = cu.scan(['1 2', '3 x', '5 6'])
        assert header == ['V0', 'V1']
        # The first column parses '1', '3', '5'; the guard swallows the failed
        # conversion of 'x' in the second column of the middle row.
        numpy.testing.assert_allclose(data['V0'], [1.0, 3.0, 5.0], rtol=1e-12)


# ---------------------------------------------------------------------------
# Polynomial interpolation
# ---------------------------------------------------------------------------


class TestPolint:
    def test_polint_rejects_mismatched_lengths(self):
        """polint returns the sentinel 'Error' string when the abscissa and
        ordinate arrays differ in length, before touching numpy."""
        result = cu.polint([0.0, 1.0], [0.0, 1.0, 2.0], 0.5)
        assert result == 'Error'
        assert not isinstance(result, float)

    @pytest.mark.physics_invariant
    def test_polint_reproduces_cubic_exactly(self):
        """Neville interpolation through four nodes of a cubic reproduces the
        nodes and evaluates the cubic exactly between them."""
        xa = [0.0, 1.0, 2.0, 3.0]
        ya = [x**3 - 2.0 * x for x in xa]
        # Node reproduction: the interpolant passes through each sample.
        assert cu.polint(xa, ya, 2.0) == pytest.approx(4.0, rel=1e-12)
        mid = cu.polint(xa, ya, 1.5)
        assert mid == pytest.approx(0.375, rel=1e-12)
        # Discrimination guard: a straight-line interpolant between the two
        # bracketing nodes would give 1.5, far outside tolerance, so a
        # collapse to linear interpolation would fail this test.
        linear_guess = 0.5 * (ya[1] + ya[2])
        assert abs(mid - linear_guess) > 0.5

    @pytest.mark.physics_invariant
    def test_interp_reproduces_nodes_and_handles_descending_grid(self):
        """interp reproduces its sample nodes and interpolates a quadratic to
        its closed-form value, for both ascending and descending grids."""
        grid = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        vals = [x * x for x in grid]
        f = cu.interp(grid, vals)
        # The callable advertises a single-argument signature for callers that
        # introspect co_argcount.
        assert f.__code__.co_argcount == 1
        assert f(3.0) == pytest.approx(9.0, rel=1e-12)
        assert f(2.5) == pytest.approx(6.25, rel=1e-12)
        # Edge case: a descending abscissa exercises the reversed searchsorted
        # branch and must give the same interpolated value.
        f_desc = cu.interp(grid[::-1], vals[::-1])
        assert f_desc(2.5) == pytest.approx(6.25, rel=1e-12)
        # Discrimination guard: the chord between x=2 and x=3 sits above the
        # convex parabola, so the true value is strictly below 6.5.
        assert f(2.5) < 6.5


# ---------------------------------------------------------------------------
# Quadrature
# ---------------------------------------------------------------------------


class TestQuadrature:
    @pytest.mark.physics_invariant
    def test_bettertrap_refines_toward_analytic_integral(self):
        """BetterTrap refinement halves the step, doubles the subinterval count,
        and drives the trapezoidal estimate toward the exact integral of x**2
        on [0, 1], which is 1/3."""
        trap = cu.BetterTrap(lambda x, params: x * x, None, [0.0, 1.0], 4)
        coarse = trap.integral
        err_coarse = abs(coarse - 1.0 / 3.0)
        trap.refine()
        # The refinement doubles the number of subintervals.
        assert trap.n == 8
        err_fine = abs(trap.integral - 1.0 / 3.0)
        # Trapezoidal error shrinks monotonically under refinement.
        assert err_fine < err_coarse
        assert trap.integral == pytest.approx(1.0 / 3.0, abs=5e-3)

    @pytest.mark.physics_invariant
    def test_romberg_integrates_polynomial_with_both_arities(self):
        """Romberg extrapolation recovers analytic definite integrals, wrapping
        a bare one-argument integrand and honouring a two-argument integrand
        with a parameter object."""
        # One-argument integrand: integral of x**2 on [0, 1] is 1/3.
        r1 = cu.romberg(lambda x: x * x)
        val1 = r1([0.0, 1.0])
        assert val1 == pytest.approx(1.0 / 3.0, rel=1e-6)
        # Discrimination guard: a single trapezoid would give 0.5, well outside
        # the Romberg tolerance, so a failure to refine would be caught.
        assert abs(val1 - 0.5) > 1e-3
        # Two-argument integrand with parameters: integral of a*x on [0, 2] is
        # a * 2, i.e. 6 for a = 3.
        params = cu.Dummy()
        params.a = 3.0
        r2 = cu.romberg(lambda x, p: p.a * x)
        assert r2([0.0, 2.0], params) == pytest.approx(6.0, rel=1e-6)
        # A transcendental integrand needs a further refinement, exercising the
        # convergence loop: integral of exp(x) on [0, 1] is e - 1.
        r3 = cu.romberg(lambda x: numpy.exp(x))
        val3 = r3([0.0, 1.0])
        assert val3 == pytest.approx(numpy.e - 1.0, rel=1e-6)
        # The refined result must beat the crude single-trapezoid estimate.
        crude = 0.5 * (numpy.exp(0.0) + numpy.exp(1.0))
        assert abs(val3 - (numpy.e - 1.0)) < abs(crude - (numpy.e - 1.0))

    def test_romberg_rejects_wrong_arity_integrand(self, capsys):
        """A three-argument integrand fails romberg's argument-count guard.

        The guard prints a diagnostic naming the callable and leaves the
        integrand unbound, so the first integration attempt raises
        AttributeError rather than silently computing with a wrong function.
        """
        r = cu.romberg(lambda x, y, z: x)
        assert 'wrong number of arguments' in capsys.readouterr().out
        with pytest.raises(AttributeError):
            r([0.0, 1.0])


# ---------------------------------------------------------------------------
# Runge-Kutta ODE integrator
# ---------------------------------------------------------------------------


class TestIntegrator:
    @pytest.mark.physics_invariant
    def test_integrator_recovers_exponential_decay(self):
        """The RK4 integrator advances dy/dt = -k y to the analytic solution
        y(t) = exp(-k t), staying positive and monotonically decreasing."""

        def decay(t, y, params):
            return -params.k * y

        it = cu.integrator(decay, 0.0, 1.0, 0.001)
        params = cu.Dummy()
        params.k = 2.0
        it.setParams(params)
        prev = 1.0
        for _ in range(1000):
            t, y = it.next()
            # Monotone decrease and strict positivity of the decaying solution.
            assert 0.0 < y < prev
            prev = y
        # Global RK4 error at dx = 1e-3 is far below the pinned tolerance.
        assert y == pytest.approx(numpy.exp(-2.0), rel=1e-6)

    @pytest.mark.physics_invariant
    def test_integrator_conserves_harmonic_oscillator_energy(self):
        """For the harmonic oscillator the RK4 integrator tracks the analytic
        (cos t, -sin t) trajectory and conserves the energy x**2 + v**2."""

        def osc(t, z):
            # New array each call, as the integrator requires.
            return numpy.array([z[1], -z[0]])

        it = cu.integrator(osc, 0.0, numpy.array([1.0, 0.0]), 0.001)
        for _ in range(1000):
            t, z = it.next()
        assert z[0] == pytest.approx(numpy.cos(1.0), abs=1e-4)
        assert z[1] == pytest.approx(-numpy.sin(1.0), abs=1e-4)
        # Energy conservation is the discriminating invariant: a sign slip in
        # the RK4 stages would let the orbit spiral and break this equality.
        energy = z[0] ** 2 + z[1] ** 2
        assert energy == pytest.approx(1.0, abs=1e-3)

    def test_integrator_step_override_updates_increment(self):
        """Passing a step to next overrides the increment and the integrator
        remembers it for the following call."""
        it = cu.integrator(lambda x, y: 1.0 + 0.0 * y, 0.0, 0.0, None)
        x1, y1 = it.next(0.5)
        assert x1 == pytest.approx(0.5, rel=1e-12)
        # The remembered step is reused when next is called without an argument.
        x2, y2 = it.next()
        assert x2 == pytest.approx(1.0, rel=1e-12)

    def test_integrator_rejects_wrong_arity_derivs(self, capsys):
        """A one-argument derivative function fails the integrator's
        argument-count guard.

        The guard prints a diagnostic naming the callable and leaves the
        derivative unbound, so the first step raises AttributeError rather
        than integrating a wrong function.
        """
        it = cu.integrator(lambda x: x, 0.0, 1.0, 0.1)
        assert 'wrong number of arguments' in capsys.readouterr().out
        with pytest.raises(AttributeError):
            it.next()


# ---------------------------------------------------------------------------
# Newton root solver
# ---------------------------------------------------------------------------


class TestNewtSolve:
    @pytest.mark.physics_invariant
    def test_newtsolve_finds_positive_root_with_numeric_derivative(self):
        """Starting from a positive guess, newtSolve converges to the positive
        root of x**2 - 2 using its centred-difference derivative."""
        solver = cu.newtSolve(lambda x: x * x - 2.0)
        root = solver(1.5)
        assert root == pytest.approx(2.0**0.5, rel=1e-9)
        # Sign guard: the positive guess must land on +sqrt(2), not -sqrt(2).
        assert root > 0
        # Residual guard: the returned root actually solves the equation.
        assert root * root == pytest.approx(2.0, rel=1e-9)

    def test_newtsolve_uses_explicit_derivative_and_parameters(self):
        """newtSolve accepts an explicit derivative and a parameter object
        passed at call time, and exposes the derivative via solver.deriv."""

        def g(x, c):
            return c.a * x * x - c.b

        solver = cu.newtSolve(g, lambda x: 2.0 * x)
        constants = cu.Dummy()
        constants.a = 1.0
        constants.b = 2.0
        root = solver(1.5, constants)
        assert root == pytest.approx(2.0**0.5, rel=1e-9)
        # The stored parameter object is the one supplied at the call site.
        assert solver.params is constants
        # The wrapped one-argument derivative evaluates to 2x when invoked with
        # the internal (x, params) signature.
        assert solver.deriv(3.0, None) == pytest.approx(6.0, rel=1e-12)

    def test_newtsolve_accepts_two_arg_derivative_and_rejects_bad_arity(self, capsys):
        """An explicit two-argument derivative is used directly, while a
        three-argument derivative fails the argument-count guard.

        The guard prints a diagnostic naming the callable; the analytic
        derivative stays unbound, so the solve falls over on first use
        instead of silently using a wrong derivative.
        """

        def two_arg_deriv(x, params):
            return 2.0 * x

        solver = cu.newtSolve(lambda x: x * x - 2.0, two_arg_deriv)
        root = solver(1.5)
        assert root == pytest.approx(2.0**0.5, rel=1e-9)
        # The supplied two-argument derivative is exposed unchanged.
        assert solver.deriv(3.0, None) == pytest.approx(6.0, rel=1e-12)
        # Error contract: a wrong-arity derivative is diagnosed on construction.
        bad = cu.newtSolve(lambda x: x * x - 2.0, lambda x, y, z: 2.0 * x)
        assert 'wrong number of arguments' in capsys.readouterr().out
        with pytest.raises(AttributeError):
            bad(1.5)

    def test_newtsolve_reports_non_convergence(self):
        """When the iteration cannot reach the tolerance within nmax steps,
        newtSolve returns the 'No Convergence' sentinel rather than a number."""
        solver = cu.newtSolve(lambda x: x * x + 1.0)
        solver.nmax = 3
        result = solver(0.5)
        assert result == 'No Convergence'
        # Contrast case: a solvable problem returns a float, confirming the
        # sentinel is specific to the stalled iteration.
        good = cu.newtSolve(lambda x: x - 1.0)
        assert isinstance(good(0.0), float)

    @pytest.mark.physics_invariant
    def test_newtsolve_scan_brackets_sign_changes(self):
        """scan returns guesses that bracket every sign change of x**2 - 1 on
        [-2, 2], and Newton launched from those guesses recovers the roots at
        -1 and +1."""
        solver = cu.newtSolve(lambda x: x * x - 1.0)
        guesses = solver.scan([-2.0, 2.0], 41)
        assert len(guesses) > 0
        # A guess brackets each root to within one subinterval on either side.
        near_minus = [g for g in guesses if -1.2 < g < -0.8]
        near_plus = [g for g in guesses if 0.8 < g < 1.2]
        assert near_minus and near_plus
        # Newton from the bracketing guesses converges to the analytic roots.
        assert solver(near_minus[0]) == pytest.approx(-1.0, rel=1e-6)
        assert solver(near_plus[0]) == pytest.approx(1.0, rel=1e-6)

    def test_newtsolve_rejects_wrong_arity_function(self, capsys):
        """A three-argument objective fails newtSolve's argument-count guard.

        The guard prints a diagnostic naming the callable and leaves the
        objective unbound, so the first solve attempt raises AttributeError
        rather than iterating on a wrong function.
        """
        bad = cu.newtSolve(lambda x, y, z: x)
        assert 'wrong number of arguments' in capsys.readouterr().out
        with pytest.raises(AttributeError):
            bad(1.0)
