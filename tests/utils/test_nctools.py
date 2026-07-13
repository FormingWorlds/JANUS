"""Tests for src/janus/utils/nctools.py.

nctools.py writes SOCRATES-compatible netCDF input files: surface-albedo
weight files (``ncout_surf``, ``ncout_spectral_surf``), 2D and 3D
pressure-level field files (``ncout2d``, ``ncout3d``), prescribed
optical-property files (``ncout_opt_prop``), the ``.view`` writer
(``ncout_view``), and the low-level ``write_dim`` / ``write_var`` /
``create_cdf`` helpers. This file exercises:

* Round-trip integrity: values, dimension sizes, and unit/title attributes
  survive a write-then-read cycle.
* The array-size dispatch in each writer: scalar broadcast, exact-match, and
  the size-mismatch RuntimeError contract.
* The optional-attribute guards in ``write_var`` and the name-from-filename
  fallback in ``ncout2d`` / ``ncout3d``.
* The pressure-level sort that reorders an unsorted level axis.

These are pure netCDF I/O utilities (not one of the six physics sources), so
the assertions pin file structure and value round-trips rather than a physical
invariant. See docs/How-to/test.md.
"""

import numpy as np
import pytest
from netCDF4 import Dataset

import janus.utils.nctools as nctools

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


# ---------------------------------------------------------------------------
# Low-level helpers: write_dim / write_var / create_cdf
# ---------------------------------------------------------------------------


def test_write_var_sets_units_only_for_string_attributes(tmp_path):
    """write_var attaches units/title only when they are passed as strings.

    write_var is the low-level variable writer. Given string units and title it
    must store them as netCDF attributes; given None (the SOCRATES callers pass
    None for unlabelled fields) it must skip the attribute rather than store the
    literal string "None". The None case is the edge that exercises both
    attribute guards. Both variables must round-trip their values.
    """
    fpath = tmp_path / 'vars.nc'
    ds = nctools.create_cdf(str(fpath))
    ds.createDimension('x', 3)
    labelled = np.array([1.0, 2.0, 3.0], dtype='f4')
    plain = np.array([4.0, 5.0, 6.0], dtype='f4')
    nctools.write_var(ds, labelled, 'labelled', 'f4', ('x',), 'K', 'TEMPERATURE')
    nctools.write_var(ds, plain, 'plain', 'f4', ('x',), None, None)
    ds.close()

    with Dataset(str(fpath)) as check:
        np.testing.assert_allclose(check['labelled'][:], labelled, rtol=1e-6)
        np.testing.assert_allclose(check['plain'][:], plain, rtol=1e-6)
        assert check['labelled'].units == 'K'
        assert check['labelled'].title == 'TEMPERATURE'
        # The None units/title must not have been stored as attributes.
        assert 'units' not in check['plain'].ncattrs()
        assert 'title' not in check['plain'].ncattrs()


def test_write_dim_round_trips_axis_values_and_metadata(tmp_path):
    """write_dim creates a coordinate dimension, its variable, and metadata.

    write_dim underlies every axis in the SOCRATES netCDF files. A multi-point
    longitude axis must produce a dimension of the right length with its values
    and units/title attributes stored, and a degenerate single-point latitude
    axis (the edge case) must be handled as a length-1 dimension rather than
    collapsing to a scalar.
    """
    fpath = tmp_path / 'dim.nc'
    ds = nctools.create_cdf(str(fpath))
    lon = np.array([0.0, 90.0, 180.0, 270.0], dtype='f4')
    nctools.write_dim(ds, lon.size, lon, 'lon', 'f4', 'lon', 'degree', 'LONGITUDE')
    lat = np.array([45.0], dtype='f4')  # single-point axis
    nctools.write_dim(ds, lat.size, lat, 'lat', 'f4', 'lat', 'degree', 'LATITUDE')
    ds.close()

    with Dataset(str(fpath)) as check:
        assert check.dimensions['lon'].size == 4
        assert check.dimensions['lat'].size == 1
        np.testing.assert_allclose(check['lon'][:], lon, rtol=1e-6)
        assert check['lon'].units == 'degree'
        assert check['lon'].title == 'LONGITUDE'


# ---------------------------------------------------------------------------
# ncout_surf
# ---------------------------------------------------------------------------


def test_ncout_surf_scalar_and_full_array(tmp_path):
    """ncout_surf broadcasts a scalar albedo and stores a full array as-is.

    A scalar surface albedo must be broadcast to every (basis, lat, lon) cell;
    a full (basis, lat, lon) array must round-trip unchanged. The scalar path
    uses a single-latitude column (an edge case for the broadcast) and the
    array path uses two basis functions with distinct values, so a transposed
    or mis-shaped write would be caught.
    """
    lon = np.array([0.0, 120.0, 240.0])
    fscalar = tmp_path / 'uniform.surf'
    nctools.ncout_surf(str(fscalar), lon, np.array([0.0]), np.array([1]), 0.3)
    with Dataset(str(fscalar)) as ds:
        alb = np.asarray(ds['alb'][:])
        assert alb.shape == (1, 1, 3)
        np.testing.assert_allclose(alb, 0.3, rtol=1e-6)

    farr = tmp_path / 'varied.surf'
    basis = np.array([1, 2])
    lat = np.array([-30.0, 30.0])
    albvals = np.arange(basis.size * lat.size * lon.size, dtype=float).reshape(
        basis.size, lat.size, lon.size
    )
    nctools.ncout_surf(str(farr), lon, lat, basis, albvals)
    with Dataset(str(farr)) as ds:
        stored = np.asarray(ds['alb'][:])
        assert stored.shape == (2, 2, 3)
        np.testing.assert_allclose(stored, albvals, rtol=1e-6)


def test_ncout_surf_rejects_mismatched_albedo(tmp_path):
    """ncout_surf raises when the albedo count matches no expected size.

    The albedo array must be either a scalar or exactly (basis * lat * lon)
    long. An array of an intermediate length (5 here) satisfies neither branch
    and must raise rather than writing a malformed file. This is the size
    dispatch's error contract.
    """
    fpath = tmp_path / 'bad.surf'
    with pytest.raises(RuntimeError, match='dont match') as exc:
        nctools.ncout_surf(
            str(fpath),
            np.array([0.0, 90.0, 180.0]),
            np.array([0.0]),
            np.array([1]),
            np.arange(5, dtype=float),
        )
    # The error reports both the offending count (5) and the expected size (3).
    assert 5 in exc.value.args
    assert 3 in exc.value.args


# ---------------------------------------------------------------------------
# ncout_spectral_surf
# ---------------------------------------------------------------------------


def test_ncout_spectral_surf_scalar_and_per_band(tmp_path):
    """ncout_spectral_surf writes a per-band surface-albedo weight file.

    A scalar albedo is broadcast across the single-band, single-latitude grid
    (the degenerate edge case), and a full (bands, basis, lat, lon) array is
    stored as-is with three bands. Both must round-trip, and the band dimension
    must carry the requested number of bands so a dropped band axis is caught.
    """
    lon = np.array([0.0, 90.0, 180.0])
    fscalar = tmp_path / 'uniform_spec.surf'
    nctools.ncout_spectral_surf(str(fscalar), lon, np.array([0.0]), 1, 0.2)
    with Dataset(str(fscalar)) as ds:
        alb = np.asarray(ds['alb'][:])
        assert alb.shape == (1, 1, 1, 3)
        np.testing.assert_allclose(alb, 0.2, rtol=1e-6)
        assert ds.dimensions['bands'].size == 1

    farr = tmp_path / 'perband_spec.surf'
    lat = np.array([0.0])
    lon2 = np.array([0.0, 90.0])
    bands = 3
    albvals = np.arange(bands * 1 * lat.size * lon2.size, dtype=float).reshape(
        bands, 1, lat.size, lon2.size
    )
    nctools.ncout_spectral_surf(str(farr), lon2, lat, bands, albvals)
    with Dataset(str(farr)) as ds:
        stored = np.asarray(ds['alb'][:])
        assert stored.shape == (3, 1, 1, 2)
        np.testing.assert_allclose(stored, albvals, rtol=1e-6)


def test_ncout_spectral_surf_rejects_mismatched_albedo(tmp_path):
    """ncout_spectral_surf raises on an albedo count matching no expected size.

    The albedo must be scalar or exactly (lon * lat * bands) long. A length-5
    array against a two-band, single-lat, two-lon grid (expected 4) matches
    neither branch and must raise. This is the error contract for the spectral
    surface writer.
    """
    fpath = tmp_path / 'bad_spec.surf'
    with pytest.raises(RuntimeError, match='dont match') as exc:
        nctools.ncout_spectral_surf(
            str(fpath),
            np.array([0.0, 90.0]),
            np.array([0.0]),
            2,
            np.arange(5, dtype=float),
        )
    # The error reports the offending count (5) and the expected size (4 = 2*1*2).
    assert 5 in exc.value.args
    assert 4 in exc.value.args


# ---------------------------------------------------------------------------
# ncout2d
# ---------------------------------------------------------------------------


def test_ncout2d_scalar_broadcast_and_named_field(tmp_path):
    """ncout2d writes a single-level field with an explicit variable name.

    A scalar value is broadcast across the (lat, lon) grid, and an exact-length
    array is stored in order. The single-latitude grid is the edge case for the
    broadcast. Both are written under the supplied variable name and must
    round-trip with the requested units attribute.
    """
    lon = np.array([0.0, 90.0, 180.0])
    lat = np.array([0.0])
    fscalar = tmp_path / 'uniform2d.nc'
    nctools.ncout2d(
        str(fscalar), lon, lat, 1365.0, name='sol', longname='Solar Irradiance', units='W m-2'
    )
    with Dataset(str(fscalar)) as ds:
        vals = np.asarray(ds['sol'][:])
        assert vals.shape == (1, 3)
        np.testing.assert_allclose(vals, 1365.0, rtol=1e-6)
        assert ds['sol'].units == 'W m-2'

    farr = tmp_path / 'varied2d.nc'
    field = np.array([10.0, 20.0, 30.0])
    nctools.ncout2d(str(farr), lon, lat, field, name='temp', units='K')
    with Dataset(str(farr)) as ds:
        np.testing.assert_allclose(np.asarray(ds['temp'][:]).ravel(), field, rtol=1e-6)
        assert ds['temp'].units == 'K'


def test_ncout2d_derives_name_from_filename(tmp_path):
    """ncout2d derives the variable name from the file extension when unnamed.

    When ``name`` is not a string, ncout2d takes the field name from the text
    after the first dot in the file path (VULCAN/SOCRATES convention). A file
    ending in ``.stoa`` must yield a variable named ``stoa``. Passing an
    exact-length field also exercises the (lat, lon) match branch.
    """
    lon = np.array([0.0, 90.0, 180.0])
    lat = np.array([0.0])
    fpath = tmp_path / 'field.stoa'
    nctools.ncout2d(str(fpath), lon, lat, np.array([1.0, 2.0, 3.0]))
    with Dataset(str(fpath)) as ds:
        assert 'stoa' in ds.variables
        derived = [v for v in ds.variables if v not in ('lon', 'lat')]
        assert derived == ['stoa']
        np.testing.assert_allclose(
            np.asarray(ds['stoa'][:]).ravel(), np.array([1.0, 2.0, 3.0]), rtol=1e-6
        )


def test_ncout2d_rejects_mismatched_field(tmp_path):
    """ncout2d raises when the field length matches neither scalar nor lat*lon.

    A field that is neither a scalar nor exactly (lat * lon) long is ambiguous
    and must raise rather than reshaping arbitrarily. This is the size
    dispatch's error contract.
    """
    fpath = tmp_path / 'bad2d.nc'
    with pytest.raises(RuntimeError, match='ncout2d') as exc:
        nctools.ncout2d(
            str(fpath),
            np.array([0.0, 90.0, 180.0]),
            np.array([0.0]),
            np.arange(5, dtype=float),
            name='x',
        )
    # The error reports the offending count (5) and the expected size (3 = 3*1).
    assert 5 in exc.value.args
    assert 3 in exc.value.args


# ---------------------------------------------------------------------------
# ncout3d
# ---------------------------------------------------------------------------


def test_ncout3d_scalar_field_and_single_level(tmp_path):
    """ncout3d broadcasts a scalar over all levels and handles a single level.

    A scalar value must fill every (plev, lat, lon) cell, and an unsorted
    pressure axis must be reordered so the stored levels ascend. A single-level
    column (the edge case) skips the sort entirely and stores exactly one
    level. Both must round-trip the constant field.
    """
    lon = np.array([0.0, 120.0])
    lat = np.array([0.0])
    fmulti = tmp_path / 'scalar3d.nc'
    p_unsorted = np.array([1.0e5, 1.0e4, 1.0e3])
    nctools.ncout3d(
        str(fmulti), lon, lat, p_unsorted, 250.0, name='t', longname='Temperature', units='K'
    )
    with Dataset(str(fmulti)) as ds:
        vals = np.asarray(ds['t'][:])
        assert vals.shape == (3, 1, 2)
        np.testing.assert_allclose(vals, 250.0, rtol=1e-6)
        plev = np.asarray(ds['plev'][:])
        assert np.all(np.diff(plev) > 0.0)  # reordered to ascending pressure
        assert ds['t'].units == 'K'

    fsingle = tmp_path / 'single3d.nc'
    nctools.ncout3d(str(fsingle), lon, lat, np.array([5.0e4]), 300.0, name='t', units='K')
    with Dataset(str(fsingle)) as ds:
        assert ds.dimensions['plev'].size == 1
        np.testing.assert_allclose(np.asarray(ds['t'][:]), 300.0, rtol=1e-6)


def test_ncout3d_full_field_reorders_and_derives_defaults(tmp_path):
    """ncout3d reshapes a full field, reorders unsorted levels, and defaults.

    A full (level * lat * lon) field must be reshaped onto the grid, and an
    unsorted pressure axis reordered so each level's value follows its
    pressure. With ``name``/``longname``/``units`` all None, ncout3d derives the
    name from the filename and stores no units string. The per-level values
    (10, 30, 20 at pressures 1e3, 1e5, 1e4) must end up as (10, 20, 30) once the
    levels ascend, which discriminates a reorder that moves pressures but not
    data.
    """
    lon = np.array([0.0, 120.0])
    lat = np.array([0.0])
    p_unsorted = np.array([1.0e3, 1.0e5, 1.0e4])
    level_vals = np.array([10.0, 30.0, 20.0])  # one value per level, in p order
    field = np.repeat(level_vals, lat.size * lon.size)  # full level*lat*lon field
    fpath = tmp_path / 'temperature.tfield'
    nctools.ncout3d(
        str(fpath), lon, lat, p_unsorted, field, name=None, longname=None, units=None
    )
    with Dataset(str(fpath)) as ds:
        assert 'tfield' in ds.variables  # name derived from extension
        plev = np.asarray(ds['plev'][:])
        np.testing.assert_allclose(plev, np.array([1.0e3, 1.0e4, 1.0e5]), rtol=1e-6)
        # After the sort the column at (lat0, lon0) is monotone with pressure.
        np.testing.assert_allclose(
            np.asarray(ds['tfield'][:, 0, 0]), np.array([10.0, 20.0, 30.0]), rtol=1e-6
        )
        # A None units becomes the literal string 'None'; a None longname stores
        # no title attribute at all.
        assert ds['tfield'].units == 'None'
        assert 'title' not in ds['tfield'].ncattrs()


def test_ncout3d_rejects_mismatched_field(tmp_path):
    """ncout3d raises when the field size matches no expected layout.

    A field that is neither scalar, nor level-length, nor (lon * lat * level)
    long cannot be placed on the grid and must raise. This is the error
    contract of the 3D size dispatch.
    """
    fpath = tmp_path / 'bad3d.nc'
    with pytest.raises(RuntimeError, match='ncout3d') as exc:
        nctools.ncout3d(
            str(fpath),
            np.array([0.0, 120.0]),
            np.array([0.0]),
            np.array([1.0e5, 1.0e4, 1.0e3]),
            np.arange(7, dtype=float),
            name='x',
        )
    # The error reports the offending count (7) and the expected size (6 = 2*1*3).
    assert 7 in exc.value.args
    assert 6 in exc.value.args


def test_ncout3d_per_level_column_tiled_across_grid(tmp_path):
    """ncout3d tiles a per-level column across the (lat, lon) grid.

    When the field length equals the number of levels, each level's single
    value is written to every (lat, lon) cell of that level. With an
    already-sorted pressure axis (the edge case that skips the reorder) the
    tiled column at any (lat, lon) must equal the input per-level array, and a
    second horizontal point must hold the identical column.
    """
    lon = np.array([0.0, 120.0])
    lat = np.array([0.0])
    p = np.array([1.0e3, 1.0e4, 1.0e5])  # ascending: no reorder
    column = np.array([220.0, 250.0, 300.0])  # one value per level
    fpath = tmp_path / 'levels.t'
    nctools.ncout3d(str(fpath), lon, lat, p, column, name='t', units='K')
    with Dataset(str(fpath)) as ds:
        np.testing.assert_allclose(np.asarray(ds['t'][:, 0, 0]), column, rtol=1e-6)
        np.testing.assert_allclose(np.asarray(ds['t'][:, 0, 1]), column, rtol=1e-6)
        np.testing.assert_allclose(np.asarray(ds['plev'][:]), p, rtol=1e-6)


# ---------------------------------------------------------------------------
# ncout_opt_prop
# ---------------------------------------------------------------------------


def test_ncout_opt_prop_full_field_writes_and_sorts(tmp_path):
    """ncout_opt_prop writes full absorption/scattering fields and a scalar phf.

    Full (band * level * lat * lon) absorption and scattering arrays are placed
    directly on the grid, a scalar phase function is broadcast, and an unsorted
    pressure axis is reordered to ascend. The stored fields must be finite and
    within the input range, the phase function uniform at the scalar value (the
    broadcast edge case), and the pressure axis ascending after the sort.
    """
    lon = np.array([0.0, 120.0])
    lat = np.array([0.0])
    p = np.array([1.0e5, 1.0e3, 1.0e4])  # unsorted
    bands = 2
    levels = p.size
    size = lon.size * lat.size * levels * bands
    absp = np.linspace(1.0, 2.0, size)
    scat = np.linspace(0.1, 0.2, size)
    fpath = tmp_path / 'field.op_soot'
    nctools.ncout_opt_prop(str(fpath), lon, lat, p, bands, absp, scat, 0.5)
    with Dataset(str(fpath)) as ds:
        abs_field = np.asarray(ds['abs'][:])
        assert abs_field.shape == (bands, levels, 1, 2)
        assert np.all(np.isfinite(abs_field))
        assert abs_field.min() >= 1.0 - 1e-9 and abs_field.max() <= 2.0 + 1e-9
        np.testing.assert_allclose(np.asarray(ds['phf'][:]), 0.5, rtol=1e-6)
        plev = np.asarray(ds['plev'][:])
        assert np.all(np.diff(plev) > 0.0)


def test_ncout_opt_prop_per_band_level_arrays(tmp_path):
    """ncout_opt_prop tiles per-(band, level) properties across the grid.

    When absorption and scattering are given as (band, level) arrays they are
    replicated over every (lat, lon) cell. With an already-sorted pressure axis
    (the edge case that skips the reorder) the tiled column at (lat0, lon0) must
    equal the input per-band-level array exactly.
    """
    lon = np.array([0.0, 120.0])
    lat = np.array([0.0])
    p = np.array([1.0e3, 1.0e4, 1.0e5])  # already ascending
    bands = 2
    levels = p.size
    absp = np.arange(bands * levels, dtype=float).reshape(bands, levels)
    scat = (np.arange(bands * levels, dtype=float) + 100.0).reshape(bands, levels)
    fpath = tmp_path / 'perband.op_soot'
    nctools.ncout_opt_prop(str(fpath), lon, lat, p, bands, absp, scat, 0.5)
    with Dataset(str(fpath)) as ds:
        np.testing.assert_allclose(np.asarray(ds['abs'][:, :, 0, 0]), absp, rtol=1e-6)
        np.testing.assert_allclose(np.asarray(ds['scat'][:, :, 0, 0]), scat, rtol=1e-6)
        assert ds.dimensions['band'].size == bands


def test_ncout_opt_prop_full_phase_function(tmp_path):
    """ncout_opt_prop accepts a full-size phase-function array.

    A phase function sized (lon * lat * level * band) exercises the full-field
    branch of the phf dispatch, which reshapes onto the (band, moment, level,
    lat, lon) grid. The stored phase function must have that five-dimensional
    shape with a single moment and stay finite.
    """
    lon = np.array([0.0, 120.0])
    lat = np.array([0.0])
    p = np.array([1.0e3, 1.0e4, 1.0e5])
    bands = 2
    levels = p.size
    size = lon.size * lat.size * levels * bands
    absp = np.linspace(1.0, 2.0, size)
    scat = np.linspace(0.1, 0.2, size)
    phf = np.linspace(0.0, 1.0, size)
    fpath = tmp_path / 'fullphf.op_soot'
    nctools.ncout_opt_prop(str(fpath), lon, lat, p, bands, absp, scat, phf)
    with Dataset(str(fpath)) as ds:
        phf_field = np.asarray(ds['phf'][:])
        assert phf_field.shape == (bands, 1, levels, 1, 2)
        assert np.all(np.isfinite(phf_field))


def test_ncout_opt_prop_rejects_mismatched_arrays(tmp_path):
    """ncout_opt_prop raises for each optical-property array that mis-sizes.

    Absorption, scattering, and phase-function arrays are each checked against
    the (band, level) and full-grid sizes. A mis-sized absorption array raises
    first; with valid absorption, a mis-sized scattering array raises next; with
    both valid, a mis-sized phase function raises. Each guard names its array.
    """
    lon = np.array([0.0, 120.0])
    lat = np.array([0.0])
    p = np.array([1.0e3, 1.0e4, 1.0e5])
    bands = 2
    levels = p.size
    good = np.linspace(1.0, 2.0, lon.size * lat.size * levels * bands)
    fpath = tmp_path / 'badop.op_soot'
    with pytest.raises(RuntimeError, match='absp'):
        nctools.ncout_opt_prop(
            str(fpath), lon, lat, p, bands, np.arange(5, dtype=float), good, 0.5
        )
    with pytest.raises(RuntimeError, match='scat'):
        nctools.ncout_opt_prop(
            str(fpath), lon, lat, p, bands, good, np.arange(5, dtype=float), 0.5
        )
    with pytest.raises(RuntimeError, match='phf'):
        nctools.ncout_opt_prop(
            str(fpath), lon, lat, p, bands, good, good, np.arange(5, dtype=float)
        )


# ---------------------------------------------------------------------------
# ncout_view
# ---------------------------------------------------------------------------


def test_ncout_view_reaches_undefined_level_axis(tmp_path):
    """ncout_view processes its viewing angles but cannot write the level axis.

    ncout_view broadcasts scalar polar/azimuth/level angles (and, in the array
    form, stores them directly), then attempts to write a 'level' dimension.
    That write references an undefined level count, so a fully valid call
    reaches it and raises NameError. Both the scalar-input path and the
    array-input path (the edge case exercising the exact-size branches) hit the
    same undefined-axis failure.
    """
    lon = np.array([0.0, 120.0])
    lat = np.array([0.0])
    with pytest.raises(NameError, match='levels'):
        nctools.ncout_view(
            str(tmp_path / 'scalar.view'), lon, lat, np.array([1]), np.array([1]), 0.5, 0.5, 0.5
        )
    # Array-sized polar/azimuth/level inputs take the exact-size branches and
    # still reach the same undefined level axis.
    direction = np.array([1, 2])
    level = np.array([1, 2])
    n_full = lon.size * lat.size * direction.size
    with pytest.raises(NameError, match='levels'):
        nctools.ncout_view(
            str(tmp_path / 'array.view'),
            lon,
            lat,
            direction,
            level,
            np.zeros(n_full),
            np.zeros(n_full),
            np.zeros(2),
        )


def test_ncout_view_rejects_mismatched_angles(tmp_path):
    """ncout_view raises when a viewing-angle array matches no expected size.

    Each of the polar, azimuth, and viewing-level arrays must be a scalar or
    match its expected length. A polar array of an intermediate length raises
    first; with a scalar polar angle, a mis-sized azimuth array raises; with
    both scalar, a mis-sized viewing-level array raises. Each guard reports the
    ncout_view size mismatch.
    """
    lon = np.array([0.0, 120.0])
    lat = np.array([0.0])
    direction = np.array([1])
    level = np.array([1, 2])
    base = str(tmp_path / 'bad.view')
    with pytest.raises(RuntimeError, match='ncout_view'):
        nctools.ncout_view(
            base, lon, lat, direction, level, np.arange(3, dtype=float), 0.5, 0.5
        )
    with pytest.raises(RuntimeError, match='ncout_view'):
        nctools.ncout_view(
            base, lon, lat, direction, level, 0.5, np.arange(3, dtype=float), 0.5
        )
    with pytest.raises(RuntimeError, match='ncout_view'):
        nctools.ncout_view(
            base, lon, lat, direction, level, 0.5, 0.5, np.arange(5, dtype=float)
        )
