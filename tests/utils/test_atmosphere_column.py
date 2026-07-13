"""Tests for src/janus/utils/atmosphere_column.py.

atmosphere_column.py defines the ``atmos`` class, the container that carries
the evolving 1D column state (pressure, temperature, composition, radiation
fluxes) and the planetary constants. This file exercises:

* Constructor input validation: exactly one of ``vol_mixing`` / ``vol_partial``
  must be supplied, the surface pressure must be positive, and the core
  volatiles H2O/CO2/N2 must all be present.
* Partial-pressure initialisation: the surface pressure is the sum of the
  partial pressures and the mixing ratios normalise to unity.
* Surface-gravity construction (inverse-square GM/r^2), the pressure/band grid
  invariants, and the no-band and cloud-enabled branches.
* The in-place setters (surface temperature/pressure, planet properties,
  volatiles, radiative-skin tropopause temperature).
* The PT-profile text writer (pressure-unit scaling, descending order, and the
  unknown-unit error contract).
* The netCDF writer round-trip: values, dimensions, units, and metadata survive
  a write-then-read cycle, and the hydrostatic blow-up flag is recorded.

atmosphere_column.py is a state container rather than one of the six physics
sources, so most assertions pin structure and round-trips; tests that assert a
genuine physical invariant (inverse-square gravity, the radiative-skin law)
carry the physics_invariant marker. See docs/How-to/test.md.
"""

import os
from unittest import mock

import numpy as np
import pytest
from netCDF4 import Dataset

import janus.utils.phys as phys
from janus import __version__
from janus.utils.atmosphere_column import atmos

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Three ascending band edges in nm give two bands; a fourth widens the grid.
BAND_EDGES = np.array([100.0, 500.0, 1000.0, 10000.0])  # nm
VOL_MIXING = {'H2O': 0.5, 'CO2': 0.3, 'N2': 0.2}

M_EARTH = 5.972e24
R_EARTH = 6.371e6


def _make_atmos(**overrides):
    """Build a valid ``atmos`` object with Earth-like defaults for testing."""
    kwargs = dict(
        T_surf=1200.0,
        P_surf=1.0e5,
        P_top=1.0,
        pl_radius=R_EARTH,
        pl_mass=M_EARTH,
        band_edges=BAND_EDGES,
        vol_mixing=dict(VOL_MIXING),
    )
    kwargs.update(overrides)
    return atmos(**kwargs)


def _prepare_small_column(atm, ncell, nbands):
    """Resize an ``atmos`` object's arrays to a small synthetic column.

    write_ncdf reads cell-centre arrays of length ``ncell``, cell-edge arrays
    of length ``ncell + 1``, and spectral flux arrays of shape (nbands,
    ncell+1). Production regrids the 10000-level adiabat down to the save grid
    before writing; this mirrors that end state at a tiny size so the netCDF
    write stays fast.
    """
    nedge = ncell + 1
    # Input convention for integrate_heights: index 0 = top of atmosphere
    # (low pressure), index -1 = surface, so the pressure array ascends.
    atm.p = np.logspace(3, 5, ncell)  # Pa, ascending toward the surface
    atm.pl = np.logspace(3, 5, nedge)
    atm.tmp = np.linspace(1200.0, 400.0, ncell)  # K
    atm.tmpl = np.linspace(1200.0, 400.0, nedge)
    atm.mu = np.full(ncell, 0.029)  # kg/mol, air-like mean molar mass
    atm.grav_z = np.full(ncell, atm.grav_s)
    atm.cp = np.full(ncell, 1000.0)  # J/(kg K)
    atm.re = np.zeros(ncell)
    atm.lwm = np.zeros(ncell)
    atm.clfr = np.zeros(ncell)
    atm.net_heating = np.zeros(ncell)
    for vol in atm.vol_list:
        atm.x_gas[vol] = np.full(ncell, atm.vol_list[vol])
        atm.x_cond[vol] = np.zeros(ncell)
        atm.p_vol[vol] = np.full(ncell, atm.ps * atm.vol_list[vol])
        atm.pl_vol[vol] = np.full(nedge, atm.ps * atm.vol_list[vol])
    for name in (
        'LW_flux_down',
        'LW_flux_net',
        'SW_flux_up',
        'SW_flux_down',
        'SW_flux_net',
        'flux_up_total',
        'flux_down_total',
        'net_flux',
    ):
        setattr(atm, name, np.zeros(nedge))
    # Distinctive up-LW profile so the round-trip check discriminates a
    # transposed or mis-sliced write.
    atm.LW_flux_up = np.linspace(300.0, 50.0, nedge)
    for name in (
        'LW_spectral_flux_up',
        'LW_spectral_flux_down',
        'SW_spectral_flux_up',
        'SW_spectral_flux_down',
    ):
        setattr(atm, name, np.zeros((nbands, nedge)))


# ---------------------------------------------------------------------------
# Constructor validation
# ---------------------------------------------------------------------------


def test_atmos_requires_exactly_one_volatile_specification():
    """The constructor rejects both zero and two volatile specifications.

    Exactly one of ``vol_mixing`` / ``vol_partial`` defines the composition.
    Passing neither is under-determined and passing both is contradictory;
    both are the constructor's documented error contract. The two messages
    differ ("Neither" vs "Both"), which discriminates the two guard branches.
    """
    with pytest.raises(Exception, match='Neither'):
        atmos(1200.0, 1.0e5, 1.0, R_EARTH, M_EARTH, BAND_EDGES, vol_mixing={}, vol_partial={})
    with pytest.raises(Exception, match='Both'):
        atmos(
            1200.0,
            1.0e5,
            1.0,
            R_EARTH,
            M_EARTH,
            BAND_EDGES,
            vol_mixing={'H2O': 1.0},
            vol_partial={'H2O': 1.0e5},
        )


def test_atmos_rejects_nonpositive_surface_pressure():
    """A non-positive surface pressure with mixing-ratio input is rejected.

    When the composition is set by mixing ratio the surface pressure is taken
    directly and must be strictly positive (P > 0 Pa everywhere is a hard
    invariant). Zero pressure sits exactly on the boundary and a negative
    pressure is unphysical; both must raise rather than reaching the grid.
    """
    with pytest.raises(Exception, match='positive'):
        _make_atmos(P_surf=0.0)  # boundary: exactly zero
    with pytest.raises(Exception, match='positive'):
        _make_atmos(P_surf=-1.0e4)  # unphysical negative


def test_atmos_requires_core_volatiles():
    """Missing any of the required H2O/CO2/N2 volatiles raises.

    The radiative-convective solver assumes the three core volatiles are
    present. A composition that omits one (here N2) is incomplete and must
    raise rather than silently proceeding with a partial gas list.
    """
    with pytest.raises(Exception, match='required volatiles'):
        _make_atmos(vol_mixing={'H2O': 0.6, 'CO2': 0.4})
    # A single-species (H2O-only) column is likewise rejected: an edge case
    # that also omits CO2 and N2.
    with pytest.raises(Exception, match='required volatiles'):
        _make_atmos(vol_mixing={'H2O': 1.0})


def test_atmos_partial_pressures_set_surface_and_normalise():
    """Partial-pressure input sets P_surf as the sum and normalises ratios.

    When the composition is given as partial pressures, the surface pressure
    is their sum and each mixing ratio is the partial pressure over that sum,
    so the ratios add to one. This is the vol_partial branch, distinct from the
    mixing-ratio branch. Uses an unequal split so a dropped normalisation would
    leave the ratios summing to something other than unity.
    """
    partials = {'H2O': 3.0e4, 'CO2': 1.0e4, 'N2': 1.0e4}  # Pa
    atm = atmos(
        1200.0, 9.9e9, 1.0, R_EARTH, M_EARTH, BAND_EDGES, vol_mixing={}, vol_partial=partials
    )
    # P_surf is overridden by the sum of partials, NOT the 9.9e9 placeholder.
    assert atm.ps == pytest.approx(5.0e4, rel=1e-12)
    assert atm.vol_list['H2O'] == pytest.approx(0.6, rel=1e-12)
    ratio_sum = sum(atm.vol_list[v] for v in ('H2O', 'CO2', 'N2'))
    assert ratio_sum == pytest.approx(1.0, rel=1e-12)
    # Surface partial pressure of each gas is reconstructed as ps * ratio.
    assert atm.p_vol['H2O'][0] == pytest.approx(3.0e4, rel=1e-12)


# ---------------------------------------------------------------------------
# Grid, gravity, and optional-branch construction
# ---------------------------------------------------------------------------


@pytest.mark.physics_invariant
def test_atmos_construction_sets_gravity_and_band_grid():
    """Construction fixes inverse-square surface gravity and a valid band grid.

    Surface gravity is GM/r^2, positive and matching the Earth reference. The
    spectral grid derived from N band edges must hold N-1 positive band widths
    and monotone-increasing band centres. The surface pressure/temperature seed
    the base of the column. Doubling the planet mass at construction must
    double gravity, which discriminates a wrong power on the mass.
    """
    atm = _make_atmos()
    assert atm.grav_s == pytest.approx(phys.G * M_EARTH / R_EARTH**2, rel=1e-12)
    assert atm.grav_s > 0.0
    heavier = _make_atmos(pl_mass=2.0 * M_EARTH)
    assert heavier.grav_s == pytest.approx(2.0 * atm.grav_s, rel=1e-12)
    # Band grid: 4 edges -> 3 bands, all widths positive, centres ascending.
    assert atm.nbands == len(BAND_EDGES) - 1
    assert np.all(atm.band_widths > 0.0)
    assert np.all(np.diff(atm.band_centres) > 0.0)
    # Column base seeded with the surface state.
    assert atm.p[0] == pytest.approx(1.0e5, rel=1e-12)
    assert atm.tmp[0] == pytest.approx(1200.0, rel=1e-12)


def test_atmos_without_band_edges_skips_radiation_arrays():
    """A single band edge leaves the spectral grid and flux arrays unset.

    ``bands_set`` is True only when more than one band edge is supplied. With a
    single edge (a degenerate, no-band column) the constructor must skip the
    radiation-array allocation entirely; querying a flux array or the band
    count then raises AttributeError. This is the false branch of the band
    guard.
    """
    atm = _make_atmos(band_edges=np.array([100.0]))
    assert atm.bands_set is False
    assert not hasattr(atm, 'LW_flux_up')
    assert not hasattr(atm, 'nbands')


def test_atmos_cloud_branch_allocates_cloud_state():
    """Enabling clouds allocates the cloud-scheme arrays and retains droplets.

    With ``do_cloud`` True the constructor selects the mixed-overlap cloud
    scheme and allocates per-level effective-radius, liquid-water, and cloud-
    fraction arrays. When the condensate-retention fraction ``alpha_cloud`` is
    above the rain-out floor the droplet inputs are kept rather than zeroed,
    which is the alpha_cloud >= 1e-20 branch. Uses a non-trivial droplet radius
    so a stray zeroing would be caught.
    """
    atm = _make_atmos(do_cloud=True, alpha_cloud=0.5, re=1.2e-5, lwm=3.0e-4, clfr=0.4)
    assert atm.cloud_scheme == 2  # ip_cloud_mix_max
    assert len(atm.re) == atm.nlev_save
    assert len(atm.clfr) == atm.nlev_save
    # alpha_cloud above the floor keeps the droplet inputs.
    assert atm.effective_radius == pytest.approx(1.2e-5, rel=1e-12)
    assert atm.cloud_fraction == pytest.approx(0.4, rel=1e-12)


# ---------------------------------------------------------------------------
# In-place setters
# ---------------------------------------------------------------------------


def test_set_surface_temperature_and_pressure_update_base():
    """The surface setters push the new value into the base of the column.

    setSurfaceTemperature and setSurfacePressure update the scalar surface
    state and the first entry of the centre and edge arrays together, so the
    base of the column stays consistent. Uses distinct new values so a setter
    that forgot the edge array would leave tmpl[0] / pl[0] stale.
    """
    atm = _make_atmos()
    atm.setSurfaceTemperature(1850.0)
    assert atm.ts == pytest.approx(1850.0, rel=1e-12)
    assert atm.tmp[0] == pytest.approx(1850.0, rel=1e-12)
    assert atm.tmpl[0] == pytest.approx(1850.0, rel=1e-12)
    atm.setSurfacePressure(2.5e5)
    assert atm.ps == pytest.approx(2.5e5, rel=1e-12)
    assert atm.p[0] == pytest.approx(2.5e5, rel=1e-12)
    assert atm.pl[0] == pytest.approx(2.5e5, rel=1e-12)


@pytest.mark.physics_invariant
def test_set_planet_properties_recomputes_inverse_square_gravity():
    """Updating planet mass and radius recomputes GM/r^2 surface gravity.

    setPlanetProperties overwrites the mass and radius and recomputes surface
    gravity, propagating it into the base of the gravity profile. Halving the
    radius must quadruple gravity (inverse-square), which discriminates an
    inverse-linear bug; the recomputed value must also match GM/r^2 exactly.
    """
    atm = _make_atmos()
    g_before = atm.grav_s
    atm.setPlanetProperties(pl_radius=R_EARTH / 2.0, pl_mass=M_EARTH)
    assert atm.grav_s == pytest.approx(phys.G * M_EARTH / (R_EARTH / 2.0) ** 2, rel=1e-12)
    assert atm.grav_s == pytest.approx(4.0 * g_before, rel=1e-12)
    assert atm.grav_s > 3.0 * g_before  # rules out the inverse-linear 2x result
    assert atm.grav_z[0] == pytest.approx(atm.grav_s, rel=1e-12)


def test_set_volatiles_normalises_and_floors_water():
    """setVolatiles renormalises the mixing ratios and floors water.

    setVolatiles accepts unnormalised mixing ratios, divides by their sum so
    they add to one, floors H2O to 1e-20 to avoid NaNs, and refreshes the
    surface partial pressures. Passing an unnormalised set with zero water
    exercises both the renormalisation and the water floor edge case.
    """
    atm = _make_atmos()
    atm.setVolatiles({'H2O': 0.0, 'CO2': 6.0, 'N2': 4.0})  # sum 10, not 1
    assert atm.vol_list['CO2'] == pytest.approx(0.6, rel=1e-12)
    assert atm.vol_list['N2'] == pytest.approx(0.4, rel=1e-12)
    # Water floored to 1e-20 rather than left at zero.
    assert atm.vol_list['H2O'] == pytest.approx(1e-20, rel=1e-6)
    # Surface partial pressure refreshed from the new ratio.
    assert atm.p_vol['CO2'][0] == pytest.approx(atm.ps * 0.6, rel=1e-12)


@pytest.mark.physics_invariant
def test_set_tropopause_temperature_radiative_skin():
    """The tropopause temperature is the radiative-skin fraction of equilibrium.

    setTropopauseTemperature computes the planetary equilibrium temperature
    from the instellation and Bond albedo, then scales it by 0.5**0.25 (the
    grey radiative-skin factor). The result must be positive and strictly below
    the equilibrium temperature; the 0.5**0.25 factor (~0.841) discriminates a
    dropped skin scaling. A zero-instellation column has a zero skin
    temperature, the limit case.
    """
    atm = _make_atmos()
    atm.instellation = 1361.0
    atm.setTropopauseTemperature()
    t_eqm = (atm.instellation * atm.inst_sf * (1.0 - atm.albedo_pl) / phys.sigma) ** 0.25
    assert atm.trppT == pytest.approx(t_eqm * 0.5**0.25, rel=1e-12)
    assert atm.trppT > 0.0
    assert atm.trppT < t_eqm  # radiative skin is cooler than equilibrium
    # Limit case: no instellation -> no skin temperature.
    atm.instellation = 0.0
    atm.setTropopauseTemperature()
    assert atm.trppT == pytest.approx(0.0, abs=1e-30)


# ---------------------------------------------------------------------------
# PT-profile text writer
# ---------------------------------------------------------------------------


def test_write_pt_applies_pressure_unit_scaling(tmp_path):
    """write_PT scales pressure by the requested unit and descends in pressure.

    write_PT saves a two-column pressure/temperature table with the pressure
    converted to the requested unit and the rows reversed. The bar (1e-5) and
    dyne/cm2 (1e5) conversions sit ten orders of magnitude apart, so a wrong
    scale factor for either is far outside tolerance. The temperature column is
    unit-invariant and must round-trip unchanged.
    """
    atm = _make_atmos()
    atm.p = np.array([1.0e5, 5.0e4, 1.0e4])  # Pa, descending
    atm.tmp = np.array([1000.0, 600.0, 300.0])  # K
    scales = {'Pa': 1.0, 'bar': 1.0e-5, 'dyne/cm2': 1.0e5, 'atm': 1.01325e5}
    for unit, factor in scales.items():
        fpath = tmp_path / f'pt_{unit.replace("/", "_")}.tsv'
        atm.write_PT(filename=str(fpath), punit=unit)
        data = np.loadtxt(str(fpath), skiprows=2)
        expected_p = (atm.p * factor)[::-1]  # source multiplies Pa by the scale factor
        np.testing.assert_allclose(data[:, 0], expected_p, rtol=1e-4)
        np.testing.assert_allclose(data[:, 1], atm.tmp[::-1], rtol=1e-6)


def test_write_pt_rejects_unknown_pressure_unit(tmp_path):
    """write_PT raises on an unrecognised pressure unit instead of writing.

    The pressure-unit dispatch has a closed set of options; an unknown unit
    ("kbar" here) falls through to the default guard and raises rather than
    writing a file with an undefined scale. No file must be produced. A valid
    unit on the same object still writes, confirming the guard is specific to
    the bad unit and not a general failure.
    """
    atm = _make_atmos()
    atm.p = np.array([1.0e5, 1.0e4])
    atm.tmp = np.array([900.0, 300.0])
    bad = tmp_path / 'bad.tsv'
    with pytest.raises(Exception, match='Unrecognised pressure unit'):
        atm.write_PT(filename=str(bad), punit='kbar')
    assert not bad.exists()
    good = tmp_path / 'good.tsv'
    atm.write_PT(filename=str(good), punit='Pa')
    assert good.exists()


# ---------------------------------------------------------------------------
# netCDF writer
# ---------------------------------------------------------------------------


def test_write_ncdf_round_trips_column(tmp_path, monkeypatch):
    """write_ncdf serialises the column and metadata to a readable netCDF file.

    The writer records the SOCRATES version (read from RAD_DIR/version), the
    JANUS version, the scalar surface state, the per-level pressure/temperature
    and flux profiles, and the band-edge grid converted from nm to m. Reading
    the file back must recover every value, the units attributes, and the
    dimension sizes. Writing twice to the same path exercises the pre-existing-
    file removal branch, and enabling the contribution-function output adds the
    contfunc variable (an optional-branch edge case). The nm->m band-edge
    factor (1e-9) discriminates a dropped unit conversion.
    """
    (tmp_path / 'version').write_text('sfx-test-1\n')
    monkeypatch.setenv('RAD_DIR', str(tmp_path))
    atm = _make_atmos()
    ncell = 6
    _prepare_small_column(atm, ncell, atm.nbands)
    atm.instellation = 1361.0
    atm.trppP = 1234.0
    atm.has_contfunc = True
    atm.cff = np.zeros((atm.nbands, ncell))
    fpath = tmp_path / 'atm.nc'
    atm.write_ncdf(str(fpath))
    atm.write_ncdf(str(fpath))  # second write removes and rebuilds the file
    with Dataset(str(fpath)) as ds:
        assert ds.dimensions['nlev_c'].size == ncell
        assert ds.dimensions['nlev_l'].size == ncell + 1
        assert ds.dimensions['nbands'].size == atm.nbands
        np.testing.assert_allclose(ds['p'][:], atm.p, rtol=1e-5)
        np.testing.assert_allclose(ds['tmp'][:], atm.tmp, rtol=1e-5)
        np.testing.assert_allclose(ds['fl_U_LW'][:], atm.LW_flux_up, rtol=1e-5)
        np.testing.assert_allclose(ds['bandmin'][:], atm.band_edges[:-1] * 1e-9, rtol=1e-6)
        assert float(ds['tmp_surf'][...]) == pytest.approx(atm.ts, rel=1e-5)
        assert float(ds['instellation'][...]) == pytest.approx(1361.0, rel=1e-5)
        assert ds['p'].units == 'Pa'
        assert ds['fl_U_LW'].units == 'W m-2'
        assert ds.SOCRATES_version == 'sfx-test-1'
        assert ds.JANUS_version == __version__
        assert ds['contfunc'].shape == (atm.nbands, ncell)


def test_write_ncdf_flags_failed_height_integration(tmp_path, monkeypatch):
    """write_ncdf records a hydrostatic blow-up and stores a bounded dummy grid.

    When the hydrostatic integration diverges (an absurdly hot column drives
    the height past the 1e8 m guard), write_ncdf must store height_ok = 'n' and
    a bounded 0..1000 m dummy altitude grid rather than propagating the runaway
    value. This is the error-contract branch of the height integration as seen
    through the writer.
    """
    (tmp_path / 'version').write_text('sfx-test-2\n')
    monkeypatch.setenv('RAD_DIR', str(tmp_path))
    atm = _make_atmos()
    ncell = 6
    _prepare_small_column(atm, ncell, atm.nbands)
    atm.tmp = np.full(ncell, 1.0e12)  # absurd temperature forces divergence
    fpath = tmp_path / 'atm_fail.nc'
    atm.write_ncdf(str(fpath))
    with Dataset(str(fpath)) as ds:
        flag = np.asarray(ds['height_ok'][:]).tobytes()
        assert b'n' in flag
        assert np.all(np.isfinite(ds['z'][:]))
        assert float(np.max(ds['z'][:])) <= 1000.0 + 1e-6


def test_write_ncdf_username_falls_back_when_getlogin_fails(tmp_path, monkeypatch):
    """write_ncdf resolves the username via pwd when os.getlogin raises.

    In detached or containerised runs os.getlogin() raises OSError (there is no
    controlling terminal). write_ncdf must fall back to the password database
    (pwd.getpwuid) rather than propagating the error, and still produce a
    complete, readable file. The recorded username must be a non-empty string
    matching the fallback lookup, which is the error-contract branch of the
    metadata block.
    """
    import pwd

    (tmp_path / 'version').write_text('sfx-test-3\n')
    monkeypatch.setenv('RAD_DIR', str(tmp_path))
    atm = _make_atmos()
    ncell = 5
    _prepare_small_column(atm, ncell, atm.nbands)
    fpath = tmp_path / 'atm_user.nc'
    expected_user = pwd.getpwuid(os.getuid()).pw_name
    with mock.patch('os.getlogin', side_effect=OSError('no controlling terminal')):
        atm.write_ncdf(str(fpath))
    with Dataset(str(fpath)) as ds:
        assert ds.username == expected_user
        assert len(ds.username) > 0
