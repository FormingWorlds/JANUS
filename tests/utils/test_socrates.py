"""Tests for src/janus/utils/socrates.py.

Exercises the SOCRATES driver with the binary, the netCDF writers, and the
result files mocked: input validation guards, the gas-overlap dispatch, the
cloud and verbose call-sequence variants, flux assembly and its
conservation identities, the contribution-function branch, and the output
cleaner. See docs/How-to/test.md.
"""

from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pytest

import janus.utils.socrates as jus

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

NB = 3  # spectral bands
NL = 4  # levels


def _stub_atm(**over):
    """Atmosphere stand-in with physically plausible driver inputs."""
    nlev = NL
    atm = SimpleNamespace(
        bands_set=True,
        nbands=NB,
        cp=np.full(nlev, 1.2e3),  # J K-1 kg-1, air-like
        x_gas={'H2O': np.full(nlev, 0.01), 'N2': np.full(nlev, 0.99)},
        vol_list={'H2O': 1, 'N2': 1},
        overlap_type=2,
        do_cloud=False,
        albedo_s=0.1,
        ts=288.0,
        ps=1.0e5,
        zenith_angle=48.19,
        instellation=1361.0,
        albedo_pl=0.3,
        inst_sf=0.25,
        p=np.logspace(5, 3, nlev),
        pl=np.logspace(5, 3, nlev + 1),
        tmp=np.linspace(288.0, 210.0, nlev),
        tmpl=np.linspace(288.0, 205.0, nlev + 1),
        mu=np.full(nlev, 0.028),  # kg mol-1
        planet_radius=6.371e6,
        grav_s=9.81,
    )
    for k, v in over.items():
        setattr(atm, k, v)
    return atm


def _flux_cube(base):
    """A (band, level, lat, lon) cube with distinct per-band values."""
    cube = np.zeros((NB, NL + 1, 1, 1))
    for b in range(NB):
        cube[b, :, 0, 0] = base * (b + 1) + np.arange(NL + 1)
    return cube


class _FakeDataset:
    def __init__(self, variables):
        self.variables = variables
        self.closed = False

    def close(self):
        self.closed = True


def _install_socrates_mocks(monkeypatch, tmp_path, cubes, with_cff=False):
    """Mock the writers, the binary, and the result datasets in tmp_path."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'star.sf').write_text('SPECTRAL')
    (tmp_path / 'star.sf_k').write_text('K')
    if with_cff:
        (tmp_path / 'currentlw.cff').write_text('marker')

    for writer in ('ncout_surf', 'ncout2d', 'ncout3d'):
        monkeypatch.setattr(jus.nctools, writer, MagicMock(), raising=True)

    run_calls = []

    def fake_run(seq, check=True, stderr=None, stdout=None):
        run_calls.append({'seq': list(seq), 'stdout': stdout})

    monkeypatch.setattr(jus.subprocess, 'run', fake_run)

    datasets = {
        'currentsw.vflx': {'vflx': cubes['vflx']},
        'currentsw.sflx': {'sflx': cubes['vflx'] * 0.5},
        'currentsw.dflx': {'dflx': cubes['vflx'] * 0.5},
        'currentsw.uflx': {'uflx': cubes['uflx_sw']},
        'currentsw.hrts': {'hrts': cubes['hrts_sw']},
        'currentlw.dflx': {'dflx': cubes['dflx_lw']},
        'currentlw.uflx': {'uflx': cubes['uflx_lw']},
        'currentlw.hrts': {'hrts': cubes['hrts_lw']},
        'currentlw.cff': {'cff': cubes['hrts_lw'] * 2.0},
    }
    monkeypatch.setattr(
        jus.net, 'Dataset', lambda fname: _FakeDataset(datasets[fname]), raising=True
    )
    return run_calls


def _cubes():
    return {
        'vflx': _flux_cube(20.0),
        'uflx_sw': _flux_cube(10.0),
        'dflx_lw': _flux_cube(30.0),
        'uflx_lw': _flux_cube(40.0),
        'hrts_sw': _flux_cube(1.0),
        'hrts_lw': _flux_cube(2.0),
    }


def test_radcompsoc_input_guards():
    """Unloaded bands, negative heat capacity, negative mixing ratios, and
    a bad overlap choice all abort before any SOCRATES call.

    The first three are hard physical-validity guards (an unset spectral
    grid, cp <= 0, x < 0 are all unphysical inputs); the overlap check is
    the documented ValueError contract for its 2..8 integer range.
    """
    dirs = {'output': '/nonexistent/'}

    with pytest.raises(Exception, match='bands have not been loaded'):
        jus.radCompSoc(_stub_atm(bands_set=False), dirs, recalc=False)

    with pytest.raises(SystemExit):
        jus.radCompSoc(_stub_atm(cp=np.array([-1.0])), dirs, recalc=False)

    with pytest.raises(SystemExit):
        jus.radCompSoc(_stub_atm(x_gas={'H2O': np.array([-0.1])}), dirs, recalc=False)

    with pytest.raises(ValueError, match='Invalid overlap choice'):
        jus.radCompSoc(_stub_atm(overlap_type=99), dirs, recalc=False)


@pytest.mark.physics_invariant
def test_radcompsoc_flux_assembly_and_conservation(tmp_path, monkeypatch):
    """Band-summed fluxes match the mocked cubes and satisfy the total-flux
    identities.

    The driver must sum each flux cube over bands at every level, and the
    assembled totals must obey flux_up_total = SW_up + LW_up and
    net_flux = up_total - down_total by construction; distinct per-cube
    magnitudes discriminate any up/down or SW/LW mix-up. The overlap value
    8 must map to the SOCRATES flag string '8 0'.
    """
    cubes = _cubes()
    run_calls = _install_socrates_mocks(monkeypatch, tmp_path, cubes)
    atm = _stub_atm(overlap_type=8)
    dirs = {'output': str(tmp_path) + '/'}

    out = jus.radCompSoc(atm, dirs, recalc=False, rscatter=False)

    expected_sw_up = np.sum(cubes['uflx_sw'], axis=0)[:, 0, 0]
    expected_lw_up = np.sum(cubes['uflx_lw'], axis=0)[:, 0, 0]
    expected_sw_dn = np.sum(cubes['vflx'], axis=0)[:, 0, 0]
    expected_lw_dn = np.sum(cubes['dflx_lw'], axis=0)[:, 0, 0]
    np.testing.assert_allclose(out.SW_flux_up, expected_sw_up, rtol=1e-12)
    np.testing.assert_allclose(out.LW_flux_up, expected_lw_up, rtol=1e-12)
    np.testing.assert_allclose(out.SW_flux_down, expected_sw_dn, rtol=1e-12)
    np.testing.assert_allclose(out.LW_flux_down, expected_lw_dn, rtol=1e-12)
    # Conservation identities of the assembled totals.
    np.testing.assert_allclose(out.flux_up_total, expected_sw_up + expected_lw_up, rtol=1e-12)
    np.testing.assert_allclose(
        out.net_flux,
        (expected_sw_up + expected_lw_up) - (expected_sw_dn + expected_lw_dn),
        rtol=1e-12,
    )
    # Heating rates present because recalc is False.
    np.testing.assert_allclose(
        out.SW_heating, np.sum(cubes['hrts_sw'], axis=0)[:, 0, 0], rtol=1e-12
    )
    # TOA heating: instellation scaled by Bond albedo, the instellation
    # scale factor, and the zenith-angle cosine.
    expected_toah = 1361.0 * (1.0 - 0.3) * 0.25 * np.cos(np.deg2rad(48.19))
    assert out.toa_heating == pytest.approx(expected_toah, rel=1e-12)
    assert out.toa_heating > 0.0

    # Four subprocess calls: SW run, SW move, LW run, LW move; overlap 8
    # arrives as the compound flag and output is suppressed by default.
    assert len(run_calls) == 4
    sw_seq = ' '.join(run_calls[0]['seq'])
    assert '-g 8 0' in sw_seq
    assert run_calls[0]['stdout'] is jus.subprocess.DEVNULL


def test_radcompsoc_cloud_verbose_and_recalc_variants(tmp_path, monkeypatch):
    """Cloudy columns write the condensate files and cloud flags, verbose
    mode streams output, and recalc skips the heating rates.

    The cloud branch must write effective radius, liquid water mass, and
    cloud fraction via the netCDF writers and pass the cloud scheme
    numbers to the binary; verbose mode must not attach DEVNULL; a recalc
    call must leave no heating attributes on the returned object.
    """
    cubes = _cubes()
    run_calls = _install_socrates_mocks(monkeypatch, tmp_path, cubes)
    atm = _stub_atm(
        do_cloud=True,
        re=np.full(NL, 1.0e-5),
        lwm=np.full(NL, 1.0e-4),
        clfr=np.full(NL, 0.5),
        cloud_scheme=2,
        cloud_representation=1,
        droplet_type=5,
        solver=16,
    )
    dirs = {'output': str(tmp_path) + '/'}

    out = jus.radCompSoc(atm, dirs, recalc=True, verbose=True)

    # The three condensate profiles were written.
    written_vars = [c.args[5] for c in jus.nctools.ncout3d.call_args_list]
    for varname in ('re', 'lwm', 'clfr'):
        assert varname in written_vars
    # Cloud parameters reach the call sequence and verbose streams output.
    sw_seq = ' '.join(run_calls[0]['seq'])
    assert ' -C ' in sw_seq and ' -K ' in sw_seq
    assert str(atm.solver) in run_calls[0]['seq']
    assert run_calls[0]['stdout'] is None
    # recalc=True: heating rates are not populated.
    assert not hasattr(out, 'SW_heating')
    assert not hasattr(out, 'net_heating')


def test_radcompsoc_contribution_function_branch(tmp_path, monkeypatch):
    """A contribution-function file switches on the cff readback.

    When SOCRATES leaves a currentlw.cff file the driver must flag it,
    load the per-band contribution function, and mirror the spectral LW
    flux; without the file has_contfunc stays False.
    """
    cubes = _cubes()
    _install_socrates_mocks(monkeypatch, tmp_path, cubes, with_cff=True)
    atm = _stub_atm()
    dirs = {'output': str(tmp_path) + '/'}

    out = jus.radCompSoc(atm, dirs, recalc=False)

    assert out.has_contfunc is True
    np.testing.assert_allclose(out.cff, 2.0 * cubes['hrts_lw'][:, :, 0, 0], rtol=1e-12)
    np.testing.assert_allclose(out.LW_flux_up_i, cubes['uflx_lw'][:, :, 0, 0], rtol=1e-12)


def test_natural_sort_and_clean_output_dir(tmp_path):
    """natural_sort orders numerically and CleanOutputDir removes exactly
    the SOCRATES scratch files.

    The cleaner globs current*.* and profile.*; a spectral file and an
    unrelated log must survive while every scratch file disappears. The
    empty-directory call is the degenerate no-op case.
    """
    assert jus.natural_sort(['b10', 'b2', 'a5']) == ['a5', 'b2', 'b10']

    scratch = ['currentsw.vflx', 'currentlw.uflx', 'profile.t', 'profile.nml']
    keep = ['star.sf', 'run.log']
    for name in scratch + keep:
        (tmp_path / name).write_text('x')

    jus.CleanOutputDir(str(tmp_path))
    for name in scratch:
        assert not (tmp_path / name).exists()
    for name in keep:
        assert (tmp_path / name).exists()

    # Second pass on the already-clean directory must be a silent no-op.
    jus.CleanOutputDir(str(tmp_path))
    assert (tmp_path / 'star.sf').exists()
