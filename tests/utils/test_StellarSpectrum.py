"""Tests for src/janus/utils/StellarSpectrum.py.

Exercises the stellar-spectrum preparation for SOCRATES: unit conversion
and file format of the written spectrum, the down-sampling path with its
bin-count contracts, and the prep_spec insertion wrapper with the binary
mocked. See docs/How-to/test.md.
"""

import numpy as np
import pytest

from janus.utils import StellarSpectrum as ss

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def _read_back(star_file):
    """Parse the written spectrum block into wavelength/flux arrays."""
    lines = star_file.read_text().splitlines()
    i0 = lines.index('*BEGIN_DATA') + 1
    i1 = lines.index('*END')
    rows = [tuple(map(float, ln.split())) for ln in lines[i0:i1]]
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1]


@pytest.mark.physics_invariant
def test_prepare_spectrum_units_format_and_length_contract(tmp_path):
    """The written spectrum converts nm to m and erg s-1 cm-2 nm-1 to W m-3.

    Both conversions are pure powers of ten (1e-9 and 1e6), so reading the
    file back and comparing against the inputs pins them independently; a
    swapped or missing factor moves every row by many orders of magnitude.
    Mismatched input lengths are the documented error contract.
    """
    n = 600  # above the short-spectrum warning threshold
    wl = np.linspace(100.0, 3000.0, n)  # nm
    fl = np.full(n, 2.5)  # erg s-1 cm-2 nm-1
    star = tmp_path / 'star.txt'
    ss.PrepareStellarSpectrum(wl, fl, str(star))

    wl_m, fl_w = _read_back(star)
    assert wl_m.size == n
    np.testing.assert_allclose(wl_m, wl * 1e-9, rtol=1e-6)
    np.testing.assert_allclose(fl_w, fl * 1e6, rtol=1e-6)
    # Scale guards: metre-scale wavelengths for nm input, positive flux.
    assert np.all(wl_m < 1e-5) and np.all(wl_m > 1e-8)
    assert np.all(fl_w > 0.0)

    with pytest.raises(Exception, match='different lengths'):
        ss.PrepareStellarSpectrum(wl, fl[:-1], str(star))


def test_prepare_spectrum_rebins_long_input(tmp_path, caplog):
    """An overlong spectrum is down-sampled onto a log-wavelength grid.

    Input beyond the SOCRATES bin ceiling must come out with exactly
    nbins_max rows, preserve the wavelength span, and reproduce a smooth
    input flux after interpolation; the small-bin-count warning fires for
    a target under 500 bins.
    """
    n = 100_100  # beyond the 99999-bin ceiling
    wl = np.logspace(2.0, 4.0, n)  # nm
    fl = np.full(n, 3.0)  # constant flux is reproduced exactly by pchip
    star = tmp_path / 'rebinned.txt'
    ss.PrepareStellarSpectrum(wl, fl, str(star), nbins_max=400)

    assert 'small' in caplog.text.lower()
    wl_m, fl_w = _read_back(star)
    assert wl_m.size == 400
    assert wl_m[0] == pytest.approx(wl[0] * 1e-9, rel=1e-6)
    assert wl_m[-1] == pytest.approx(wl[-1] * 1e-9, rel=1e-6)
    assert np.all(np.diff(wl_m) > 0.0)
    np.testing.assert_allclose(fl_w, 3.0e6, rtol=1e-8)

    # Requesting more bins than SOCRATES accepts is the other contract.
    with pytest.raises(Exception, match='Too many bins'):
        ss.PrepareStellarSpectrum(wl, fl, str(star), nbins_max=100_000)


def test_insert_stellar_spectrum_replaces_outputs_and_calls_prep_spec(
    tmp_path, monkeypatch, caplog
):
    """Insertion copies both spectral files, replaces stale outputs, and
    feeds prep_spec the expected script.

    Pre-existing outputs must be overwritten by the fresh copies; the
    prep_spec stdin script must reference the star file between the insert
    commands; a nonzero return code logs a warning instead of raising.
    """
    orig = tmp_path / 'orig.sf'
    orig.write_text('SPECTRAL DATA')
    (tmp_path / 'orig.sf_k').write_text('K DATA')
    star = tmp_path / 'star.txt'
    star.write_text('STAR')
    outdir = tmp_path / 'out'
    outdir.mkdir()
    # Stale outputs that must be replaced.
    (outdir / 'star.sf').write_text('stale')
    (outdir / 'star.sf_k').write_text('stale')

    calls = {}

    class _Proc:
        returncode = 1

    def fake_run(cmd, stdout=None, input=None, encoding=None):
        calls['cmd'] = cmd
        calls['input'] = input
        return _Proc()

    monkeypatch.setattr(ss.subprocess, 'run', fake_run)
    ss.InsertStellarSpectrum(str(orig), str(star), str(outdir))

    assert (outdir / 'star.sf').read_text() == 'SPECTRAL DATA'
    assert (outdir / 'star.sf_k').read_text() == 'K DATA'
    assert calls['cmd'] == ['prep_spec']
    assert str(star) in calls['input']
    assert 'prep_spec returned with code 1' in caplog.text
