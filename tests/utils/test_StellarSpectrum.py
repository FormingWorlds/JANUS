"""Tests for src/janus/utils/StellarSpectrum.py.

Exercises the stellar-spectrum preparation for SOCRATES: unit conversion
and file format of the written spectrum, the down-sampling path with its
bin-count contracts, and the prep_spec insertion wrapper with a
stand-in for the binary. See docs/How-to/test.md.
"""

import os
import shutil
import time
from pathlib import Path
from types import SimpleNamespace

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


FAKE_PREP_SPEC = r"""#!/bin/bash
# Stand-in for SOCRATES prep_spec: reads a name into 132 characters, as prep_spec.f90 does.
{ read -r name; pwd; echo "$name"; ls -a; cat; } > "$FAKE_LOG"
name=${name:0:132}
case "$FAKE_MODE" in
  sleep) exec sleep 600 ;;
  flood) while :; do echo '+++ Invalid input: Please re-enter.' >&2; done ;;
  rc2) echo '*** Error: bad block' >&2; seq 40; exit 2 ;;
  noblock) echo '*** Error: no block' >&2; exit 0 ;;
esac
[ -f "$name" ] || { echo 'A new spectral file will be created.' >&2; exit 0; }
printf '*BLOCK: TYPE =    2: VERSION =    1\n*BLOCK: TYPE =    6: VERSION =    1\n' >> "$name"
echo written >> "${name}_k"
"""


@pytest.fixture
def insert(tmp_path, monkeypatch):
    """Inputs for InsertStellarSpectrum and a fake prep_spec first on PATH."""
    bindir = tmp_path / 'bin'
    bindir.mkdir()
    tool = bindir / 'prep_spec'
    tool.write_text(FAKE_PREP_SPEC)
    tool.chmod(0o755)
    monkeypatch.setenv('PATH', f'{bindir}{os.pathsep}{os.environ["PATH"]}')
    monkeypatch.setenv('FAKE_LOG', str(tmp_path / 'fake.log'))
    monkeypatch.setenv('FAKE_MODE', 'ok')
    orig = tmp_path / 'orig.sf'
    orig.write_text('SPECTRAL DATA\n')
    (tmp_path / 'orig.sf_k').write_text('K DATA\n')
    star = tmp_path / 'star.txt'
    star.write_text('Star\n*BEGIN_DATA\n 1.0e-7 1.0e6\n*END\n')
    out = tmp_path / 'out'
    out.mkdir()
    return SimpleNamespace(
        orig=str(orig), star=str(star), out=out, log=tmp_path / 'fake.log', tmp=tmp_path
    )


def _folder(root, length):
    """Make a folder whose path, with its trailing slash, has exactly ``length`` characters."""
    path = str(root)
    path = os.path.join(path, 'd' * (length - 2 - len(path)))
    assert len(path) + 1 == length
    os.makedirs(path, exist_ok=True)
    return path


def test_insert_gives_prep_spec_short_names_and_a_lock(insert):
    """prep_spec reads names into 132 characters and loops without a lock file.

    Whatever the length of the output and star paths, it must see relative names in a
    working folder that holds radiation_code.lock, and the answers must stay in order.
    """
    out = _folder(insert.tmp, 180)
    star = os.path.join(_folder(insert.tmp / 'stars', 170), 'star.txt')
    shutil.copyfile(insert.star, star)
    ss.InsertStellarSpectrum(insert.orig, star, out)

    cwd, first, *rest = insert.log.read_text().splitlines()
    assert first == 'star.sf' and len(first) <= 132
    assert 'radiation_code.lock' in rest and 'star.sf_k' in rest
    assert cwd.startswith(os.path.realpath(out)) and cwd != os.path.realpath(out)
    answers = rest[rest.index('a') :]
    assert answers[:6] == ['a', '6', 'n', 'T', '100 4000', '250']
    assert answers[6:9] == ['2', 'n', 'star.txt']
    assert answers[9:11] == ['y', '-1']
    assert Path(out, 'star.sf').read_text().endswith('VERSION =    1\n')
    assert Path(out, 'star.sf_k').read_text() == 'K DATA\nwritten\n'
    assert sorted(os.listdir(out)) == ['star.sf', 'star.sf_k']


def test_insert_replaces_stale_outputs_and_resolves_a_relative_star(insert, monkeypatch):
    """Old outputs go; a relative star path and PATH entry still work after the cwd change."""
    (insert.out / 'star.sf').write_text('stale')
    (insert.out / 'star.sf_k').write_text('stale')
    monkeypatch.chdir(insert.tmp)
    monkeypatch.setenv('PATH', f'bin{os.pathsep}{os.environ["PATH"]}')
    ss.InsertStellarSpectrum(insert.orig, 'star.txt', str(insert.out))
    assert (insert.out / 'star.sf').read_text().startswith('SPECTRAL DATA')
    assert 'stale' not in (insert.out / 'star.sf_k').read_text()


@pytest.mark.parametrize('length, stops', [(190, False), (191, True)])
def test_insert_stops_on_a_folder_above_the_socrates_limit(insert, length, stops):
    """Cl_run_cdf holds the folder and mcica_data in 200 bytes: 190 passes, 191 stops."""
    out = _folder(insert.tmp, length) + '/'
    Path(out, 'star.sf').write_text('old')
    if stops:
        with pytest.raises(ValueError, match=f'{length} bytes'):
            ss.InsertStellarSpectrum(insert.orig, insert.star, out)
        assert not insert.log.exists(), 'prep_spec must not start'
        assert Path(out, 'star.sf').read_text() == 'old', 'no side effect before the stop'
    else:
        ss.InsertStellarSpectrum(insert.orig, insert.star, out)
        assert Path(out, 'star.sf').read_text().startswith('SPECTRAL DATA')


@pytest.mark.parametrize('body', [None, 'no data marker\n'])
def test_insert_stops_on_a_star_file_prep_spec_would_loop_on(insert, body):
    """Without *BEGIN_DATA in the star file prep_spec reads on without end, even with a lock."""
    star = insert.tmp / 'bad.txt'
    if body is not None:
        star.write_text(body)
    with pytest.raises(ValueError, match='BEGIN_DATA'):
        ss.InsertStellarSpectrum(insert.orig, str(star), str(insert.out))
    assert not insert.log.exists(), 'prep_spec must not start'
    assert os.listdir(insert.out) == []


@pytest.mark.parametrize(
    'mode, message',
    [
        ('rc2', r'code 2: \*\*\* Error: bad block'),
        ('noblock', r'block[\s\S]*no block'),
    ],
)
def test_insert_raises_when_prep_spec_writes_no_star(insert, monkeypatch, mode, message):
    """A return code of 0 proves nothing: prep_spec can skip a block and still exit 0.

    Both failures name the last output, leave no star.sf that a later call could read
    and no working folder.
    """
    (insert.out / 'star.sf').write_text('stale')
    (insert.out / 'star.sf_k').write_text('stale')
    monkeypatch.setenv('FAKE_MODE', mode)
    with pytest.raises(RuntimeError, match=message):
        ss.InsertStellarSpectrum(insert.orig, insert.star, str(insert.out))
    assert os.listdir(insert.out) == []


def test_insert_counts_the_trailing_separator_of_the_folder(insert):
    """PROTEUS passes the folder with a slash and other callers without: both count the same."""
    with pytest.raises(ValueError, match='191 bytes'):
        ss.InsertStellarSpectrum(insert.orig, insert.star, _folder(insert.tmp, 191))
    assert not insert.log.exists()
    ss.InsertStellarSpectrum(insert.orig, insert.star, _folder(insert.tmp / 'b', 190))
    assert insert.log.exists()


def test_insert_counts_bytes_of_a_folder_name_with_non_ascii_characters(insert):
    """The limit in the Fortran string is in bytes: 190 characters here are above 190 bytes."""
    out = insert.tmp / ('\u00e9' * (188 - len(str(insert.tmp))))
    out.mkdir()
    assert len(os.path.join(out, '')) == 190
    with pytest.raises(ValueError, match=' bytes'):
        ss.InsertStellarSpectrum(insert.orig, insert.star, str(out))
    assert not insert.log.exists()


def test_insert_accepts_a_binned_star_file(insert, tmp_path):
    """prep_spec reads both data markers, *BEGIN_DATA and *BEGIN_BINNED_DATA."""
    star = tmp_path / 'binned.txt'
    star.write_text('Star\n*BEGIN_BINNED_DATA\n 1.0e-7 2.0e-7 1.0e6\n*END\n')
    ss.InsertStellarSpectrum(insert.orig, str(star), str(insert.out))
    assert (insert.out / 'star.sf').read_text().count('*BLOCK') == 2
    assert insert.log.exists()


def test_insert_needs_both_blocks(insert, monkeypatch):
    """A file with the thermal block and no stellar block fails, and the other way round."""
    for line, missing in (('6', '2'), ('2', '6')):
        tool = insert.tmp / 'bin' / 'prep_spec'
        tool.write_text(
            '#!/bin/bash\ncat > /dev/null\n'
            f'printf "*BLOCK: TYPE =    {line}: X\\n*BLOCK: TYPE =   1{missing}: X\\n"'
            ' >> star.sf\n'
        )
        with pytest.raises(RuntimeError, match='block'):
            ss.InsertStellarSpectrum(insert.orig, insert.star, str(insert.out))


@pytest.mark.parametrize('mode', ['sleep', 'flood'])
def test_insert_stops_a_stuck_prep_spec_without_growth(insert, monkeypatch, mode):
    """A hang or an endless error loop ends at the time limit and keeps only the last lines."""
    monkeypatch.setenv('FAKE_MODE', mode)
    monkeypatch.setattr(ss, 'PREP_SPEC_TIMEOUT', 1.0)
    start = time.monotonic()
    with pytest.raises(RuntimeError, match='time limit') as err:
        ss.InsertStellarSpectrum(insert.orig, insert.star, str(insert.out))
    assert time.monotonic() - start < 20.0
    assert len(str(err.value)) < 5000
    assert os.listdir(insert.out) == []


def test_insert_names_socrates_when_prep_spec_is_absent(insert, monkeypatch):
    """A SOCRATES tree without the tool gives a stop that names RAD_DIR, before any change."""
    monkeypatch.setenv('PATH', str(insert.tmp / 'empty'))
    (insert.out / 'star.sf').write_text('old')
    with pytest.raises(RuntimeError, match='RAD_DIR'):
        ss.InsertStellarSpectrum(insert.orig, insert.star, str(insert.out))
    assert (insert.out / 'star.sf').read_text() == 'old'
