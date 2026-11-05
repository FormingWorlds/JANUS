"""InsertStellarSpectrum with the real prep_spec.

Each call runs in a child process group under a time limit, with its output dropped, so a
loop in prep_spec fails the test by the limit, without growth on disk. The child imports the
code of this tree and exits with its own code if it did not. See docs/How-to/test.md.
"""

import os
import signal
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from helpers import DownloadSpectralFiles, spectral_file_dir

pytestmark = [pytest.mark.smoke, pytest.mark.timeout(60)]

LIMIT_S = 30
SRC = Path(__file__).resolve().parents[1] / 'src'
NO_STAR_BLOCK = 3
OTHER_TREE = 4

CHILD = f"""
import sys
import janus
from janus.utils import StellarSpectrum as ss
if not janus.__file__.startswith(sys.argv[1]):
    sys.exit({OTHER_TREE})
try:
    ss.InsertStellarSpectrum(*sys.argv[2:])
except RuntimeError as err:
    sys.exit({NO_STAR_BLOCK} if 'wavelengths are incorrectly ordered' in str(err) else 1)
"""


@pytest.fixture(scope='module', autouse=True)
def _oak():
    DownloadSpectralFiles('Oak')


def _insert(star, out):
    """Run the insert on Oak 318 in a child; return its code, or None at the time limit."""
    orig = str(spectral_file_dir('Oak') / 'Oak.sf')
    proc = subprocess.Popen(
        [sys.executable, '-c', CHILD, str(SRC / 'janus'), orig, str(star), str(out)],
        env={**os.environ, 'PYTHONPATH': str(SRC)},
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
    )
    try:
        return proc.wait(LIMIT_S)
    except subprocess.TimeoutExpired:
        return None
    finally:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()


def _star(path, end):
    """Write a star file from 1e-10 m up to ``end``."""
    rows = ''.join(f' {w:.7e} 1.0e+06\n' for w in np.geomspace(1e-10, end, 80))
    path.write_text(f'Star\n WAVELENGTH IRRADIANCE\n (m) (W/m3)\n*BEGIN_DATA\n{rows}*END\n ')
    return path


def test_long_output_path_gives_the_same_spectral_file_as_a_short_one(tmp_path):
    """A 190-byte folder is above the 132 of prep_spec and inside the limit from Cl_run_cdf:
    the call ends and writes the bytes of a short-path call, with the star in them.
    """
    star = _star(tmp_path / 'star.txt', 1.01e-4)
    short = tmp_path / 'short'
    short.mkdir()
    assert _insert(star, short) == 0

    long = str(tmp_path / 'long')
    long = os.path.join(long, 'd' * (188 - len(long)))
    os.makedirs(long)
    assert len(long) + 1 == 190
    assert _insert(star, long) == 0
    for name in ('star.sf', 'star.sf_k'):
        with open(short / name, 'rb') as a, open(os.path.join(long, name), 'rb') as b:
            assert a.read() == b.read()
    text = (short / 'star.sf').read_text()
    assert '*BLOCK: TYPE =    2' in text and '*BLOCK: TYPE =    6' in text
    assert sorted(os.listdir(long)) == ['star.sf', 'star.sf_k']


def test_star_ending_on_an_inner_band_edge_stops_the_insert(tmp_path):
    """A star that ends exactly on an inner band edge, 1e-4 m on Oak 318: prep_spec reports an
    error, exits 0 and leaves out block 2. The error line of prep_spec must reach the message.
    """
    out = tmp_path / 'out'
    out.mkdir()
    assert _insert(_star(tmp_path / 'star.txt', 1e-4), out) == NO_STAR_BLOCK
    assert os.listdir(out) == []


def test_the_time_limit_leaves_the_fetch_out(request):
    """A mirror read may take 60 s by itself, so a limit of 60 s on the fixture ends the fetch
    before it can reach the next mirror. The limit must hold for the test call only.
    """
    assert request.node.get_closest_marker('timeout').kwargs == {'func_only': True}
