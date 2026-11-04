"""InsertStellarSpectrum with the real prep_spec on a path above its 132-character limit.

The call runs in its own process group under a time limit, with its output dropped, so a
loop in prep_spec fails this test by the limit, without growth on disk. See docs/How-to/test.md.
"""

import os
import shutil
import signal
import subprocess
import sys
import textwrap

import pytest
from helpers import DownloadSpectralFiles, spectral_file_dir

pytestmark = [pytest.mark.smoke, pytest.mark.timeout(60)]

LIMIT_S = 30

CHILD = textwrap.dedent(
    """
    import sys
    from janus.utils import StellarSpectrum as ss
    try:
        ss.InsertStellarSpectrum(*sys.argv[1:])
    except ValueError:
        sys.exit(STOPPED)
    """
).replace('STOPPED', '3')


def _insert(orig, star, out):
    """Run the insert in a child process group; return its code, or None at the time limit."""
    proc = subprocess.Popen(
        [sys.executable, '-c', CHILD, orig, star, out],
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


def test_long_output_path_gives_the_same_spectral_file_as_a_short_one(tmp_path):
    """A 190-byte folder is above the 132 of prep_spec and inside the limit from Cl_run_cdf:
    the call ends and writes the bytes of a short-path call, with the star in them. At 191
    bytes the call stops with the ValueError of the limit and writes nothing.
    """
    DownloadSpectralFiles('Oak')
    orig = str(spectral_file_dir('Oak') / 'Oak.sf')
    star = tmp_path / 'star.txt'
    rows = ''.join(f' {w:.7e} 1.0e+06\n' for w in [1e-10 * 1.5**i for i in range(80)])
    star.write_text(f'Star\n WAVELENGTH IRRADIANCE\n (m) (W/m3)\n*BEGIN_DATA\n{rows}*END\n ')

    short = tmp_path / 'short'
    short.mkdir()
    assert _insert(orig, str(star), str(short)) == 0

    def folder(length):
        path = str(tmp_path / 'long')
        while len(path) + 1 < length:
            path = os.path.join(path, 'd' * min(60, max(1, length - len(path) - 2)))
        os.makedirs(path)
        assert len(path) + 1 == length
        return path

    long = folder(190)
    assert _insert(orig, str(star), long) == 0
    for name in ('star.sf', 'star.sf_k'):
        with open(short / name, 'rb') as a, open(os.path.join(long, name), 'rb') as b:
            assert a.read() == b.read()
    text = (short / 'star.sf').read_text()
    assert '*BLOCK: TYPE =    2' in text and '*BLOCK: TYPE =    6' in text
    assert sorted(os.listdir(long)) == ['star.sf', 'star.sf_k']

    shutil.rmtree(tmp_path / 'long')
    over = folder(191)
    assert _insert(orig, str(star), over) == 3
    assert os.listdir(over) == []
