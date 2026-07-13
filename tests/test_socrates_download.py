"""Tests for src/janus/socrates.py.

The file carries a unique basename because tests/utils/test_socrates.py
mirrors the SOCRATES driver module of the same name.

Exercises the SOCRATES download helper with the network and build mocked:
chunked download to disk, HTTP error propagation, archive extraction,
permission setting, build invocation, and symlink management for repeated
installs. See docs/How-to/test.md.
"""

import zipfile
from unittest.mock import MagicMock, patch

import pytest
import requests

import janus.socrates as jsoc

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def _mock_response(chunks, status_ok=True):
    """A requests.get context manager yielding the given byte chunks."""
    resp = MagicMock()
    resp.__enter__.return_value = resp
    resp.__exit__.return_value = False
    resp.headers = {'Content-Length': str(sum(len(c) for c in chunks))}
    resp.iter_content.return_value = iter(chunks)
    if not status_ok:
        resp.raise_for_status.side_effect = requests.HTTPError('404 Client Error')
    return resp


def test_download_writes_chunks_and_propagates_http_error(tmp_path, monkeypatch):
    """The chunked download reassembles the payload; an HTTP error aborts
    before any write.

    The reassembled file must equal the concatenated chunks byte for byte;
    on a failing status nothing may be written because raise_for_status
    fires before the output file is opened.
    """
    monkeypatch.chdir(tmp_path)
    chunks = [b'fortran ', b'radiative ', b'transfer']

    with patch.object(jsoc.requests, 'get', return_value=_mock_response(chunks)):
        jsoc._download(url='https://example.invalid/s.zip', filename='s.zip')
    assert (tmp_path / 's.zip').read_bytes() == b''.join(chunks)

    bad = _mock_response([b'x'], status_ok=False)
    with patch.object(jsoc.requests, 'get', return_value=bad):
        with pytest.raises(requests.HTTPError):
            jsoc._download(url='https://example.invalid/bad.zip', filename='bad.zip')
    assert not (tmp_path / 'bad.zip').exists()


def _make_socrates_zip(path, subdir):
    """A minimal SOCRATES source archive with the files the installer touches."""
    with zipfile.ZipFile(path, 'w') as zf:
        zf.writestr(f'{subdir}/make/mkdep', '#!/bin/sh\n')
        zf.writestr(f'{subdir}/configure', 'echo configure\n')
        zf.writestr(f'{subdir}/build_code', 'echo build\n')
        zf.writestr(f'{subdir}/version', '2311.4\n')


@pytest.mark.parametrize('ref', ['main', 'abc123'], ids=['default-branch', 'pinned-hash'])
def test_download_socrates_extracts_builds_and_links(tmp_path, monkeypatch, ref):
    """The installer extracts the archive, marks mkdep executable, runs the
    two build steps in the tree, and points the SOCRATES symlink at it.

    The requested git ref selects both the download URL and the extracted
    subdirectory name; asserting both catches a ref mixed into only one of
    the two places. Re-running must replace the existing symlink (the
    missing_ok unlink branch) instead of failing.
    """
    monkeypatch.chdir(tmp_path)
    data_dir = tmp_path / 'sdata'
    monkeypatch.setattr(jsoc, 'SOCRATES_DATA_DIR', data_dir, raising=True)

    urls = []

    def fake_download(*, url, filename):
        urls.append(url)
        _make_socrates_zip(filename, f'SOCRATES-{ref}')

    run_calls = []
    monkeypatch.setattr(jsoc, '_download', fake_download)
    monkeypatch.setattr(
        jsoc.sp, 'run', lambda cmd, cwd=None, **kw: run_calls.append((cmd, cwd))
    )

    jsoc.download_socrates(ref=ref)
    # Second run exercises symlink replacement.
    jsoc.download_socrates(ref=ref)

    expected_path = 'refs/heads/main' if ref == 'main' else ref
    assert all(u.endswith(f'/{expected_path}.zip') for u in urls)
    target = data_dir / f'SOCRATES-{ref}'
    assert (target / 'configure').exists()
    # mkdep must carry the executable bits set by _set_permissions.
    assert (target / 'make' / 'mkdep').stat().st_mode == 33261
    # Both build steps run inside the extracted tree, in order.
    assert [c[0] for c in run_calls[:2]] == [['bash', 'configure'], ['bash', 'build_code']]
    assert all(cwd == target for _, cwd in run_calls)
    link = data_dir / 'SOCRATES'
    assert link.is_symlink()
    assert link.resolve() == target.resolve()
