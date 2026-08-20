"""Tests for the OSF retry helper in src/janus/utils/data.py.

Exercises `_osf_retry` and its use inside `download_folder` and
`DownloadStellarSpectra` with the OSF client mocked: retrying a transient
listing or write failure, giving up once the retry budget is spent, retrying
a listing that fails partway through iteration, and not retrying a
non-transient error. See docs/How-to/test.md.
"""

import re
import time

import pytest
import requests

import janus.utils.data as data_module
from janus.utils.data import OSF_RETRY_ATTEMPTS, _osf_retry, download_folder

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

_OSF_502 = 'Response has status code 502 not (200,)'
_OSF_404 = 'Response has status code 404 not (200,)'


class _FakeFile:
    """File stub whose write_to() fails a fixed number of times before succeeding."""

    def __init__(self, path, fail_times=0, content=b'file contents'):
        self.path = path
        self._fail_times = fail_times
        self._content = content
        self.calls = 0

    def write_to(self, fileobj):
        self.calls += 1
        if self.calls <= self._fail_times:
            raise RuntimeError(_OSF_502)
        fileobj.write(self._content)


class _FlakyStorage:
    """Storage stub whose `.files` listing fails a fixed number of times before succeeding."""

    def __init__(self, files, fail_times=0):
        self._files = files
        self._fail_times = fail_times
        self.calls = 0

    @property
    def files(self):
        self.calls += 1
        if self.calls <= self._fail_times:
            raise RuntimeError(_OSF_502)
        return iter(self._files)


class _PartialListingStorage:
    """Storage stub whose `.files` is a fresh generator each access, like the
    real (lazy, paginated) osfclient `Storage.files`. It fails partway
    through iteration, not at access time, for a fixed number of accesses.
    """

    def __init__(self, files, fail_after, fail_times=0):
        self._files = files
        self._fail_after = fail_after
        self._fail_times = fail_times
        self.attempts = 0

    @property
    def files(self):
        self.attempts += 1
        will_fail = self.attempts <= self._fail_times

        def _gen():
            for i, f in enumerate(self._files):
                if will_fail and i == self._fail_after:
                    raise RuntimeError(_OSF_502)
                yield f

        return _gen()


@pytest.fixture(autouse=True)
def _no_sleep(monkeypatch):
    # The retry helper sleeps between attempts; skip the real delay in tests.
    monkeypatch.setattr(time, 'sleep', lambda seconds: None)


def test_download_folder_retries_transient_listing_failure(tmp_path):
    """A transient failure listing storage.files is retried until it succeeds."""
    file = _FakeFile('/Oak/spectrum.sf')
    storage = _FlakyStorage([file], fail_times=OSF_RETRY_ATTEMPTS - 1)

    download_folder(storage=storage, folders=['Oak'], data_dir=tmp_path)

    assert storage.calls == OSF_RETRY_ATTEMPTS
    assert (tmp_path / 'Oak' / 'spectrum.sf').read_bytes() == b'file contents'


def test_download_folder_retries_transient_write_failure(tmp_path):
    """A transient failure writing a single file is retried until it succeeds."""
    file = _FakeFile('/Oak/spectrum.sf', fail_times=OSF_RETRY_ATTEMPTS - 1)
    storage = _FlakyStorage([file])

    download_folder(storage=storage, folders=['Oak'], data_dir=tmp_path)

    assert file.calls == OSF_RETRY_ATTEMPTS
    assert (tmp_path / 'Oak' / 'spectrum.sf').read_bytes() == b'file contents'


def test_download_folder_gives_up_after_retry_budget(tmp_path):
    """A failure that outlasts the retry budget propagates instead of looping forever."""
    storage = _FlakyStorage([], fail_times=OSF_RETRY_ATTEMPTS + 1)

    with pytest.raises(RuntimeError, match=re.escape(_OSF_502)):
        download_folder(storage=storage, folders=['Oak'], data_dir=tmp_path)

    assert storage.calls == OSF_RETRY_ATTEMPTS


def test_download_folder_retries_listing_failure_mid_iteration(tmp_path):
    """A listing that fails partway through iteration is retried from scratch.

    The real osfclient Storage.files is a lazy, paginated generator that can
    fail partway through, not a property that fails before yielding
    anything. Confirm the retry-from-scratch on `list(storage.files)` still
    produces the complete, correct file set once a later attempt succeeds.
    """
    files = [_FakeFile('/Oak/a.sf'), _FakeFile('/Oak/b.sf'), _FakeFile('/Oak/c.sf')]
    storage = _PartialListingStorage(files, fail_after=1, fail_times=OSF_RETRY_ATTEMPTS - 1)

    download_folder(storage=storage, folders=['Oak'], data_dir=tmp_path)

    assert storage.attempts == OSF_RETRY_ATTEMPTS
    for name in ('a.sf', 'b.sf', 'c.sf'):
        assert (tmp_path / 'Oak' / name).read_bytes() == b'file contents'


def test_osf_retry_does_not_retry_non_transient_status():
    """A 404 is not a transient status, so `_osf_retry` raises after one attempt."""
    calls = []

    def _raise_404():
        calls.append(1)
        raise RuntimeError(_OSF_404)

    with pytest.raises(RuntimeError, match=re.escape(_OSF_404)):
        _osf_retry(_raise_404)

    assert len(calls) == 1


def test_osf_retry_retries_network_error():
    """A connection error is retried until the call succeeds within budget."""
    calls = []

    def _flaky():
        calls.append(1)
        if len(calls) < OSF_RETRY_ATTEMPTS:
            raise requests.exceptions.ConnectionError('connection reset')
        return 'ok'

    result = _osf_retry(_flaky)

    assert result == 'ok'
    assert len(calls) == OSF_RETRY_ATTEMPTS


def test_download_stellar_spectra_retries_project_and_storage_lookup(tmp_path, monkeypatch):
    """Transient failures resolving the OSF project or its storage are retried."""
    project_calls = []
    storage_calls = []

    class _FakeProject:
        def storage(self, name):
            storage_calls.append(1)
            if len(storage_calls) < OSF_RETRY_ATTEMPTS:
                raise RuntimeError(_OSF_502)
            return object()

    class _FakeOSF:
        def project(self, project_id):
            project_calls.append(1)
            if len(project_calls) < OSF_RETRY_ATTEMPTS:
                raise RuntimeError(_OSF_502)
            return _FakeProject()

    monkeypatch.setattr(data_module, 'OSF', lambda: _FakeOSF())
    monkeypatch.setattr(data_module, 'GetFWLData', lambda: tmp_path)
    monkeypatch.setattr(data_module, 'download_folder', lambda **kwargs: None)

    data_module.DownloadStellarSpectra()

    assert len(project_calls) == OSF_RETRY_ATTEMPTS
    assert len(storage_calls) == OSF_RETRY_ATTEMPTS
