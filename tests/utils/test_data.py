"""Tests for src/janus/utils/data.py.

Exercises the FWL data-download plumbing with the OSF client mocked: folder
filtering and target-path construction in download_folder, the
skip-if-present logic of the stellar and spectral entry points, folder-list
selection per dataset name, and the unknown-name error contract.
See docs/How-to/test.md.
"""

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from janus.utils import data as jdata

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


class _FakeFile:
    """OSF file stub that records what gets written where."""

    def __init__(self, path, payload=b'flux-table'):
        self.path = path
        self.payload = payload

    def write_to(self, handle):
        handle.write(self.payload)


def _storage(files):
    return SimpleNamespace(files=files)


def test_download_folder_filters_and_writes_targets(tmp_path):
    """Only files under the requested folders are written, at mirrored paths.

    The OSF paths carry a leading slash and nested directories; the target
    must reproduce the nesting under data_dir and the payload must survive
    the write. A file outside the requested folders must not appear.
    """
    files = [
        _FakeFile('/Named/sun.txt', b'solar spectrum'),
        _FakeFile('/Named/sub/hd97658.txt', b'k-dwarf spectrum'),
        _FakeFile('/Other/ignore.txt', b'unrelated'),
    ]
    jdata.download_folder(storage=_storage(files), folders=['Named'], data_dir=tmp_path)

    assert (tmp_path / 'Named' / 'sun.txt').read_bytes() == b'solar spectrum'
    assert (tmp_path / 'Named' / 'sub' / 'hd97658.txt').read_bytes() == b'k-dwarf spectrum'
    assert not (tmp_path / 'Other').exists()


def test_stellar_spectra_download_skips_when_present(tmp_path, monkeypatch):
    """DownloadStellarSpectra downloads once and skips when data exist.

    The presence check is on the Named subfolder: absent means the folder
    list is fetched, present means no OSF file is touched. The module-level
    FWL_DATA_DIR constant is frozen at import, so it is patched directly.
    """
    monkeypatch.setattr(jdata, 'FWL_DATA_DIR', tmp_path, raising=True)
    files = [_FakeFile('/Named/sun.txt')]

    with patch.object(jdata, 'OSF') as mock_osf:
        mock_osf.return_value.project.return_value.storage.return_value = _storage(files)
        jdata.DownloadStellarSpectra()
    assert (tmp_path / 'stellar_spectra' / 'Named' / 'sun.txt').exists()

    # Second call: folder exists, so no file may be rewritten.
    marker = tmp_path / 'stellar_spectra' / 'Named' / 'sun.txt'
    marker.write_bytes(b'unchanged')
    with patch.object(jdata, 'OSF') as mock_osf:
        mock_osf.return_value.project.return_value.storage.return_value = _storage(
            [_FakeFile('/Named/sun.txt', b'would overwrite')]
        )
        jdata.DownloadStellarSpectra()
    assert marker.read_bytes() == b'unchanged'


@pytest.mark.parametrize(
    'fname,nband,expected',
    [
        ('Dayspring', 4096, ['Dayspring/4096']),
        ('Oak', 256, ['Oak']),
    ],
    ids=['banded-dataset', 'flat-dataset'],
)
def test_spectral_download_folder_selection(tmp_path, monkeypatch, fname, nband, expected):
    """The folder list follows the dataset naming scheme.

    Banded datasets (Dayspring) select a resolution subfolder from nband;
    flat datasets (Oak) ignore nband entirely. The selected folders arrive
    verbatim at download_folder, which pins the dispatch.
    """
    monkeypatch.setattr(jdata, 'FWL_DATA_DIR', tmp_path, raising=True)
    seen = {}

    def fake_download_folder(*, storage, folders, data_dir):
        seen['folders'] = list(folders)
        seen['data_dir'] = Path(data_dir)

    monkeypatch.setattr(jdata, 'download_folder', fake_download_folder)
    with patch.object(jdata, 'OSF', MagicMock()):
        jdata.DownloadSpectralFiles(fname=fname, nband=nband)

    assert seen['folders'] == expected
    assert seen['data_dir'] == tmp_path / 'spectral_files'


def test_spectral_download_basic_list_skip_and_error(tmp_path, monkeypatch):
    """No name selects the basic list, present folders are skipped, and an
    unknown name raises.

    With every basic folder already on disk the download must not run at
    all; the ValueError branch is the documented contract for a typo in the
    dataset name.
    """
    monkeypatch.setattr(jdata, 'FWL_DATA_DIR', tmp_path, raising=True)
    called = []
    monkeypatch.setattr(jdata, 'download_folder', lambda **kw: called.append(kw['folders']))

    # All basic-list folders present: nothing to download.
    for folder in jdata.basic_list:
        (tmp_path / 'spectral_files' / folder).mkdir(parents=True)
    with patch.object(jdata, 'OSF', MagicMock()):
        jdata.DownloadSpectralFiles()
    assert called == []

    # Empty name with one folder missing: exactly the missing one is fetched.
    (tmp_path / 'spectral_files' / 'Oak').rmdir()
    with patch.object(jdata, 'OSF', MagicMock()):
        jdata.DownloadSpectralFiles()
    assert called == [['Oak']]

    with pytest.raises(ValueError, match='Unrecognised folder name'):
        with patch.object(jdata, 'OSF', MagicMock()):
            jdata.DownloadSpectralFiles(fname='NotADataset')

    assert jdata.GetFWLData() == tmp_path.absolute()
