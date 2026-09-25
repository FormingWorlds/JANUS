"""Tests for src/janus/utils/data.py.

Exercises the FWL data-download plumbing with the OSF client mocked: folder
filtering and target-path construction in download_folder, the
skip-if-present logic of the stellar and spectral entry points, folder-list
selection per dataset name, and the unknown-name error contract.
See docs/How-to/test.md.
"""

import sys
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
    list is fetched, present means OSF is not contacted at all, so a cached
    tree needs no network. The module-level FWL_DATA_DIR constant is frozen at
    import, so it is patched directly.
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
    mock_osf.assert_not_called()


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
    with patch.object(jdata, 'OSF', MagicMock()) as mock_osf:
        jdata.DownloadSpectralFiles()
    assert called == []
    # A cached tree needs no network: OSF is not contacted.
    mock_osf.assert_not_called()

    # Empty name with one folder missing: exactly the missing one is fetched.
    (tmp_path / 'spectral_files' / 'Oak').rmdir()
    with patch.object(jdata, 'OSF', MagicMock()):
        jdata.DownloadSpectralFiles()
    assert called == [['Oak']]

    with pytest.raises(ValueError, match='Unrecognised folder name'):
        with patch.object(jdata, 'OSF', MagicMock()):
            jdata.DownloadSpectralFiles(fname='NotADataset')

    assert jdata.GetFWLData() == tmp_path.absolute()


def _cache_module():
    """Load tools/nightly_data_cache.py, skipping when its inputs are absent."""
    pytest.importorskip('mors')
    pytest.importorskip('fwl_io', minversion='26.8.31')
    import importlib.util

    path = Path(__file__).parents[2] / 'tools' / 'nightly_data_cache.py'
    spec = importlib.util.spec_from_file_location('nightly_data_cache', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_manifest(drc: Path, *, record: str, checksum: str) -> Path:
    """Write a one-dataset manifest and its registry, and return the manifest."""
    manifest = drc / 'mors_manifest.toml'
    manifest.write_text(
        '[star.tracks.baraffe_2015]\n'
        'name = "Baraffe tracks"\n'
        f'zenodo = "10.5281/zenodo.{record}"\n'
        'required_by = ["mors"]\n',
        encoding='utf-8',
    )
    registry = drc / 'star.tracks.baraffe_2015.registry.txt'
    registry.write_text(
        f'BHAC15-M0p010.txt md5:{checksum}\nBHAC15-M0p015.txt md5:{"b" * 32}\n',
        encoding='utf-8',
    )
    return manifest


def test_nightly_cache_key_moves_with_the_dataset_and_not_otherwise(monkeypatch, tmp_path):
    """The cache key tracks the record pin and the checksums, and nothing else.

    A key that stops moving with the data does not fail anything: the nightly
    stays green and silently refetches every run, which is the whole defect
    this key exists to prevent. Both directions are pinned here.
    """
    mod = _cache_module()
    import mors.data

    calls = []

    def key_for(record: str, checksum: str, name: str = 'Baraffe tracks') -> str:
        # A fresh directory per call, so no two cases can share a manifest and
        # agree for that reason rather than on their inputs.
        calls.append(record)
        drc = tmp_path / f'pins-{len(calls)}'
        drc.mkdir()
        manifest = _write_manifest(drc, record=record, checksum=checksum)
        if name != 'Baraffe tracks':
            manifest.write_text(
                manifest.read_text(encoding='utf-8').replace('Baraffe tracks', name),
                encoding='utf-8',
            )
        monkeypatch.setattr(mors.data, 'manifest_path', lambda: manifest, raising=True)
        return mod.resolve_key()

    baseline = key_for('15729114', 'a' * 32)

    # Re-pinning the record moves the data into a new r<record-id> directory.
    assert key_for('15729115', 'a' * 32) != baseline
    # Re-syncing the files changes the checksums without moving the directory.
    assert key_for('15729114', 'c' * 32) != baseline
    # A cosmetic manifest edit leaves the tree alone, so it must not cost a
    # full refetch of a dataset that has only one upstream mirror.
    assert key_for('15729114', 'a' * 32, name='BHAC15 tracks') == baseline
    # Same inputs, same key: a steady-state night has to hit its own entry.
    assert key_for('15729114', 'a' * 32) == baseline
    # An empty digest would equal the workflow's restore-key prefix, exact-hit
    # its own entry and freeze the tree with nothing reporting it.
    assert baseline.startswith(mod.KEY_PREFIX)
    assert len(baseline) > len(mod.KEY_PREFIX)


def test_nightly_workflow_derives_its_key_and_keeps_a_restore_prefix():
    """The workflow reads the resolved key and falls back to the shared prefix."""
    mod = _cache_module()
    workflow = (Path(__file__).parents[2] / '.github' / 'workflows' / 'nightly.yml').read_text(
        encoding='utf-8'
    )

    import re

    assert 'tools/nightly_data_cache.py key' in workflow
    assert 'key: ${{ steps.cachekey.outputs.key }}' in workflow
    # ANY literal, not just the one this replaced: actions/cache never rewrites
    # an entry whose key it hits, so a literal of any value freezes the tree.
    assert not re.search(rf'key:\s*{re.escape(mod.KEY_PREFIX)}\S', workflow)
    # Without the prefix a moved key starts cold and refetches every dataset,
    # including the OSF ones the key does not track.
    assert 'restore-keys:' in workflow
    assert f'\n            {mod.KEY_PREFIX}\n' in workflow
    # Fetch and check run on every night, in that order, before the tests.
    steps = workflow.split('      - name: ')
    for command in ('nightly_data_cache.py fetch', 'nightly_data_cache.py check'):
        (step,) = [s for s in steps if command in s]
        assert '\n        if:' not in step
    fetch = workflow.index('tools/nightly_data_cache.py fetch')
    assert (
        fetch
        < workflow.index('tools/nightly_data_cache.py check')
        < workflow.index('pytest -m')
    )


def test_cache_key_refuses_to_resolve_when_the_pins_are_missing(monkeypatch, tmp_path):
    """Resolving without a manifest or registry stops the job with a diagnostic."""
    mod = _cache_module()
    import mors.data

    missing = tmp_path / 'absent' / 'mors_manifest.toml'
    monkeypatch.setattr(mors.data, 'manifest_path', lambda: missing, raising=True)
    with pytest.raises(mod.ResolutionError, match='does not exist'):
        mod.resolve_key()
    assert mod.main(['key']) == 1

    drc = tmp_path / 'no-registry'
    drc.mkdir()
    manifest = _write_manifest(drc, record='15729114', checksum='a' * 32)
    (drc / 'star.tracks.baraffe_2015.registry.txt').unlink()
    monkeypatch.setattr(mors.data, 'manifest_path', lambda: manifest, raising=True)
    with pytest.raises(mod.ResolutionError, match='registry'):
        mod.resolve_key()

    empty = tmp_path / 'no-dataset'
    empty.mkdir()
    (empty / 'mors_manifest.toml').write_text('', encoding='utf-8')
    monkeypatch.setattr(
        mors.data, 'manifest_path', lambda: empty / 'mors_manifest.toml', raising=True
    )
    with pytest.raises(mod.ResolutionError, match='declares no dataset'):
        mod.resolve_key()


def test_restore_check_counts_registry_files_not_directories(monkeypatch, tmp_path):
    """The restore check counts intact registry files, not directory existence.

    A directory that exists but is short of its registry is exactly what a
    frozen cache looks like, so presence alone must not pass, and neither may
    a file whose contents differ from its pinned checksum.
    """
    import hashlib

    mod = _cache_module()
    monkeypatch.setattr(mod, 'OSF_DATA', ())
    import mors.data

    drc = tmp_path / 'pins'
    drc.mkdir()
    manifest = _write_manifest(drc, record='15729114', checksum=hashlib.md5(b'x').hexdigest())
    registry = drc / 'star.tracks.baraffe_2015.registry.txt'
    registry.write_text(
        registry.read_text(encoding='utf-8').replace('b' * 32, hashlib.md5(b'y').hexdigest()),
        encoding='utf-8',
    )
    monkeypatch.setattr(mors.data, 'manifest_path', lambda: manifest, raising=True)

    root = tmp_path / 'fwl_data'
    target = root / 'star' / 'tracks' / 'baraffe_2015' / 'r15729114'
    target.mkdir(parents=True)

    # An empty but existing directory is the frozen-cache case.
    assert mod.check_restored(root) == [('star/tracks/baraffe_2015/r15729114', 0, 2, 'intact')]
    assert mod.main(['check', '--data-root', str(root)]) == 1

    (target / 'BHAC15-M0p010.txt').write_text('x', encoding='utf-8')
    assert mod.check_restored(root) == [('star/tracks/baraffe_2015/r15729114', 1, 2, 'intact')]
    assert mod.main(['check', '--data-root', str(root)]) == 1

    # Present with the wrong contents is not intact.
    (target / 'BHAC15-M0p015.txt').write_text('x', encoding='utf-8')
    assert mod.check_restored(root) == [('star/tracks/baraffe_2015/r15729114', 1, 2, 'intact')]
    assert mod.main(['check', '--data-root', str(root)]) == 1

    (target / 'BHAC15-M0p015.txt').write_text('y', encoding='utf-8')
    assert mod.check_restored(root) == [('star/tracks/baraffe_2015/r15729114', 2, 2, 'intact')]
    assert mod.main(['check', '--data-root', str(root)]) == 0


def _write_archive_manifest(drc: Path, *, extract: bool) -> tuple[Path, bytes]:
    """Write a one-archive manifest and registry; return it and the tarball."""
    import hashlib
    import io
    import tarfile

    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode='w:gz') as tar:
        for name, data in (('fs255_grid/0p1.dat', b'track-a'), ('fs255_grid/0p2.dat', b'b')):
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    tarball = buf.getvalue()
    manifest = drc / 'mors_manifest.toml'
    manifest.write_text(
        '[star.tracks.spada_2013]\n'
        'name = "Spada tracks"\n'
        'zenodo = "10.5281/zenodo.15729101"\n'
        + ('extract = "tar"\n' if extract else '')
        + 'required_by = ["mors"]\n',
        encoding='utf-8',
    )
    (drc / 'star.tracks.spada_2013.registry.txt').write_text(
        f'fs255_grid.tar.gz md5:{hashlib.md5(tarball).hexdigest()}\n', encoding='utf-8'
    )
    return manifest, tarball


def test_fetch_extracts_an_archive_dataset_that_check_then_accepts(
    monkeypatch, tmp_path, capsys
):
    """An archive dataset is fetched, unpacked, and checked by its members.

    fwl-io discards the archive once it is unpacked, so a check that looks for
    the archive by name fails on a complete tree. Fetching runs through the
    fwl-io shared cache here, with the network switched off.
    """
    mod = _cache_module()
    monkeypatch.setattr(mod, 'OSF_DATA', ())
    import mors.data

    drc = tmp_path / 'pins'
    drc.mkdir()
    manifest, tarball = _write_archive_manifest(drc, extract=True)
    monkeypatch.setattr(mors.data, 'manifest_path', lambda: manifest, raising=True)
    rel_dir = 'star/tracks/spada_2013/r15729101'
    cache = tmp_path / 'cache'
    (cache / rel_dir).mkdir(parents=True)
    (cache / rel_dir / 'fs255_grid.tar.gz').write_bytes(tarball)
    monkeypatch.setenv('FWL_DATA_CACHE', str(cache))
    monkeypatch.setenv('FWL_IO_OFFLINE', '1')

    root = tmp_path / 'fwl_data'
    # Nothing fetched yet: one missing item, never zero of zero.
    assert mod.check_restored(root) == [(rel_dir, 0, 1, 'present')]
    assert mod.main(['check', '--data-root', str(root)]) == 1

    assert mod.main(['fetch', '--data-root', str(root)]) == 0
    assert (root / rel_dir / 'fs255_grid' / '0p1.dat').read_bytes() == b'track-a'
    assert not (root / rel_dir / 'fs255_grid.tar.gz').exists()
    assert mod.check_restored(root) == [(rel_dir, 2, 2, 'present')]
    assert mod.main(['check', '--data-root', str(root)]) == 0

    # A member lost after extraction makes the tree incomplete again.
    (root / rel_dir / 'fs255_grid' / '0p2.dat').unlink()
    assert mod.check_restored(root) == [(rel_dir, 1, 2, 'present')]
    capsys.readouterr()
    assert mod.main(['check', '--data-root', str(root)]) == 1
    assert 'missing files the registry pins or unpacked archive members' in (
        capsys.readouterr().err
    )


def test_fetch_and_check_fail_loudly(monkeypatch, tmp_path, capsys):
    """A failed download and an empty registry both stop the job with a diagnostic."""
    mod = _cache_module()
    import mors.data

    manifest, _ = _write_archive_manifest(tmp_path, extract=True)
    monkeypatch.setattr(mors.data, 'manifest_path', lambda: manifest, raising=True)
    monkeypatch.delenv('FWL_DATA_CACHE', raising=False)
    monkeypatch.setenv('FWL_IO_OFFLINE', '1')
    assert mod.main(['fetch', '--data-root', str(tmp_path)]) == 1
    assert 'error: reading or fetching the data failed' in capsys.readouterr().err

    # A registry with no entries would be 0 of 0; fwl-io refuses it outright.
    (tmp_path / 'star.tracks.spada_2013.registry.txt').write_text('', encoding='utf-8')
    assert mod.main(['check', '--data-root', str(tmp_path)]) == 1
    assert 'empty registry' in capsys.readouterr().err


def test_fetch_downloads_the_osf_data_the_tests_read_and_check_counts_it(
    monkeypatch, tmp_path, capsys
):
    """The OSF spectral file and stellar spectra the tests open are fetched into the data
    root and checked by sha256; a partial or corrupt folder, which the JANUS downloader
    would skip, is replaced, so no test has to download."""
    import hashlib

    mod = _cache_module()
    import mors.data

    pin = hashlib.sha256(b'x').hexdigest()
    monkeypatch.setattr(
        mod, 'OSF_DATA', tuple((d, dict.fromkeys(f, pin), c) for d, f, c in mod.OSF_DATA)
    )

    jdata = mod._janus_data()
    monkeypatch.setattr(mod, '_janus_data', lambda: jdata)
    manifest = _write_manifest(tmp_path, record='15729114', checksum='a' * 32)
    monkeypatch.setattr(mors.data, 'manifest_path', lambda: manifest, raising=True)
    monkeypatch.setattr(mod, '_fetchers', lambda root: [])
    root = tmp_path / 'fwl_data'
    assert mod.check_restored(root) == [
        ('spectral_files/Oak', 0, 2, 'intact'),
        ('stellar_spectra/Named', 0, 1, 'intact'),
    ]
    partial = root / 'stellar_spectra' / 'Named'
    partial.mkdir(parents=True)
    (partial / 'stale.txt').write_text('x')

    def _download(folder, files):
        def write():
            assert jdata.GetFWLData() == root
            if (root / folder).exists():  # as the JANUS downloader does
                return
            for name in files:
                (root / folder / name).parent.mkdir(parents=True, exist_ok=True)
                (root / folder / name).write_text('x')

        return write

    monkeypatch.setattr(
        jdata,
        'DownloadSpectralFiles',
        lambda fname: _download('spectral_files/Oak', ('318/Oak.sf', '318/Oak.sf_k'))(),
    )
    monkeypatch.setattr(
        jdata, 'DownloadStellarSpectra', _download('stellar_spectra/Named', ('sun.txt',))
    )
    assert mod.main(['fetch', '--data-root', str(root)]) == 0
    assert [r[1:3] for r in mod.check_restored(root)] == [(2, 2), (1, 1)]
    assert mod.main(['check', '--data-root', str(root)]) == 0
    (root / 'stellar_spectra' / 'Named' / 'sun.txt').write_text('truncated')
    assert mod.main(['check', '--data-root', str(root)]) == 1
    assert mod.main(['fetch', '--data-root', str(root)]) == 0
    assert mod.main(['check', '--data-root', str(root)]) == 0
    (root / 'stellar_spectra' / 'Named' / 'sun.txt').unlink()
    assert mod.main(['check', '--data-root', str(root)]) == 1


def test_key_moves_with_the_janus_downloader_source(monkeypatch, tmp_path):
    """The OSF project ids live in janus.utils.data, so its source is part of the key, as
    are the pinned sha256 of the OSF files."""
    mod = _cache_module()
    import mors.data

    manifest = _write_manifest(tmp_path, record='15729114', checksum='a' * 32)
    monkeypatch.setattr(mors.data, 'manifest_path', lambda: manifest, raising=True)
    before = mod.resolve_key()
    source = tmp_path / 'data.py'
    source.write_text("project_id = 'other'\n")
    monkeypatch.setattr(mod, 'JANUS_DATA_PY', source)
    moved = mod.resolve_key()
    assert moved != before
    assert mod.resolve_key() == moved
    pins = tuple((d, dict.fromkeys(f, '0' * 64), c) for d, f, c in mod.OSF_DATA)
    monkeypatch.setattr(mod, 'OSF_DATA', pins)
    assert mod.resolve_key() != moved


def test_key_and_osf_download_do_not_import_the_janus_package(monkeypatch, tmp_path):
    """The key and fetch steps run before SOCRATES is set up, and the janus package import
    checks for SOCRATES, so both read the downloader module by path."""
    mod = _cache_module()
    import mors.data

    manifest = _write_manifest(tmp_path, record='15729114', checksum='a' * 32)
    monkeypatch.setattr(mors.data, 'manifest_path', lambda: manifest, raising=True)
    for name in [m for m in sys.modules if m == 'janus' or m.startswith('janus.')]:
        monkeypatch.delitem(sys.modules, name)
    monkeypatch.setitem(sys.modules, 'janus', None)
    assert mod.resolve_key().startswith(mod.KEY_PREFIX)
    jdata = mod._janus_data()
    assert jdata.__file__ == str(mod.JANUS_DATA_PY)
    assert callable(jdata.DownloadSpectralFiles) and callable(jdata.DownloadStellarSpectra)


def test_key_moves_with_the_archive_kind(monkeypatch, tmp_path):
    """Unpacking a dataset changes its tree, so the key tracks the archive kind."""
    mod = _cache_module()
    import mors.data

    keys = []
    for extract in (True, False, True):
        drc = tmp_path / f'pins-{len(keys)}'
        drc.mkdir()
        manifest, _ = _write_archive_manifest(drc, extract=extract)
        monkeypatch.setattr(mors.data, 'manifest_path', lambda m=manifest: m, raising=True)
        keys.append(mod.resolve_key())
    assert keys[0] != keys[1]
    assert keys[0] == keys[2]


def test_key_moves_with_the_fwl_io_release_month(monkeypatch, tmp_path):
    """A new fwl-io month can change the tree layout, so it moves the key; a patch does not."""
    mod = _cache_module()
    import fwl_io
    import mors.data

    manifest = _write_manifest(tmp_path, record='15729114', checksum='a' * 32)
    monkeypatch.setattr(mors.data, 'manifest_path', lambda: manifest, raising=True)
    keys = {}
    for version in ('26.9.23', '26.9.30', '26.10.1'):
        monkeypatch.setattr(fwl_io, '__version__', version)
        keys[version] = mod.resolve_key()
    assert keys['26.9.23'] == keys['26.9.30']
    assert keys['26.9.23'] != keys['26.10.1']


def test_key_command_writes_the_output_line_the_workflow_reads(monkeypatch, tmp_path):
    """The key subcommand emits exactly the GITHUB_OUTPUT line the cache step consumes.

    This is the one contract wiring the script to the workflow. A malformed
    line leaves steps.cachekey.outputs.key empty, and an empty cache key means
    the restore-key prefix always wins, so the tree either freezes or refetches
    every night with nothing failing.
    """
    import re

    mod = _cache_module()
    out = tmp_path / 'gh_output'
    monkeypatch.setenv('GITHUB_OUTPUT', str(out))

    assert mod.main(['key']) == 0
    lines = out.read_text(encoding='utf-8').splitlines()
    assert len(lines) == 1
    assert re.fullmatch(rf'key={re.escape(mod.KEY_PREFIX)}[0-9a-f]{{64}}', lines[0])
    assert lines[0].split('=', 1)[1] == mod.resolve_key()

    # Appended, not truncated: a second call must not lose the first line.
    assert mod.main(['key']) == 0
    assert len(out.read_text(encoding='utf-8').splitlines()) == 2


def test_key_is_independent_of_registry_and_manifest_ordering(monkeypatch, tmp_path):
    """Reordering the registry lines or the manifest tables leaves the key alone.

    Neither order is part of the data, so a re-sync that shuffles them must
    not cost a full refetch of a dataset with one upstream mirror.
    """
    mod = _cache_module()
    import mors.data

    def key_for(lines: str, tag: str) -> str:
        drc = tmp_path / tag
        drc.mkdir()
        manifest = drc / 'mors_manifest.toml'
        manifest.write_text(
            '[star.tracks.baraffe_2015]\n'
            'name = "Baraffe tracks"\n'
            'zenodo = "10.5281/zenodo.15729114"\n'
            'required_by = ["mors"]\n',
            encoding='utf-8',
        )
        (drc / 'star.tracks.baraffe_2015.registry.txt').write_text(lines, encoding='utf-8')
        monkeypatch.setattr(mors.data, 'manifest_path', lambda: manifest, raising=True)
        return mod.resolve_key()

    first = f'BHAC15-M0p010.txt md5:{"a" * 32}\nBHAC15-M0p015.txt md5:{"b" * 32}\n'
    reversed_lines = f'BHAC15-M0p015.txt md5:{"b" * 32}\nBHAC15-M0p010.txt md5:{"a" * 32}\n'
    assert key_for(first, 'forward') == key_for(reversed_lines, 'reverse')
    # Swapping which file carries which checksum IS a data change and must move it.
    swapped = f'BHAC15-M0p010.txt md5:{"b" * 32}\nBHAC15-M0p015.txt md5:{"a" * 32}\n'
    assert key_for(swapped, 'swapped') != key_for(first, 'forward-again')

    # The order of the tables in the manifest is not part of the data either.
    keys = []
    for tag in ('baraffe-first', 'spada-first'):
        drc = tmp_path / tag
        drc.mkdir()
        manifest = _write_manifest(drc, record='15729114', checksum='a' * 32)
        baraffe = manifest.read_text(encoding='utf-8')
        spada = _write_archive_manifest(drc, extract=True)[0].read_text(encoding='utf-8')
        tables = (baraffe, spada) if tag == 'baraffe-first' else (spada, baraffe)
        manifest.write_text('\n'.join(tables), encoding='utf-8')
        monkeypatch.setattr(mors.data, 'manifest_path', lambda m=manifest: m, raising=True)
        keys.append(mod.resolve_key())
    assert keys[0] == keys[1]
    assert keys[0] != key_for(first, 'baraffe-only')


def test_check_refuses_to_run_without_a_data_root(monkeypatch, capsys):
    """With no root given, check reports it rather than inspecting the working directory.

    Path('') is Path('.'), so a guard applied after the Path is built cannot
    fire and would silently count files in whatever directory the job sits in.
    """
    mod = _cache_module()
    monkeypatch.delenv('FWL_DATA', raising=False)

    with pytest.raises(mod.ResolutionError, match='no data root'):
        mod._cmd_check(SimpleNamespace(data_root=None))
    assert mod.main(['check']) == 1
    # fetch must not download into the working directory either.
    assert mod.main(['fetch']) == 1
    assert 'no data root' in capsys.readouterr().err

    # An empty string is the same case and must not fall through to the CWD.
    monkeypatch.setenv('FWL_DATA', '')
    with pytest.raises(mod.ResolutionError, match='no data root'):
        mod._cmd_check(SimpleNamespace(data_root=None))
