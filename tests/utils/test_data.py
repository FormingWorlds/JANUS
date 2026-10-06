"""Tests for src/janus/utils/data.py and tools/nightly_data_cache.py.

Exercises the fwl-io data plumbing without the network: the spectral-file and
stellar-spectra keys and directories, the downloads each entry point asks
fwl-io for, the unknown-name error contract, and the nightly cache tool.
See docs/How-to/test.md.
"""

from pathlib import Path
from types import SimpleNamespace

import pytest

from janus.utils import data as jdata

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


class _Fetcher:
    """Records the datasets the data module asks fwl-io to fetch."""

    fetched = []

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.target_dir = kwargs['data_root'] / kwargs['subdir']

    def fetch_all(self):
        self.fetched.append(self.kwargs['subdir'])


@pytest.fixture
def fetches(monkeypatch, tmp_path):
    """Route janus data fetches to a recorder below tmp_path and return the log."""
    import fwl_io

    _Fetcher.fetched = []
    monkeypatch.setattr(jdata, 'FWL_DATA_DIR', tmp_path)
    monkeypatch.setattr(fwl_io, 'create_fetcher', lambda **kw: _Fetcher(**kw))
    return _Fetcher.fetched


@pytest.mark.parametrize(
    ('fname', 'nband', 'expected'),
    [
        ('Dayspring', 48, ['atmos_clim/spectral_files/dayspring/48']),
        ('Frostflow', 4096, ['atmos_clim/spectral_files/frostflow/4096']),
        ('Oak', 256, ['atmos_clim/spectral_files/oak/318']),
        (
            '',
            16,
            [
                'atmos_clim/spectral_files/dayspring/256',
                'atmos_clim/spectral_files/frostflow/256',
                'atmos_clim/spectral_files/oak/318',
            ],
        ),
    ],
)
def test_spectral_download_asks_fwl_io_for_the_named_dataset(fetches, fname, nband, expected):
    """A group with several band counts fetches the one asked for; Oak has only 318.

    The empty name fetches the default list whatever nband says, and nothing else.
    """
    jdata.DownloadSpectralFiles(fname=fname, nband=nband)
    assert fetches == expected
    assert len(set(fetches)) == len(fetches)


def test_unknown_spectral_file_raises_before_any_fetch(fetches):
    """A group or band count the manifest does not declare raises and fetches nothing."""
    for fname, nband in (('Dayspring', 100), ('NotADataset', 256)):
        with pytest.raises(ValueError, match='No spectral file'):
            jdata.DownloadSpectralFiles(fname=fname, nband=nband)
    assert fetches == []


def test_stellar_spectra_and_directories_resolve_through_fwl_io(fetches, tmp_path):
    """The named spectra are fetched by key, and the helpers return the dataset
    directories that hold sun.txt and Oak.sf."""
    jdata.DownloadStellarSpectra()
    assert fetches == ['star/spectra/named']
    assert jdata.stellar_spectra_dir() == tmp_path / 'star/spectra/named'
    assert jdata.spectral_file_dir('Oak') == tmp_path / 'atmos_clim/spectral_files/oak/318'
    assert jdata.GetFWLData() == tmp_path.absolute()


def test_fetches_carry_the_manifest_pins_of_the_dataset(fetches, monkeypatch):
    """Each fetch gets the Zenodo record, the DataverseNL mirror and the registry of
    its own manifest entry."""
    import fwl_io

    seen = []
    monkeypatch.setattr(
        fwl_io, 'create_fetcher', lambda **kw: seen.append(kw) or _Fetcher(**kw)
    )
    jdata.DownloadSpectralFiles('Oak')
    ds = jdata._shared_datasets()['atmos_clim.spectral_files.oak.318']
    assert seen[0]['zenodo'] == ds.zenodo == '10.5281/zenodo.15743843'
    assert seen[0]['dataverse'] == ds.dataverse and seen[0]['dataverse'].startswith('10.34894/')
    assert seen[0]['registry'] == ds.registry() and 'Oak.sf' in seen[0]['registry']


def test_directories_are_the_version_directories_fwl_io_fills(monkeypatch, tmp_path):
    """With the real fwl-io fetcher, the helpers return the r<record-id> directory."""
    monkeypatch.setattr(jdata, 'FWL_DATA_DIR', tmp_path)
    oak = jdata.spectral_file_dir('Oak')
    named = jdata.stellar_spectra_dir()
    assert oak == tmp_path / 'atmos_clim/spectral_files/oak/318/r15743843'
    assert named.parent == tmp_path / 'star/spectra/named' and named.name.startswith('r')
    assert not oak.exists()


def test_band_count_and_group_errors(fetches, caplog):
    """A single-band group says when it overrides the band count; an unknown group or
    a multi-band group without a count raises before any fetch."""
    import logging

    with caplog.at_level(logging.INFO, logger='fwl.janus.utils.data'):
        assert jdata.spectral_file_key('Oak', 4096).endswith('oak.318')
    assert 'one band count; using 318' in caplog.text
    with pytest.raises(ValueError, match="No spectral file group 'Mallard'"):
        jdata.DownloadSpectralFiles('Mallard')
    with pytest.raises(ValueError, match=r"band counts declared for Dayspring: \['16'"):
        jdata.spectral_file_dir('Dayspring')
    assert fetches == []


def _cache_module():
    """Load tools/nightly_data_cache.py, skipping when its inputs are absent."""
    pytest.importorskip('mors')
    pytest.importorskip('fwl_io', minversion='26.10.6')
    import importlib.util

    path = Path(__file__).parents[2] / 'tools' / 'nightly_data_cache.py'
    spec = importlib.util.spec_from_file_location('nightly_data_cache', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_nightly_shared_keys_are_the_datasets_the_test_helpers_read(tmp_path):
    """SHARED_KEYS is the Oak spectral file and the named spectra, and the tool builds
    their fetchers at the version directories, so the key and the cache cover them."""
    mod = _cache_module()
    assert set(mod.SHARED_KEYS) == {jdata.spectral_file_key('Oak'), jdata.STELLAR_SPECTRA_NAMED}
    rel = [f.rel_dir for f in mod._fetchers(tmp_path)]
    assert 'atmos_clim/spectral_files/oak/318/r15743843' in rel
    assert any(r.startswith('star/spectra/named/r') for r in rel)


def test_nightly_key_tracks_the_shared_keys_and_refuses_an_unknown_one(monkeypatch, capsys):
    """Dropping a shared dataset moves the key; a key the shared manifest lacks is refused."""
    mod = _cache_module()
    both = mod.resolve_key()
    monkeypatch.setattr(mod, 'SHARED_KEYS', ('star.spectra.named',))
    assert mod.resolve_key() != both
    monkeypatch.setattr(mod, 'SHARED_KEYS', ('atmos_clim.spectral_files.none.1',))
    with pytest.raises(mod.ResolutionError, match='shared manifest declares no'):
        mod.resolve_key()


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
    # Without the prefix a moved key starts cold and refetches every dataset.
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
    monkeypatch.setattr(mod, 'SHARED_KEYS', ())
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
    monkeypatch.setattr(mod, 'SHARED_KEYS', ())
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
