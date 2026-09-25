#!/usr/bin/env python3
"""Cache key and restore check for the FWL data tree the nightly caches.

Three subcommands, all used by ``.github/workflows/nightly.yml``::

    python tools/nightly_data_cache.py key
    python tools/nightly_data_cache.py fetch
    python tools/nightly_data_cache.py check

``key`` prints ``key=<value>`` for ``GITHUB_OUTPUT``. The value carries a
digest of the track-data layout JANUS resolves through its ``fwl-mors``
dependency: for every dataset the installed manifest declares, the
directory fwl-io places it in and the per-file checksums the registry
pins, plus the archive kind of a dataset shipped as one archive and the
fwl-io release (year.month) that lays out the tree. It also carries the
source of ``janus.utils.data``, which pins the OSF projects of the
spectral file and stellar spectra the tests read. It therefore moves
when the data or its layout moves and stays put otherwise.

Every declared dataset counts, not only the one JANUS fetches today. A
dataset the manifest gains later costs at most one extra refetch, where
narrowing the digest to a named subset would leave a dataset JANUS starts
consuming untracked, which is the failure worth avoiding.

``fetch`` downloads every dataset the key covers that is not already in
place, so a tree saved under the key holds all of it and no test
downloads. ``check`` verifies that tree without the network: every
registry file present with its pinned checksum, for an archive dataset
every member its extraction recorded, and for the OSF data the files the
tests open, the last two by name.

Both subcommands fail with a diagnostic rather than degrade: an empty or
partial digest would collide with the workflow's restore-key prefix and
freeze the cached tree with nothing reporting it.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import shutil
import sys
import tempfile
from pathlib import Path

KEY_PREFIX = 'fwl-data-'

# JANUS data on OSF that tests/helpers reads: folder, files the tests open, download call.
OSF_DATA = (
    (
        'spectral_files/Oak',
        ('318/Oak.sf', '318/Oak.sf_k'),
        lambda j: j.DownloadSpectralFiles('Oak'),
    ),
    ('stellar_spectra/Named', ('sun.txt',), lambda j: j.DownloadStellarSpectra()),
)


class ResolutionError(RuntimeError):
    """The datasets that decide the cached layout could not be resolved."""


def _manifest_path() -> Path:
    """Return the dataset manifest shipped inside the installed fwl-mors.

    Returns
    -------
    Path
        Absolute path of ``mors_manifest.toml``.

    Raises
    ------
    ResolutionError
        When fwl-mors is absent, exposes no manifest, or the manifest
        file it names does not exist.
    """
    try:
        import mors.data
    except ImportError as exc:
        raise ResolutionError(
            'fwl-mors is not importable, so the Baraffe layout this key tracks '
            'cannot be resolved. Install JANUS with its dependencies before '
            'resolving the cache key.'
        ) from exc

    if not hasattr(mors.data, 'manifest_path'):
        raise ResolutionError(
            'the installed fwl-mors exposes no mors.data.manifest_path(), so the '
            'dataset pins are no longer where this script looks for them. Point it '
            'at whatever now declares the Baraffe record and checksums.'
        )

    path = Path(mors.data.manifest_path())
    if not path.is_file():
        raise ResolutionError(f'the dataset manifest fwl-mors names does not exist: {path}')
    return path


def _fetchers(data_root: Path) -> list:
    """Build one fwl-io fetcher per dataset the manifest declares.

    Parameters
    ----------
    data_root : Path
        Root the fetchers resolve their target directories below.

    Returns
    -------
    list
        Fetchers, ordered by dataset key.

    Raises
    ------
    ResolutionError
        When the manifest declares no dataset, or a declared registry
        file is missing.
    """
    from fwl_io import create_fetcher, load_manifest

    manifest = _manifest_path()
    datasets = sorted(load_manifest(manifest), key=lambda ds: ds.key)
    if not datasets:
        raise ResolutionError(
            f'{manifest} declares no dataset, so this key would track nothing.'
        )

    built = []
    for ds in datasets:
        if not Path(ds.registry_path).is_file():
            raise ResolutionError(
                f'dataset {ds.key!r} declares a registry at {ds.registry_path}, '
                'which does not exist. The checksums are half of what this key '
                'tracks, so resolving it without them would freeze the cache.'
            )
        built.append(
            create_fetcher(
                subdir=ds.subdir,
                zenodo=ds.zenodo,
                dataverse=ds.dataverse,
                registry=ds.registry_path,
                data_root=data_root,
                extract=ds.extract,
            )
        )
    return built


def resolve_key(data_root: Path | None = None) -> str:
    """Return the cache key for the FWL data tree.

    Parameters
    ----------
    data_root : Path | None
        Root passed to the fetchers. Only relative layout is hashed, so
        any writable directory gives the same key; a temporary one is
        used when this is None.

    Returns
    -------
    str
        ``fwl-data-<sha256>``.

    Raises
    ------
    ResolutionError
        When the datasets cannot be resolved, or the digest comes out
        empty and would collapse the key onto the restore-key prefix.
    """
    if data_root is None:
        with tempfile.TemporaryDirectory() as tmp:
            return resolve_key(Path(tmp))

    material = []
    for f in _fetchers(data_root):
        material.append(f'dir\t{f.rel_dir}')
        if f.extract:
            material.append(f'extract\t{f.extract}')
        for name in sorted(f.registry):
            material.append(f'file\t{name}\t{f.registry[name]}')

    if not material:
        raise ResolutionError(
            'no dataset directory or checksum was resolved, so the key would be '
            'the bare restore-key prefix and the cached tree could never be '
            'rewritten.'
        )
    from fwl_io import __version__

    material.append('fwl-io\t' + '.'.join(__version__.split('.')[:2]))
    source = Path(importlib.util.find_spec('janus.utils.data').origin).read_bytes()
    material.append('janus-osf\t' + hashlib.sha256(source).hexdigest())
    material += [f'osf\t{folder}\t{" ".join(files)}' for folder, files, _ in OSF_DATA]
    digest = hashlib.sha256('\n'.join(material).encode('utf-8')).hexdigest()
    return f'{KEY_PREFIX}{digest}'


def check_restored(data_root: Path) -> list[tuple[str, int, int, str]]:
    """Report how much of each dataset is present below ``data_root``.

    Parameters
    ----------
    data_root : Path
        Root of the restored FWL data tree.

    Returns
    -------
    list of tuple
        One ``(rel_dir, found, expected, state)`` per dataset. ``found``
        counts the files fwl-io reports sound, and ``state`` names what that
        means: ``intact`` for a plain dataset, checked by checksum, and
        ``present`` for the members of an archive dataset and for the OSF
        files the tests open, checked by name. An archive dataset with no
        extracted tree counts as its one archive, missing.

    Raises
    ------
    ResolutionError
        When the datasets cannot be resolved.
    """
    from fwl_io import check_dataset

    report = []
    for f in _fetchers(data_root):
        ds = check_dataset(f)
        state = 'intact' if ds.verifiable else 'present'
        report.append((f.rel_dir, sum(not c.faulty for c in ds.files), len(ds.files), state))
    for folder, files, _ in OSF_DATA:
        found = sum((data_root / folder / name).is_file() for name in files)
        report.append((folder, found, len(files), 'present'))
    return report


def fetch_all(data_root: Path) -> None:
    """Fetch every dataset the key covers into ``data_root``.

    Parameters
    ----------
    data_root : Path
        Root of the FWL data tree the nightly caches.

    Raises
    ------
    ResolutionError
        When the datasets cannot be resolved.
    """
    for f in _fetchers(data_root):
        f.fetch_all()
        print(f'{f.rel_dir}: fetched', file=sys.stderr)

    import janus.utils.data as jdata

    jdata.FWL_DATA_DIR = Path(data_root)
    for folder, files, download in OSF_DATA:
        path = Path(data_root) / folder
        # The JANUS downloader skips a folder that exists, so a partial one is removed.
        if path.exists() and not all((path / name).is_file() for name in files):
            shutil.rmtree(path)
        download(jdata)
        print(f'{folder}: fetched', file=sys.stderr)


def _data_root(args: argparse.Namespace) -> Path:
    # Test the argument before it becomes a Path: Path('') is Path('.'), so a
    # guard on the Path would pass and quietly use the working directory.
    given = args.data_root or os.environ.get('FWL_DATA')
    if not given:
        raise ResolutionError('no data root: pass --data-root or set FWL_DATA.')
    return Path(given)


def _cmd_key(args: argparse.Namespace) -> int:
    key = resolve_key()
    print(f'Cache key: {key}', file=sys.stderr)
    output = os.environ.get('GITHUB_OUTPUT')
    if output:
        with open(output, 'a', encoding='utf-8') as handle:
            handle.write(f'key={key}\n')
    else:
        print(f'key={key}')
    return 0


def _cmd_fetch(args: argparse.Namespace) -> int:
    fetch_all(_data_root(args))
    return 0


def _cmd_check(args: argparse.Namespace) -> int:
    incomplete = False
    for rel_dir, found, expected, state in check_restored(_data_root(args)):
        print(f'{rel_dir}: {found}/{expected} files {state}')
        if found != expected:
            incomplete = True

    if incomplete:
        print(
            'The data tree is missing files the registry pins or unpacked archive '
            'members, or holds a file with the wrong checksum, so it does not match '
            'its cache key. If it was restored on an exact-key hit, delete that '
            'cache entry: an exact hit is never re-saved, so the next run can then '
            'store a complete tree.',
            file=sys.stderr,
        )
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    """Run a subcommand and return its exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('key', help='print the cache key for the FWL data tree')
    for name, text in (
        ('fetch', 'fetch every dataset the key covers'),
        ('check', 'verify the data tree against the registry'),
    ):
        sub.add_parser(name, help=text).add_argument(
            '--data-root', default=None, help='defaults to FWL_DATA'
        )

    args = parser.parse_args(argv)
    handler = {'key': _cmd_key, 'fetch': _cmd_fetch, 'check': _cmd_check}[args.command]
    try:
        return handler(args)
    except ResolutionError as exc:
        print(f'error: {exc}', file=sys.stderr)
        return 1
    except Exception as exc:  # noqa: BLE001 -- reported, never swallowed
        # Anything fwl-io raises reaches here. Name it rather than let a
        # traceback stand in for the diagnostic this script promises.
        print(
            f'error: reading or fetching the data failed: {exc!r}. Check that Zenodo '
            'and OSF are reachable and that the installed fwl-mors, fwl-io and janus '
            'still expose the manifest, fetcher, check and downloaders this script reads.',
            file=sys.stderr,
        )
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
