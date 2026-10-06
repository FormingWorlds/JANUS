"""Reference data for JANUS: SOCRATES spectral files and stellar spectra.

Both come from the shared manifest that fwl-io ships. fwl-io fetches each dataset from its
Zenodo record, falls back to the DataverseNL mirror the manifest pins, checks every file
against the committed registry, and places the files in a version directory
``<FWL_DATA>/<key-as-path>/r<record-id>``. Callers resolve that directory through
:func:`spectral_file_dir` and :func:`stellar_spectra_dir` rather than joining it by hand.
"""

import logging
import os
from pathlib import Path

import platformdirs

log = logging.getLogger('fwl.' + __name__)

FWL_DATA_DIR = Path(os.environ.get('FWL_DATA', platformdirs.user_data_dir('fwl_data')))

STELLAR_SPECTRA_NAMED = 'star.spectra.named'

basic_list = (
    'Dayspring/256',
    'Frostflow/256',
    'Oak/318',
)


def GetFWLData() -> Path:
    """
    Get path to FWL data directory on the disk
    """
    return FWL_DATA_DIR.absolute()


def _shared_datasets() -> dict:
    """Return the datasets of the fwl-io shared manifest, keyed by manifest key."""
    from fwl_io import load_manifest
    from fwl_io.manifest import shared_manifest_path

    return {ds.key: ds for ds in load_manifest(shared_manifest_path())}


def _fetcher(key: str):
    """Build the fwl-io fetcher of one shared dataset below the FWL data directory."""
    from fwl_io import create_fetcher

    ds = _shared_datasets()[key]
    return create_fetcher(
        subdir=ds.subdir,
        zenodo=ds.zenodo,
        dataverse=ds.dataverse,
        registry=ds.registry(),
        data_root=GetFWLData(),
        extract=ds.extract,
    )


def spectral_file_key(group: str, bands: int | str | None = None, default=None) -> str:
    """Return the manifest key of a spectral file.

    Parameters
    ----------
    group : str
        Spectral file group, e.g. ``"Dayspring"`` or ``"Oak"``.
    bands : int or str, optional
        Number of bands. A group the manifest declares with one band count
        resolves to that count; a different count given here is ignored with a
        warning.
    default : int, optional
        Band count used for a group with several counts when ``bands`` is None.

    Returns
    -------
    str
        Key such as ``atmos_clim.spectral_files.dayspring.256``.

    Raises
    ------
    ValueError
        The manifest declares no such group, or a group with several band counts
        gets no count or a count it does not declare.
    """
    prefix = f'atmos_clim.spectral_files.{group.lower()}.'
    counts = sorted(
        (k.removeprefix(prefix) for k in _shared_datasets() if k.startswith(prefix)), key=int
    )
    if not counts:
        raise ValueError(f'No spectral file group {group!r} in the installed fwl-io manifest')
    if len(counts) == 1:
        if bands is not None and str(bands) != counts[0]:
            log.warning(f'{group} has only {counts[0]} bands; ignoring the requested {bands}')
        return f'{prefix}{counts[0]}'
    bands = default if bands is None else bands
    if bands is None or str(bands) not in counts:
        raise ValueError(f'{group} needs a band count, one of {counts}; got {bands}')
    return f'{prefix}{bands}'


def spectral_file_dir(group: str, bands: int | str | None = None) -> Path:
    """Return the directory that holds a spectral file, e.g. ``.../Oak.sf``.

    Resolving the path downloads nothing (call :func:`DownloadSpectralFiles`), but it
    creates the FWL data directory if it is absent.
    """
    return _fetcher(spectral_file_key(group, bands)).target_dir


def stellar_spectra_dir() -> Path:
    """Return the directory that holds the named stellar spectra, e.g. ``.../sun.txt``.

    Like :func:`spectral_file_dir`, it creates the FWL data directory if it is absent.
    """
    return _fetcher(STELLAR_SPECTRA_NAMED).target_dir


def DownloadStellarSpectra():
    """
    Download the named stellar spectra through fwl-io.
    """
    _fetcher(STELLAR_SPECTRA_NAMED).fetch_all()


def DownloadSpectralFiles(fname: str = '', nband: int | None = None):
    """
    Download spectral files through fwl-io.

    Inputs :
        - fname (optional) :    group name, e.g. "Dayspring" or "Oak"
                                if not provided download the basic list
        - nband (optional) :    number of bands = 16, 48, 256 (default), 4096
                                for a group with several band counts; a single-band
                                group such as Oak ignores it
    """
    pairs = [folder.split('/') for folder in basic_list] if not fname else [(fname, nband)]
    for key in [spectral_file_key(group, bands, default=256) for group, bands in pairs]:
        _fetcher(key).fetch_all()
