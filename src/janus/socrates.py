from __future__ import annotations

import os
import sys
import zipfile
from pathlib import Path

import click
import platformdirs
import requests

SOCRATES_DIR = Path(os.environ.get('SOCRATES', platformdirs.user_data_dir('socrates')))


def _set_permissions(drc: Path):
    """Set executable flag for scripts."""
    for p in drc.iterdir():
        if p.suffix:
            continue
        if p.stem.startswith('Mk_') or p.stem.startswith('Make'):
            continue
        # set executable flag
        p.chmod(mode=33261)


def _download(*, url: str, filename: str):
    """Download file from url."""
    chunk_size = 1024 * 8

    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        total_size = int(r.headers.get('Content-Length', 0))
        with (
            open(filename, 'wb') as f,
            click.progressbar(label='Downloading', length=total_size) as pbar,
        ):
            for chunk in r.iter_content(chunk_size=chunk_size):
                f.write(chunk)
                pbar.update(chunk_size)


def download_socrates():
    version = 'x.y.z'
    filename = 'socrates.zip'

    _download(
        url=f'https://github.com/FormingWorlds/SOCRATES/releases/download/{version}/socrates-Linux.zip',
        filename=filename,
    )

    SOCRATES_DIR.mkdir(exist_ok=True, parents=True)

    is_empty = not any(SOCRATES_DIR).iterdir()
    if not is_empty:
        raise RuntimeError(f'SOCRATES directory is not empty: {SOCRATES_DIR}')

    with zipfile.ZipFile(filename, 'r') as zip_ref:
        zip_ref.extractall(SOCRATES_DIR)

    _set_permissions(SOCRATES_DIR / 'bin')
    _set_permissions(SOCRATES_DIR / 'sbin')

    filename.unlink()

    print(f'SOCRATES downloaded to {SOCRATES_DIR}')
