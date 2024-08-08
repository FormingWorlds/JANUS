from __future__ import annotations

import os
import subprocess as sp
import zipfile
from pathlib import Path

import click
import platformdirs
import requests

SOCRATES_DATA_DIR = Path(platformdirs.user_data_dir('socrates'))
SOCRATES_DIR = Path(os.environ.get('SOCRATES', SOCRATES_DATA_DIR / 'SOCRATES'))


def _set_permissions(drc: Path):
    """Set executable flag for scripts."""
    # Set executable permissions for make script
    (drc / 'make' / 'mkdep').chmod(mode=33261)


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


def download_socrates(ref: str = 'main'):
    """Version can be 'main' or the git commit hash."""
    filename = 'socrates.zip'
    path = 'refs/heads/main' if ref == 'main' else ref

    _download(
        url=f'https://github.com/nichollsh/SOCRATES/archive/{path}.zip',
        filename=filename,
    )

    SOCRATES_DATA_DIR.mkdir(exist_ok=True, parents=True)

    subdir = f'SOCRATES-{ref}'
    target_dir = SOCRATES_DATA_DIR / subdir

    click.echo(f'Extracting to {target_dir}')

    with zipfile.ZipFile(filename, 'r') as zip_ref:
        zip_ref.extractall(SOCRATES_DATA_DIR)

    _set_permissions(drc=target_dir)

    sp.run(['bash', 'configure'], cwd=target_dir)
    sp.run(['bash', 'build_code'], cwd=target_dir)

    symlink = SOCRATES_DATA_DIR / 'SOCRATES'
    symlink.unlink(missing_ok=True)
    symlink.symlink_to(target_dir)

    print(f'SOCRATES downloaded to {target_dir}')
