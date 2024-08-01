import os
import sys
from pathlib import Path

import click
import platformdirs
import requests

import zipfile


SOCRATES_DIR = Path(os.environ.get('SOCRATES', platformdirs.user_data_dir('socrates')))


def set_rad_env():
    """Set environment variables for socrates.

    Based on:
    https://github.com/FormingWorlds/SOCRATES/blob/main/sbin/set_rad_env_tmp
    """
    if not SOCRATES_DIR.exists():
        raise RuntimeError(f'Cannot find SOCRATES in this location: {SOCRATES_DIR}')

    with open(SOCRATES_DIR / 'version') as f:
        version = f.readline()

    print(f'socrates location: {SOCRATES_DIR}')
    print(f'socrates version: {version}')

    sep = os.pathsep

    os.environ['RAD_DIR'] = str(SOCRATES_DIR)
    os.environ['RAD_BIN'] = str(SOCRATES_DIR / 'bin')
    os.environ['RAD_DATA'] = str(SOCRATES_DIR / 'data')
    os.environ['RAD_SCRIPT'] = str(SOCRATES_DIR / 'sbin')
    os.environ['LOCK_FILE'] = 'radiation_code.lock'
    os.environ['PATH'] = str(SOCRATES_DIR / 'bin') + sep + os.environ['PATH']
    os.environ['PATH'] = str(SOCRATES_DIR / 'sbin') + sep + os.environ['PATH']
    os.environ['MANPATH'] = str(SOCRATES_DIR / 'man') + sep + os.environ.get('MANPATH', '')
    sys.path.append(str(SOCRATES_DIR / 'python'))

    os.environ['LD_LIBRARY_PATH'] = 'netcdfff' + sep + os.environ.get('LD_LIBRARY_PATH', '')


def download_socrates():
    filename = 'socrates.zip'
    version = 'x.y.z'
    url = f'https://github.com/FormingWorlds/SOCRATES/releases/download/{version}/socrates-Linux.zip'

    with requests.get(url, stream=True) as r:
        total_size = int(r.headers.get('Content-Length', 0))
        chunk_size = 1024 * 8

        r.raise_for_status()

        with (
            open(filename, 'wb') as f,
            click.progressbar(label='Downloading', length=total_size) as pbar,
        ):
            for chunk in r.iter_content(chunk_size=chunk_size):
                f.write(chunk)
                pbar.update(chunk_size)

    SOCRATES_DIR.mkdir(exist_ok=True, parents=True)

    target_dir = SOCRATES_DIR / 'socrates'

    if target_dir.exists():
        raise RuntimeError(f'SOCRATES directory already exists: {SOCRATES_DIR}')

    with zipfile.ZipFile(filename, 'r') as zip_ref:
        zip_ref.extractall(SOCRATES_DIR)

    print(f'SOCRATES downloaded to {SOCRATES_DIR}')
    filename.unlink()