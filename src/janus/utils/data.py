import os
import re
import time
from pathlib import Path
import logging

import platformdirs
import requests
from osfclient.api import OSF

log = logging.getLogger("fwl."+__name__)

OSF_RETRY_ATTEMPTS = 3
OSF_RETRY_DELAYS = (15, 45)

# Most osfclient failures surface as a plain RuntimeError with a
# "...status code {N}..." message, transient or not; a 401 instead raises
# osfclient's own UnauthorizedException, which this module does not import
# or match, so it is never treated as retryable. Only retry the status
# codes that are actually worth retrying.
_OSF_RETRYABLE_STATUS = {429, 500, 502, 503, 504}
_OSF_STATUS_CODE_RE = re.compile(r'status code (\d+)')

def _is_transient_osf_error(exc: RuntimeError) -> bool:
    match = _OSF_STATUS_CODE_RE.search(str(exc))
    return bool(match) and int(match.group(1)) in _OSF_RETRYABLE_STATUS

def _osf_retry(func):
    """
    Call `func`, retrying on a transient OSF failure (a 429/5xx RuntimeError
    from osfclient, or a network-level error while streaming a response)
    with a bounded backoff, so a flaky OSF response does not fail the
    download outright. A non-transient RuntimeError (e.g. a 404, or a
    RuntimeError whose message carries no status code) is re-raised
    immediately rather than retried.
    """
    for attempt in range(OSF_RETRY_ATTEMPTS):
        try:
            return func()
        except RuntimeError as exc:
            if not _is_transient_osf_error(exc) or attempt == OSF_RETRY_ATTEMPTS - 1:
                raise
            last_exc = exc
        except requests.exceptions.RequestException as exc:
            if attempt == OSF_RETRY_ATTEMPTS - 1:
                raise
            last_exc = exc
        log.warning(f'OSF request failed (attempt {attempt + 1}/'
                    f'{OSF_RETRY_ATTEMPTS}): {last_exc!r}, retrying...')
        time.sleep(OSF_RETRY_DELAYS[min(attempt, len(OSF_RETRY_DELAYS) - 1)])

FWL_DATA_DIR = Path(os.environ.get('FWL_DATA', platformdirs.user_data_dir('fwl_data')))

log.debug(f'FWL data location: {FWL_DATA_DIR}')

basic_list = (
        "Dayspring/256",
        "Frostflow/256",
        "Legacy",
        "Mallard",
        "Oak",
        "Reach",
        )

def download_folder(*, storage, folders: list[str], data_dir: Path):
    """
    Download a specific folder in the OSF repository

    Inputs :
        - storage : OSF storage name
        - folders : folder names to download
        - data_dir : local repository where data are saved
    """
    files = _osf_retry(lambda: list(storage.files))
    for file in files:
        for folder in folders:
            if not file.path[1:].startswith(folder):
                continue
            parts = file.path.split('/')[1:]
            target = Path(data_dir, *parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            log.info(f'Downloading {file.path}...')

            def _write(file=file, target=target):
                with open(target, 'wb') as f:
                    file.write_to(f)

            _osf_retry(_write)
            break


def GetFWLData() -> Path:
    """
    Get path to FWL data directory on the disk
    """
    return FWL_DATA_DIR.absolute()

def DownloadStellarSpectra():
    """
    Download stellar spectra
    """
    #project ID of the stellar spectra on OSF
    project_id = '8r2sw'
    folder_name = 'Named'

    osf = OSF()
    project = _osf_retry(lambda: osf.project(project_id))
    storage = _osf_retry(lambda: project.storage('osfstorage'))

    data_dir = GetFWLData() / "stellar_spectra"
    data_dir.mkdir(parents=True, exist_ok=True)

    if not (data_dir / folder_name).exists():
        print(f"Downloading stellar spectra to {data_dir}")
        download_folder(storage=storage, folders=[folder_name], data_dir=data_dir)


def DownloadSpectralFiles(fname: str="",nband: int=256):
    """
    Download spectral files data

    Inputs :
        - fname (optional) :    folder name, i.e. "/Dayspring"
                                if not provided download all the basic list
        - nband (optional) :    number of band = 16, 48, 256, 4096
                                (only relevant for Dayspring, Frostflow and Honeyside)
    """
    #project ID of the spectral files on OSF
    project_id = 'vehxg'

    #Create spectral file data repository if not existing
    data_dir = GetFWLData() / "spectral_files"
    data_dir.mkdir(parents=True, exist_ok=True)

    #Link with OSF project repository
    osf = OSF()
    project = _osf_retry(lambda: osf.project(project_id))
    storage = _osf_retry(lambda: project.storage('osfstorage'))

    #If no folder specified download all basic list
    if not fname:
        folder_list = basic_list
    elif fname in ("Dayspring", "Frostflow", "Honeyside"):
        folder_list = [fname + "/" + str(nband)]
    elif fname in ("Kynesgrove","Legacy","Mallard","Oak","Reach","stellar_spectra"):
        folder_list = [fname]
    else:
        raise ValueError(f"Unrecognised folder name: {fname}")

    folders = [folder for folder in folder_list if not (data_dir / folder).exists()]

    if folders:
        print(f"Downloading SOCRATES spectral files to {data_dir}")
        download_folder(storage=storage, folders=folders, data_dir=data_dir)
