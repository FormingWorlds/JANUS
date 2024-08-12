import os
from pathlib import Path
import logging

import platformdirs
from osfclient.api import OSF

logger = logging.getLogger(__name__)

FWL_DATA_DIR = Path(os.environ.get('FWL_DATA', platformdirs.user_data_dir('fwl_data')))

logger.info(f'FWL data location: {FWL_DATA_DIR}')

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
    for file in storage.files:
        for folder in folders:
            if not file.path[1:].startswith(folder):
                continue
            parts = file.path.split('/')[1:]
            target = Path(data_dir, *parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            print(f'Downloading {file.path}...')
            with open(target, 'wb') as f:
                file.write_to(f)
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
    project = osf.project(project_id)
    storage = project.storage('osfstorage')

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
    project = osf.project(project_id)
    storage = project.storage('osfstorage')

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
