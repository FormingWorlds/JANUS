import os
from pathlib import Path
import sys

FWL_DATA = Path(os.environ.get('FWL_DATA'))
SOCRATES_DIR = FWL_DATA / 'socrates'

def set_rad_env():
    """Set environment variables for socrates.

    Based on:
    https://github.com/FormingWorlds/SOCRATES/blob/main/sbin/set_rad_env_tmp
    """
    assert SOCRATES_DIR.exists()
    
    with open(SOCRATES_DIR / 'version') as f:
        version = f.readline()

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