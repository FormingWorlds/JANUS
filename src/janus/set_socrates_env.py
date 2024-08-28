"""Import this file to set environment variables for socrates.

Based on:
https://github.com/FormingWorlds/SOCRATES/blob/main/sbin/set_rad_env_tmp
"""

from __future__ import annotations

from .socrates import SOCRATES_DIR

import os
import sys
import zipfile
from pathlib import Path
import logging

import click
import platformdirs
import requests

log = logging.getLogger(__name__)

if not SOCRATES_DIR.exists():
    raise RuntimeError(f'Cannot find SOCRATES in this location: {SOCRATES_DIR}')

with open(SOCRATES_DIR / 'version') as f:
    SOCRATES_VERSION = f.readline()

log.debug('socrates location: %s' % SOCRATES_DIR)
log.debug('socrates version: %s' % SOCRATES_VERSION)

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
