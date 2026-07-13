"""Tests for src/janus/set_socrates_env.py.

Exercises the import-time SOCRATES environment wiring: the exported
environment variables, the version string, and the hard failure when the
SOCRATES tree is absent. See docs/How-to/test.md.
"""

import importlib
import os
from pathlib import Path

import pytest

import janus.set_socrates_env as socenv
import janus.socrates

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def test_environment_variables_point_into_socrates_tree():
    """The module exports RAD_* variables rooted in an existing tree.

    Every exported path must sit inside SOCRATES_DIR and the version string
    read from the tree must be non-empty; an empty version would mean the
    wrong directory was accepted.
    """
    root = Path(os.environ['RAD_DIR'])
    assert root.exists()
    assert Path(os.environ['RAD_BIN']) == root / 'bin'
    assert Path(os.environ['RAD_DATA']) == root / 'data'
    assert Path(os.environ['RAD_SCRIPT']) == root / 'sbin'
    assert os.environ['LOCK_FILE'] == 'radiation_code.lock'
    assert len(socenv.SOCRATES_VERSION.strip()) > 0


def test_missing_socrates_tree_raises_on_import(monkeypatch, tmp_path):
    """A nonexistent SOCRATES location fails fast with RuntimeError.

    The guard runs at import time, so the module is re-executed with the
    location pointed at an empty directory entry; afterwards the module is
    restored so later tests keep a working environment.
    """
    bogus = tmp_path / 'no_socrates_here'
    monkeypatch.setattr(janus.socrates, 'SOCRATES_DIR', bogus, raising=True)
    try:
        with pytest.raises(RuntimeError, match='Cannot find SOCRATES'):
            importlib.reload(socenv)
    finally:
        monkeypatch.undo()
        importlib.reload(socenv)
    # Restored module still carries a valid version string.
    assert len(socenv.SOCRATES_VERSION.strip()) > 0
