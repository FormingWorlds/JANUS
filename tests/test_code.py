import logging
import pytest
from janus.utils.logs import setup_logger
from pathlib import Path
from janus.utils.socrates import natural_sort, CleanOutputDir
import numpy as np
import netCDF4 as net
from janus.utils import nctools
from click.testing import CliRunner
from janus.cli import cli

def test_setup_logger_levels_and_handlers():
    logger = setup_logger("INFO")
    assert logger.name.startswith("fwl.")
    assert logger.level == logging.INFO
    assert len(logger.handlers) == 1  # cleared and re-added

    logger = setup_logger("DEBUG")
    assert logger.level == logging.DEBUG

def test_setup_logger_invalid_level():
    with pytest.raises(ValueError):
        setup_logger("SILLY")

def test_natural_sort_mixed_numbers():
    names = ["file10", "file2", "file1", "file20", "file11"]
    assert natural_sort(names) == ["file1", "file2", "file10", "file11", "file20"]

def test_cli_env_outputs_env(monkeypatch, tmp_path):
    monkeypatch.setenv('RAD_DIR', str(tmp_path / 'SOCRATES'))
    result = CliRunner().invoke(cli, ['env'])
    assert result.exit_code == 0
    assert 'RAD_DIR' in result.output

