"""Tests for src/janus/cli.py.

Exercises the click command tree: option parsing and forwarding for the
three download subcommands, the bare group help path, the env report, and
the unknown-option error contract. See docs/How-to/test.md.
"""

from unittest.mock import patch

import pytest
from click.testing import CliRunner

from janus.cli import cli

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def test_download_subcommands_forward_options():
    """Each download subcommand forwards its parsed options to the worker.

    The spectral command must pass the name and the non-default band count,
    stellar takes no options, and socrates must forward a custom git ref;
    asserting the exact kwargs discriminates a swapped or dropped option.
    """
    runner = CliRunner()

    with patch('janus.utils.data.DownloadSpectralFiles') as mock_spec:
        result = runner.invoke(cli, ['download', 'spectral', '-n', 'Oak', '-b', '4096'])
    assert result.exit_code == 0
    mock_spec.assert_called_once_with(fname='Oak', nband=4096)

    with patch('janus.utils.data.DownloadStellarSpectra') as mock_stel:
        result = runner.invoke(cli, ['download', 'stellar'])
    assert result.exit_code == 0
    mock_stel.assert_called_once_with()

    with patch('janus.socrates.download_socrates') as mock_soc:
        result = runner.invoke(cli, ['download', 'socrates', '-r', 'abc123'])
    assert result.exit_code == 0
    mock_soc.assert_called_once_with(ref='abc123')


def test_spectral_band_default_and_group_help():
    """The band count defaults to 256 and the bare group prints help.

    The default matters because it selects which resolution folder is
    fetched; the bare `download` invocation exercises the group body and
    must list its subcommands rather than fail.
    """
    runner = CliRunner()

    with patch('janus.utils.data.DownloadSpectralFiles') as mock_spec:
        result = runner.invoke(cli, ['download', 'spectral'])
    assert result.exit_code == 0
    mock_spec.assert_called_once_with(fname=None, nband=256)

    # A bare group prints its usage and exits non-zero; the exact code
    # differs across click releases, so only the sign is pinned.
    result = runner.invoke(cli, ['download'])
    assert result.exit_code != 0
    for sub in ('spectral', 'stellar', 'socrates'):
        assert sub in result.output


def test_env_reports_locations_and_bad_option_fails():
    """env prints both data locations; an unknown option is rejected.

    The env command must surface both RAD_DIR and FWL_DATA so users can
    diagnose a broken environment; the unknown-option path is the CLI's
    error contract and must exit non-zero without invoking any download.
    """
    runner = CliRunner()

    from janus.socrates import SOCRATES_DIR
    from janus.utils.data import FWL_DATA_DIR

    result = runner.invoke(cli, ['env'])
    assert result.exit_code == 0
    # Pin the resolved paths themselves, not the label prose, so a reworded
    # echo string cannot hide a wrong location.
    assert str(SOCRATES_DIR) in result.output
    assert str(FWL_DATA_DIR) in result.output

    with patch('janus.utils.data.DownloadSpectralFiles') as mock_spec:
        result = runner.invoke(cli, ['download', 'spectral', '--bogus'])
    assert result.exit_code != 0
    mock_spec.assert_not_called()
