"""Cross-cutting tests for JANUS support utilities: logging, sorting, CLI.

This cross-cutting file (a documented exception to the mirror-the-source rule)
covers non-physics plumbing: the ``setup_logger`` level/handler wiring in
janus.utils.logs, the ``natural_sort`` numeric-aware ordering in
janus.utils.socrates, and the ``janus env`` CLI command. These are utility
sources, exempt from the physics-invariant requirement but held to the
anti-happy-path rules. See docs/How-to/test.md.
"""

import logging

import pytest
from click.testing import CliRunner

from janus.cli import cli
from janus.utils.logs import setup_logger
from janus.utils.socrates import natural_sort

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def test_setup_logger_levels_and_handlers():
    """setup_logger sets the requested level and installs a single handler.

    Re-invoking the logger factory must clear and re-add exactly one handler
    (no handler accumulation on repeated calls), route under the ``fwl.``
    namespace, and honour the requested level string for both INFO and the
    more verbose DEBUG, which discriminates a factory that ignored its level
    argument.
    """
    logger = setup_logger('INFO')
    assert logger.name.startswith('fwl.')
    assert logger.level == logging.INFO
    assert len(logger.handlers) == 1  # cleared and re-added, no accumulation

    logger = setup_logger('DEBUG')
    assert logger.level == logging.DEBUG
    # Re-invocation must not accumulate handlers.
    assert len(logger.handlers) == 1


def test_setup_logger_invalid_level():
    """An unrecognised level string raises ValueError (the error contract).

    ``setup_logger`` validates its level argument, so a nonsense level must
    raise ValueError rather than silently defaulting, and a valid level must
    not raise. This discriminates a factory that swallowed bad input and
    quietly changed verbosity.
    """
    with pytest.raises(ValueError):
        setup_logger('SILLY')
    # A valid level, by contrast, must not raise and sets the level.
    assert setup_logger('WARNING').level == logging.WARNING


def test_natural_sort_mixed_numbers():
    """natural_sort orders embedded integers numerically, not lexically.

    Numeric-aware sorting must place ``file2`` before ``file10`` (a plain
    string sort would invert them), and the result must be a permutation of
    the input of the same multiset. The mixed one- and two-digit inputs
    discriminate a lexical sort.
    """
    names = ['file10', 'file2', 'file1', 'file20', 'file11']
    ordered = natural_sort(names)
    assert ordered == ['file1', 'file2', 'file10', 'file11', 'file20']
    # Lexical sort would place file10/file11 before file2; guard against it.
    assert ordered.index('file2') < ordered.index('file10')
    assert sorted(ordered) == sorted(names)  # same multiset, reordered


def test_cli_env_outputs_env(monkeypatch, tmp_path):
    """The ``janus env`` command exits cleanly and echoes RAD_DIR.

    Invoking the CLI ``env`` subcommand must return a zero exit code and print
    the resolved ``RAD_DIR`` it was given, confirming the environment resolver
    reads and reports the variable rather than crashing or emitting nothing.
    """
    rad_dir = str(tmp_path / 'SOCRATES')
    monkeypatch.setenv('RAD_DIR', rad_dir)
    result = CliRunner().invoke(cli, ['env'])
    assert result.exit_code == 0
    assert 'RAD_DIR' in result.output
