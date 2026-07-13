"""Tests for src/janus/utils/logs.py.

Exercises the terminal logger factory: level validation, idempotent handler
setup, and the installed excepthook that routes uncaught exceptions and
keyboard interrupts through the logger. See docs/How-to/test.md.
"""

import logging
import sys

import pytest

from janus.utils.logs import setup_logger

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.fixture(autouse=True)
def _restore_excepthook():
    """setup_logger installs a global excepthook; restore it afterwards."""
    original = sys.excepthook
    yield
    sys.excepthook = original


def test_setup_logger_levels_and_idempotency():
    """The logger honours the requested level and re-setup does not stack
    handlers.

    Calling the factory twice must leave exactly one stream handler; an
    invalid level string is the documented error contract and must raise
    without attaching a new handler (the factory clears the old handlers
    first, so a raise leaves the logger handler-free until reconfigured).
    """
    lg = setup_logger('DEBUG')
    assert lg.level == logging.DEBUG
    assert len(lg.handlers) == 1

    lg2 = setup_logger('warning')  # case-insensitive
    assert lg2 is lg
    assert lg2.level == logging.WARNING
    assert len(lg2.handlers) == 1

    with pytest.raises(ValueError, match='Invalid log level'):
        setup_logger('LOUD')
    # The factory clears handlers before validating, so the failed call
    # leaves the logger handler-free until reconfigured.
    assert len(lg2.handlers) == 0


def test_excepthook_routes_exceptions_through_logger(monkeypatch, caplog):
    """Uncaught exceptions log as critical; KeyboardInterrupt logs as error.

    The interrupt branch must also delegate to the original system hook so
    the process still terminates the standard way; a stub records that
    delegation.
    """
    lg = setup_logger('DEBUG')
    delegated = []
    monkeypatch.setattr(
        sys, '__excepthook__', lambda *args: delegated.append(args), raising=True
    )

    with caplog.at_level(logging.DEBUG, logger=lg.name):
        try:
            raise ValueError('boom')
        except ValueError:
            sys.excepthook(*sys.exc_info())

        try:
            raise KeyboardInterrupt()
        except KeyboardInterrupt:
            sys.excepthook(*sys.exc_info())

    messages = [r.getMessage() for r in caplog.records]
    levels = [r.levelno for r in caplog.records]
    assert 'Uncaught exception' in messages
    assert logging.CRITICAL in levels
    assert 'KeyboardInterrupt' in messages
    # The interrupt branch delegated exactly once to the system hook.
    assert len(delegated) == 1
    assert delegated[0][0] is KeyboardInterrupt
