"""Shared pytest hooks for the JANUS test suite."""


def pytest_sessionstart(session):
    """Build the matplotlib font cache before the first plotting test and its time limit."""
    from matplotlib import font_manager  # noqa: F401
