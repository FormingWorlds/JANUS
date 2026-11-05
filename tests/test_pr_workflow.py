"""Tests for the data steps of tests.yaml and the time limit of the long-path tests.

Push and PR CI fetch Oak in their own step and run pytest offline. A cache prefix inside the
nightly restore-key ``fwl-data-`` would hand the nightly this small tree. See
docs/How-to/test.md.
"""

import re
from pathlib import Path

import pytest
import test_insert_stellar_long_path as long_path
import yaml
from fwl_io import OfflineDataError

from janus.utils import data as jdata

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

ROOT = '/home/runner/work/fwl_data'


def test_pr_ci_fills_the_oak_tree_before_an_offline_pytest():
    """Key, restore and fetch run in this order on one selection, one tree, before pytest."""
    path = Path(__file__).parents[1] / '.github' / 'workflows' / 'tests.yaml'
    workflow = yaml.safe_load(path.read_text(encoding='utf-8'))
    job = workflow['jobs']['test']
    steps = job['steps']
    sel = f'--key {jdata.spectral_file_key("Oak")}'
    i = next(i for i, s in enumerate(steps) if s.get('run', '').startswith('fwl-io cache key'))
    key, cache, fill = (
        {k: v.strip() if k == 'run' else v for k, v in s.items() if k != 'name'}
        for s in steps[i : i + 3]
    )
    assert key == {'id': 'oak-key', 'run': f'fwl-io cache key --prefix oak-data- {sel}'}
    assert re.fullmatch(r'actions/cache@(v[0-9]+|[0-9a-f]{40})', cache.pop('uses', ''))
    out = {k: f'${{{{ steps.oak-key.outputs.{k} }}}}' for k in ('key', 'restore-keys')}
    assert cache == {'with': {'path': ROOT, **out}}
    assert fill == {
        'timeout-minutes': 6,
        'run': f'fwl-io cache fetch {sel} --data-root {ROOT}',
    }

    tests = next(j for j, s in enumerate(steps) if 'pytest -m' in s.get('run', ''))
    assert i + 2 < tests
    assert not {'if', 'continue-on-error'} & steps[tests].keys()
    lines = [line.strip() for line in steps[tests]['run'].splitlines()]
    pytest_line = next(n for n, line in enumerate(lines) if line.startswith('pytest -m'))
    assert lines.index(f'export FWL_DATA="{ROOT}"') < pytest_line
    assert lines.index('export FWL_IO_OFFLINE=1') < pytest_line
    for scope in (workflow, job):
        assert not {'defaults', 'continue-on-error', 'if'} & scope.keys()
        assert not {'FWL_DATA', 'FWL_IO_OFFLINE'} & scope.get('env', {}).keys()


def test_offline_run_stops_with_a_named_error_when_oak_is_absent(monkeypatch, tmp_path):
    """Offline, a missing file raises OfflineDataError at once; no skip, no download."""
    monkeypatch.setenv('FWL_IO_OFFLINE', '1')
    monkeypatch.delenv('FWL_DATA_CACHE', raising=False)
    monkeypatch.setattr(jdata, 'FWL_DATA_DIR', tmp_path)
    with pytest.raises(OfflineDataError, match='offline mode is active'):
        jdata.DownloadSpectralFiles('Oak')
    assert not [p for p in tmp_path.rglob('Oak*') if p.is_file()]


def test_long_path_limit_holds_for_the_test_call_only():
    """A mirror read may take 60 s by itself, so a limit of 60 s on the fixture ends the fetch
    before it can reach the next mirror.
    """
    limit = next(m for m in long_path.pytestmark if m.name == 'timeout')
    assert limit.args == (60,)
    assert limit.kwargs == {'func_only': True}
