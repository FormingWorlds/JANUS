"""Tests for the data steps of ``.github/workflows/tests.yaml``.

Push and PR CI fetch the Oak spectral file in their own step and run pytest offline, so a
slow mirror cannot use up a test time limit, and a test that needs other data stops at once
with a named error. A cache prefix inside the nightly restore-key ``fwl-data-`` would hand
the nightly this small tree. See docs/How-to/test.md.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from janus.utils import data as jdata

yaml = pytest.importorskip('yaml')

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

ROOT = '/home/runner/work/fwl_data'


def test_pr_ci_fills_the_oak_tree_before_an_offline_pytest():
    """Key, restore and fetch run in this order on one selection, one tree, before pytest."""
    path = Path(__file__).parents[1] / '.github' / 'workflows' / 'tests.yaml'
    steps = yaml.safe_load(path.read_text(encoding='utf-8'))['jobs']['test']['steps']
    sel = f'--key {jdata.spectral_file_key("Oak")}'
    i = next(i for i, s in enumerate(steps) if s.get('run', '').startswith('fwl-io cache key'))
    key, cache, fill = steps[i : i + 3]
    assert key['run'].strip() == f'fwl-io cache key --prefix oak-data- {sel}'
    assert re.fullmatch(r'actions/cache@(v[0-9]+|[0-9a-f]{40})', cache['uses'])
    out = {k: f'${{{{ steps.{key["id"]}.outputs.{k} }}}}' for k in ('key', 'restore-keys')}
    assert cache['with'] == {'path': ROOT, **out}
    assert fill['run'].strip() == f'fwl-io cache fetch {sel} --data-root {ROOT}'
    assert not {'if', 'continue-on-error'} & fill.keys()

    tests = next(j for j, s in enumerate(steps) if 'pytest -m' in s.get('run', ''))
    assert i + 2 < tests
    assert not {'if', 'continue-on-error'} & steps[tests].keys()
    assert re.findall(r'export FWL_DATA="([^"]+)"', steps[tests]['run']) == [ROOT]
    assert re.findall(r'export FWL_IO_OFFLINE=(\S+)', steps[tests]['run']) == ['1']
