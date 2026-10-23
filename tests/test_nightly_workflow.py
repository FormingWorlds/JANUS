"""Tests for the data cache steps of ``.github/workflows/nightly.yml``.

The nightly keys, restores, checks, fetches and checks the data the tests read with
``fwl-io cache``: every fwl-mors dataset plus the Oak spectral file and the named
stellar spectra. A key over other data than the fetch fills, a fetch into another tree
than the tests read, or a check that is skipped or allowed to fail, fails nothing: the
nightly stays green and serves a stale tree. The tests pin the cache steps to their
text, up to the actions/cache version and the outer whitespace of a run, and the
selection to the datasets the test helpers read.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from janus.utils import data as jdata

yaml = pytest.importorskip('yaml')

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

REPO = Path(__file__).parents[1]

NIGHTLY_SEL = '--provider mors --key atmos_clim.spectral_files.oak.318 --key star.spectra.named'
NIGHTLY_ROOT = '/home/runner/work/fwl_data'


def test_nightly_keys_checks_fetches_and_checks_one_selection_in_one_tree():
    """The nightly keys and restores the data, checks an exact hit, then fills and checks it."""
    workflow = yaml.safe_load(
        (REPO / '.github' / 'workflows' / 'nightly.yml').read_text(encoding='utf-8')
    )
    jobs = list(workflow['jobs'].values())
    steps = [step for job in jobs for step in job['steps']]
    sel, root = NIGHTLY_SEL, NIGHTLY_ROOT
    # Fixed step bodies: a shell, an extra command or flag, or another condition would let
    # an unchecked tree be saved, or a damaged exact hit pass.
    i = next(i for i, s in enumerate(steps) if s.get('run', '').startswith('fwl-io cache key'))
    key, cache, hit, fill = (
        {k: v.strip() if k == 'run' else v for k, v in s.items() if k != 'name'}
        for s in steps[i : i + 4]
    )
    assert key == {'id': 'cachekey', 'run': f'fwl-io cache key --prefix fwl-data- {sel}'}
    out = {k: f'${{{{ steps.cachekey.outputs.{k} }}}}' for k in ('key', 'restore-keys')}
    # A major-version tag or a full commit SHA of actions/cache, never a branch such as main.
    assert re.fullmatch(r'actions/cache@(v[0-9]+|[0-9a-f]{40})', cache.pop('uses', ''))
    assert cache == {'id': 'cache-fwl-data', 'with': {'path': root, **out}}
    assert hit == {
        'if': "steps.cache-fwl-data.outputs.cache-hit == 'true'",
        'run': f'fwl-io cache check {sel} --data-root {root}',
    }
    assert fill == {
        'run': '\n'.join(
            f'fwl-io cache {c} {sel} --data-root {root}' for c in ('fetch', 'check')
        )
    }
    # The tests read the tree the steps above filled, after them, on every run.
    tests = next(j for j, s in enumerate(steps) if 'export FWL_DATA=' in s.get('run', ''))
    assert i + 3 < tests
    assert not {'if', 'continue-on-error'} & steps[tests].keys()
    assert re.findall(r'export FWL_DATA="([^"]+)"', steps[tests]['run']) == [root]
    for scope in (workflow, *jobs):
        assert not {'defaults', 'continue-on-error', 'if'} & scope.keys()
        assert 'FWL_DATA' not in scope.get('env', {})


def test_nightly_selection_is_mors_and_the_datasets_the_test_helpers_read():
    """The nightly caches every fwl-mors dataset plus the Oak file and the named spectra."""
    selection = re.findall(r'--(provider|key) (\S+)', NIGHTLY_SEL)
    assert [v for kind, v in selection if kind == 'provider'] == ['mors']
    assert {v for kind, v in selection if kind == 'key'} == {
        jdata.spectral_file_key('Oak'),
        jdata.STELLAR_SPECTRA_NAMED,
    }
