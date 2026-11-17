"""Tests for .github/scripts/docs_consistency/analyze.py.

Covers the readability step that rewords suggested fixes after the analysis: which
rewrites are used, and that a failed readability call keeps the original fixes. The
CLI call itself is mocked. See docs/How-to/test.md.
"""

import importlib
import json
import sys
from pathlib import Path

import pytest

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

SCRIPT_DIR = Path(__file__).resolve().parents[1] / '.github' / 'scripts' / 'docs_consistency'
sys.path.insert(0, str(SCRIPT_DIR))
analyze = importlib.import_module('analyze')


def fixed(fid, old='old words', new='new words'):
    """A finding with a suggested fix."""
    return {
        'id': fid,
        'type': 'inconsistency',
        'severity': 'minor',
        'doc_excerpt': f'The {old} here.',
        'code_excerpt': '    return y',
        'description': 'The doc and the code disagree.',
        'suggested_fix': {'old_text': old, 'new_text': new},
    }


def test_rewrites_replace_only_matching_unambiguous_fixes():
    """A rewrite lands on the one fix with its id; unknown, duplicated, empty and
    no-op rewrites are ignored, and findings without a fix are never touched."""
    findings = [fixed('a'), fixed('b'), fixed('dup'), fixed('dup'), {'id': 'gap'}]
    rewrites = [
        {'id': 'a', 'new_text': 'clearer words'},
        {'id': 'b', 'new_text': 'old words'},  # would undo the fix
        {'id': 'dup', 'new_text': 'which one?'},
        {'id': 'gap', 'new_text': 'no fix to reword'},
        {'id': 'missing', 'new_text': 'unknown id'},
        {'id': 'a', 'new_text': '   '},
    ]
    assert analyze.apply_rewrites(findings, rewrites) == 1
    assert findings[0]['suggested_fix']['new_text'] == 'clearer words'
    assert [f['suggested_fix']['new_text'] for f in findings[1:4]] == ['new words'] * 3
    assert 'suggested_fix' not in findings[4]
    # Edge case: nothing to reword.
    assert analyze.apply_rewrites([], rewrites) == 0


def test_readability_prompt_carries_each_fix_with_its_context():
    """The prompt lists every fix with its quotes, so the rewording can stay true."""
    prompt = analyze.build_readability_prompt([fixed('a', old='alpha=1', new='alpha=5')])
    items = json.loads(prompt[prompt.index('[') :])
    assert items == [
        {
            'id': 'a',
            'doc_excerpt': 'The alpha=1 here.',
            'code_excerpt': '    return y',
            'description': 'The doc and the code disagree.',
            'old_text': 'alpha=1',
            'new_text': 'alpha=5',
        }
    ]
    assert '{{FIXES}}' not in prompt


@pytest.mark.parametrize('readability_fails', [False, True], ids=['reworded', 'call-fails'])
def test_main_writes_findings_whether_or_not_the_readability_call_works(
    tmp_path, monkeypatch, readability_fails
):
    """The readability call only rewords: when it fails, the run still writes the
    findings with their original fixes instead of stopping."""
    calls = []

    def fake_run(prompt, schema, timeout_s):
        calls.append(schema)
        if schema is analyze.FINDINGS_SCHEMA:
            return {'findings': [fixed('a'), {'id': 'gap', 'type': 'gap'}]}
        if readability_fails:
            raise RuntimeError('claude CLI did not finish within 300 s')
        return {'fixes': [{'id': 'a', 'new_text': 'clearer words'}]}

    monkeypatch.setenv('CLAUDE_CODE_OAUTH_TOKEN', 'test')
    monkeypatch.setattr(analyze, 'REPO_ROOT', tmp_path)
    monkeypatch.setattr(analyze, 'build_prompt', lambda: 'prompt')
    monkeypatch.setattr(analyze, 'run_claude', fake_run)
    analyze.main()

    written = json.loads((tmp_path / 'analysis.json').read_text())
    expected = 'new words' if readability_fails else 'clearer words'
    assert written[0]['suggested_fix']['new_text'] == expected
    assert len(written) == 2 and calls == [analyze.FINDINGS_SCHEMA, analyze.READABILITY_SCHEMA]


def test_main_stops_when_the_analysis_call_fails(tmp_path, monkeypatch):
    """Without findings there is nothing to publish, so a failed analysis call exits
    non-zero and writes no analysis.json."""

    def fake_run(prompt, schema, timeout_s):
        raise RuntimeError('claude CLI exited 1')

    monkeypatch.setenv('CLAUDE_CODE_OAUTH_TOKEN', 'test')
    monkeypatch.setattr(analyze, 'REPO_ROOT', tmp_path)
    monkeypatch.setattr(analyze, 'build_prompt', lambda: 'prompt')
    monkeypatch.setattr(analyze, 'run_claude', fake_run)
    with pytest.raises(SystemExit) as exc:
        analyze.main()
    assert exc.value.code == 1
    assert not (tmp_path / 'analysis.json').exists()
