"""Tests for .github/scripts/docs_consistency/analyze.py.

Covers the readability step that rewords suggested fixes after the analysis: which
rewrites are used, that a rewording must keep the numbers, math, code spans, footnote
references and line breaks of the fix, and that a failed or malformed readability
reply keeps the original fixes. Also covers the error paths of the CLI call, with
subprocess.run patched. See docs/How-to/test.md.
"""

import importlib
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

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


ORIGINAL = 'T_r is 290 K, see $T_r$ [^cite-x] and `minT`.'
# Text wrapped by hand inside a `!!! note`, indented by four spaces.
NOTE = 'T_r is 290 K,\n    see $T_r$ [^cite-x] and `minT`.'


@pytest.mark.parametrize(
    'original, rewrite',
    [
        (ORIGINAL, 'T_r is set to the default.'),
        (ORIGINAL, 'T_r is 300 K, see $T_r$ [^cite-x] and `minT`.'),
        (ORIGINAL, 'T_r is 290 K, see T_r [^cite-x] and `minT`.'),
        (ORIGINAL, 'T_r is 290 K, see $T_r$ [^cite-x] and minT.'),
        (ORIGINAL, 'T_r is 290 K, see $T_r$ and `minT`.'),
        (ORIGINAL, 'T_r is 290 K.\n\nSee $T_r$ [^cite-x] and `minT`.'),
        (ORIGINAL, 'T_r is 290 K:\n- see $T_r$ [^cite-x] and `minT`.'),
        (NOTE, 'T_r is 290 K,\nsee $T_r$ [^cite-x] and `minT`.'),
        (ORIGINAL, 'T_r  is 290 K, see $T_r$ [^cite-x]\tand `minT`. '),
    ],
    ids=[
        'drops-everything',
        'number-changed',
        'math-dropped',
        'code-span-dropped',
        'footnote-dropped',
        'paragraph-break-added',
        'list-item-added',
        'line-leaves-the-note',
        'white-space-only',
    ],
)
def test_rewrite_that_loses_content_keeps_the_original_fix(original, rewrite):
    """The rewording call never sees the source code, so a rewording that changes a
    number, a formula, a code span, a footnote or the text structure could undo what
    the analysis checked; the original new_text is kept. A change of white space alone
    is not a rewording."""
    findings = [fixed('a', new=original)]
    assert analyze.apply_rewrites(findings, [{'id': 'a', 'new_text': rewrite}]) == 0
    assert findings[0]['suggested_fix']['new_text'] == original
    assert not analyze.keeps_content(original, rewrite) or rewrite.split() == original.split()


@pytest.mark.parametrize(
    'original, rewrite, indent',
    [
        (ORIGINAL, 'T_r is 290 K,\nsee $T_r$ [^cite-x] and `minT`.', ''),
        (ORIGINAL, 'T_r is 290 K,\n    see $T_r$ [^cite-x] and `minT`.', '    '),
        (NOTE, 'T_r is\n    290 K, see $T_r$\n    [^cite-x] and `minT`.', ''),
        ('T_r is\n    290 K, see $T_r$\n    [^cite-x] and `minT`.', NOTE, ''),
    ],
    ids=['line-added', 'line-added-in-note', 'note-rewrapped-longer', 'note-rewrapped-shorter'],
)
def test_rewrapping_onto_more_or_fewer_lines_is_allowed(original, rewrite, indent):
    """A slight rewording can wrap onto one line more or less; that is fine as long as
    each new line keeps an indentation the doc already uses there."""
    assert analyze.keeps_content(original, rewrite, indent)
    # Error path: the same rewrap with the indentation lost or changed is refused.
    assert not analyze.keeps_content(original, rewrite.replace('\n', '\n  '), indent)


def test_rewrite_that_reorders_kept_content_is_used():
    """Reordering and rephrasing are allowed as long as every number, formula, code
    span and footnote survives; a reply outside the schema is ignored, not raised."""
    findings = [fixed('a', new=ORIGINAL)]
    clearer = 'See $T_r$ [^cite-x] and `minT`: T_r is 290 K.'
    assert analyze.keeps_content(ORIGINAL, clearer)
    assert (
        analyze.apply_rewrites(findings, ['x', {'id': 'a'}, {'id': 'a', 'new_text': clearer}])
        == 1
    )
    assert findings[0]['suggested_fix']['new_text'] == clearer
    # Edge case: the same number twice must stay twice.
    assert not analyze.keeps_content('1 and 1', '1 and one')


def test_readability_prompt_carries_each_fix_with_its_context():
    """The prompt lists every fix with its quotes, so the rewording can stay true."""
    prompt = analyze.build_readability_prompt([fixed('a', old='alpha=1', new='alpha=5')])
    items = json.loads(prompt.split('The fixes, as JSON:', 1)[1])
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


@pytest.mark.parametrize(
    'readability', ['works', 'fails', 'bad-reply'], ids=['reworded', 'call-fails', 'bad-reply']
)
def test_main_writes_findings_whether_or_not_the_readability_call_works(
    tmp_path, monkeypatch, readability
):
    """The readability call only rewords: when it fails or replies outside its
    schema, the run still writes the findings with their original fixes instead of
    losing the analysis."""
    calls = []

    def fake_run(prompt, schema, timeout_s):
        calls.append(schema)
        if schema is analyze.FINDINGS_SCHEMA:
            return {'findings': [fixed('a'), {'id': 'gap', 'type': 'gap'}]}
        if readability == 'fails':
            raise RuntimeError('claude CLI did not finish within 300 s')
        if readability == 'bad-reply':
            return {'fixes': None}
        return {'fixes': [{'id': 'a', 'new_text': 'clearer words'}]}

    monkeypatch.setenv('CLAUDE_CODE_OAUTH_TOKEN', 'test')
    monkeypatch.setattr(analyze, 'REPO_ROOT', tmp_path)
    monkeypatch.setattr(analyze, 'build_prompt', lambda: 'prompt')
    monkeypatch.setattr(analyze, 'run_claude', fake_run)
    analyze.main()

    written = json.loads((tmp_path / 'analysis.json').read_text())
    expected = 'clearer words' if readability == 'works' else 'new words'
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


def cli_result(returncode=0, stdout=''):
    return SimpleNamespace(returncode=returncode, stdout=stdout, stderr='boom')


GOOD = {'is_error': False, 'subtype': 'success', 'structured_output': {'fixes': []}}


@pytest.mark.parametrize(
    'outcome, message',
    [
        (subprocess.TimeoutExpired('claude', 5), 'did not finish within 5 s'),
        (cli_result(returncode=2), 'exited 2'),
        (cli_result(stdout='not json'), 'did not print valid JSON'),
        (cli_result(stdout=json.dumps({**GOOD, 'is_error': True})), 'did not succeed'),
        (
            cli_result(stdout=json.dumps({**GOOD, 'structured_output': {}})),
            'no structured_output.fixes',
        ),
    ],
    ids=['timeout', 'non-zero-exit', 'non-json', 'is-error', 'no-structured-output'],
)
def test_run_claude_turns_every_cli_failure_into_runtime_error(monkeypatch, outcome, message):
    """Each way the CLI can fail ends as a RuntimeError naming the reason, which is the
    one error type main() handles for the analysis call."""

    def fake_run(*args, **kwargs):
        assert kwargs['timeout'] == 5 and kwargs['input'] == 'prompt'
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    monkeypatch.setattr(analyze.subprocess, 'run', fake_run)
    with pytest.raises(RuntimeError, match=message):
        analyze.run_claude('prompt', analyze.READABILITY_SCHEMA, 5)


def test_run_claude_returns_the_structured_output(monkeypatch):
    """A successful call returns structured_output as it is, and passes the schema to
    the CLI; an empty fix list (edge case) is still a valid reply."""
    seen = {}

    def fake_run(cmd, **kwargs):
        seen['schema'] = json.loads(cmd[cmd.index('--json-schema') + 1])
        return cli_result(stdout=json.dumps(GOOD))

    monkeypatch.setattr(analyze.subprocess, 'run', fake_run)
    assert analyze.run_claude('prompt', analyze.READABILITY_SCHEMA, 5) == {'fixes': []}
    assert seen['schema'] == analyze.READABILITY_SCHEMA


def test_one_line_fix_inside_a_note_may_wrap_only_at_the_note_indentation():
    """A one-line fix has no indentation of its own, so the indentation of the doc line
    it starts on decides where a new line may start: inside a note, at four spaces."""
    finding = fixed('a', old='old words', new=ORIGINAL)
    finding['doc_excerpt'] = '!!! note\n    The old words here.'
    unindented = 'T_r is 290 K,\nsee $T_r$ [^cite-x] and `minT`.'
    assert analyze.apply_rewrites([finding], [{'id': 'a', 'new_text': unindented}]) == 0
    indented = unindented.replace('\n', '\n    ')
    assert analyze.apply_rewrites([finding], [{'id': 'a', 'new_text': indented}]) == 1
    assert finding['suggested_fix']['new_text'] == indented
