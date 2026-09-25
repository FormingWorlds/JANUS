"""Tests for .github/scripts/docs_consistency/triage_and_apply.py.

Runs the triage step on hand-written analysis.json files in a temporary repo tree and
checks the contract the workflow relies on. See docs/How-to/test.md.
"""

import importlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

SCRIPT_DIR = Path(__file__).resolve().parents[1] / '.github' / 'scripts' / 'docs_consistency'
sys.path.insert(0, str(SCRIPT_DIR))
triage = importlib.import_module('triage_and_apply')

DOC_FILE = triage.DOC_FILES[0]
CODE_FILE = 'src/janus/utils/phys.py'

# 'alpha=1' occurs twice on purpose, so a first-occurrence replace would hit the dry
# section while the finding is about the moist one.
DOC_TEXT = (
    '## Dry\n'
    'The dry lapse uses alpha=1 here.\n'
    '## Moist\n'
    'The moist lapse uses alpha=1 too.\n'
    'Surface gamma=3 applies.\n'
)
CODE_TEXT = 'def f(x):\n    y = 2 * x\n    return y\n\n\ndef g(z):\n    return z + 1\n'

SAFE_FIX = {'old_text': 'moist lapse uses alpha=1', 'new_text': 'moist lapse uses alpha=5'}
NOTHING = {'has_fixes': 'false', 'has_inconsistency_items': 'false', 'has_gap_items': 'false'}


def finding(fid, **overrides):
    """A verified minor inconsistency on the moist passage, without a fix."""
    base = {
        'id': fid,
        'type': 'inconsistency',
        'severity': 'minor',
        'doc_file': DOC_FILE,
        'doc_excerpt': 'The moist lapse uses alpha=1 too.',
        'code_file': CODE_FILE,
        'code_lines': 'L2',
        'code_excerpt': '    y = 2 * x',
        'description': 'The doc and the code disagree.',
    }
    base.update(overrides)
    return base


def gap(fid, **overrides):
    """A verified gap with no related doc passage."""
    g = finding(fid, type='gap', doc_excerpt='', **overrides)
    del g['severity']
    return g


def outputs(**changed):
    """Expected workflow outputs: nothing to do, except the keys given."""
    return {**NOTHING, **{k: 'true' for k, v in changed.items() if v}}


@pytest.fixture
def repo(tmp_path, monkeypatch):
    """Temporary repo tree holding one reviewed doc and one reviewed source file."""
    (tmp_path / DOC_FILE).parent.mkdir(parents=True)
    (tmp_path / DOC_FILE).write_text(DOC_TEXT)
    (tmp_path / CODE_FILE).parent.mkdir(parents=True)
    (tmp_path / CODE_FILE).write_text(CODE_TEXT)
    # REPO_ROOT is a module-level constant computed at import, so patch it directly.
    monkeypatch.setattr(triage, 'REPO_ROOT', tmp_path)
    monkeypatch.setenv('GITHUB_OUTPUT', str(tmp_path / 'github_output'))
    monkeypatch.delenv('GITHUB_STEP_SUMMARY', raising=False)
    return tmp_path


def run_triage(repo, findings):
    """Run main() on the given findings and return its outputs and written files."""
    (repo / 'analysis.json').write_text(json.dumps(findings))
    (repo / 'github_output').write_text('')
    triage.main()
    lines = (repo / 'github_output').read_text().splitlines()
    return SimpleNamespace(
        outputs=dict(line.split('=', 1) for line in lines),
        doc=(repo / DOC_FILE).read_text(),
        inc=(repo / 'inconsistency_issue_body.md').read_text(),
        gaps=(repo / 'gap_issue_body.md').read_text(),
        pr=(repo / 'pr_body.md').read_text(),
    )


def post(repo, name, body):
    """Simulate the workflow posting body to the open issue: append it to its fetch file."""
    path = repo / name
    path.write_text((path.read_text() if path.exists() else '') + body)


def test_fixless_inconsistencies_open_no_pr_and_reach_the_inconsistency_issue(repo):
    """Minor inconsistencies with no applicable fix go to the inconsistency issue."""
    # Distinct code quotes: findings quoting the same code share a fingerprint and would
    # be merged into one issue entry.
    findings = [
        finding('no-fix'),
        finding(
            'fix-text-missing',
            code_excerpt='    return z + 1',
            suggested_fix={'old_text': 'beta=7', 'new_text': 'x'},
        ),
    ]
    r = run_triage(repo, findings)
    assert r.outputs == outputs(has_inconsistency_items=True)
    assert r.doc == DOC_TEXT
    assert '`no-fix`' in r.inc and '`fix-text-missing`' in r.inc
    assert '## Minor: needs a manual edit' in r.inc
    assert r.gaps == '' and 'No fixes applied' in r.pr


@pytest.mark.parametrize(
    ('old', 'new', 'reason'),
    [
        ('alpha=1', 'alpha=2', 'occurs 2 times'),
        ('gamma=3', 'gamma=3', 'identical'),
        ('gamma=3', 'gamma=4', 'does not lie inside doc_excerpt'),
        ('', 'x', 'old_text is empty'),
    ],
    ids=['old-text-not-unique', 'no-op-fix', 'fix-outside-excerpt', 'empty-old-text'],
)
def test_unsafe_fix_is_not_applied(repo, old, new, reason):
    """A fix that could edit the wrong passage, or nothing, leaves the doc untouched."""
    fix = {'suggested_fix': {'old_text': old, 'new_text': new}}
    r = run_triage(repo, [finding('unsafe', **fix)])
    assert r.doc == DOC_TEXT
    assert r.outputs == outputs(has_inconsistency_items=True)
    assert 'Suggested fix (not applied: ' in r.inc and reason in r.inc


def test_unique_fix_edits_only_the_quoted_passage(repo):
    """A unique old_text inside the excerpt is replaced; the other alpha=1 is kept."""
    r = run_triage(repo, [finding('good', suggested_fix=SAFE_FIX)])
    assert r.outputs == outputs(has_fixes=True)
    assert 'The moist lapse uses alpha=5 too.' in r.doc
    # The dry section also contains alpha=1 and must be left alone.
    assert 'The dry lapse uses alpha=1 here.' in r.doc
    assert '`good`' in r.pr and r.inc == ''


def test_second_fix_on_a_changed_passage_is_not_applied(repo):
    """Once one fix rewrites a passage, a second fix quoting the old text is refused."""
    first = finding(
        'first', suggested_fix={'old_text': 'alpha=1 too', 'new_text': 'alpha=5 too'}
    )
    second = finding(
        'second',
        code_excerpt='    return y',
        suggested_fix={'old_text': 'moist lapse', 'new_text': 'wet lapse'},
    )
    r = run_triage(repo, [first, second])
    assert 'alpha=5 too' in r.doc and 'wet lapse' not in r.doc
    assert r.outputs == outputs(has_fixes=True, has_inconsistency_items=True)
    assert 'changed by an earlier fix' in r.inc


def test_serious_inconsistency_is_filed_first_and_never_applied(repo):
    """A serious finding with a safe fix is blocked and heads the inconsistency issue.

    A serious disagreement may be a bug in the code; rewriting the doc to match the code
    would hide it, so a human decides which side is wrong.
    """
    r = run_triage(
        repo,
        [
            finding('small', code_excerpt='    return y'),
            finding('code-bug', severity='serious', suggested_fix=SAFE_FIX),
        ],
    )
    assert r.doc == DOC_TEXT
    assert r.outputs == outputs(has_inconsistency_items=True)
    assert 'not a minor finding' in r.inc and 'No fixes applied' in r.pr
    # The serious section comes before the minor one, whatever the input order.
    serious_at = r.inc.index('## Serious: possible code bug')
    assert serious_at < r.inc.index('`code-bug`') < r.inc.index('## Minor: needs a manual edit')
    assert r.inc.index('## Minor: needs a manual edit') < r.inc.index('`small`')


def test_unrated_inconsistency_is_blocked_but_not_labelled_serious(repo):
    """A finding without a severity is not auto-fixed and is filed under minor."""
    f = finding('unrated', suggested_fix=SAFE_FIX)
    del f['severity']
    r = run_triage(repo, [f])
    assert r.doc == DOC_TEXT
    assert '## Serious' not in r.inc and '## Minor: needs a manual edit' in r.inc
    assert 'not a minor finding' in r.inc and '`unrated`' in r.inc


@pytest.mark.parametrize(
    'overrides',
    [
        {'code_excerpt': '    y = 3 * x'},
        {'doc_excerpt': 'The model uses beta=7.'},
        {'code_file': '../../etc/passwd'},
        {'doc_excerpt': ''},
        {'code_excerpt': '    y = 3 * x', 'severity': 'serious'},
    ],
    ids=[
        'invented-code-quote',
        'invented-doc-quote',
        'code-file-outside-list',
        'no-doc-quote',
        'serious-but-misquoted',
    ],
)
def test_unmatched_quotes_are_filed_as_unverified_and_never_applied(repo, overrides):
    """An inconsistency whose quotes are not in the files is unverified, never fixed."""
    r = run_triage(repo, [finding('bad', suggested_fix=SAFE_FIX, **overrides)])
    assert r.doc == DOC_TEXT
    assert r.outputs == outputs(has_inconsistency_items=True)
    assert '## Unverified inconsistencies' in r.inc and '## Serious' not in r.inc
    assert 'not applied: finding unverified' in r.inc and r.gaps == ''


def test_findings_are_routed_to_the_issue_matching_their_type(repo):
    """Gaps, verified or not, go to the gap issue; inconsistencies never do."""
    r = run_triage(
        repo,
        [
            gap('real-gap', code_excerpt='    return z + 1'),
            gap('invented-gap', code_excerpt='no_such_call()'),
            finding('inc'),
        ],
    )
    assert r.outputs == outputs(has_inconsistency_items=True, has_gap_items=True)
    assert '## Documentation gaps' in r.gaps and '## Unverified gaps' in r.gaps
    assert '`real-gap`' in r.gaps and '`invented-gap`' in r.gaps and '`inc`' not in r.gaps
    assert '`inc`' in r.inc and 'gap`' not in r.inc


def test_line_number_prefixes_are_stripped_before_matching(repo):
    """A code quote copied with 'N: ' prefixes verifies, and is posted without them."""
    quoted = '2:     y = 2 * x\n3:     return y'
    r = run_triage(repo, [finding('prefixed', code_excerpt=quoted)])
    assert '## Unverified' not in r.inc
    assert '    y = 2 * x\n    return y' in r.inc
    assert '2:     y' not in r.inc


def test_gap_may_have_an_empty_doc_excerpt(repo):
    """A gap with no related doc passage is filed as a gap, not as unverified."""
    r = run_triage(repo, [gap('gap')])
    assert r.outputs == outputs(has_gap_items=True)
    assert '## Documentation gaps' in r.gaps and '## Unverified' not in r.gaps
    assert '_no related passage_' in r.gaps


def test_already_reported_findings_are_not_posted_again(repo):
    """Findings already in the open issues are skipped, even under new ids."""
    r = run_triage(repo, [finding('a'), gap('b', code_excerpt='    return z + 1')])
    assert (
        r.inc.count('docs-consistency-fp:') == 1 and r.gaps.count('docs-consistency-fp:') == 1
    )
    post(repo, 'existing_inconsistency_issue.md', r.inc)
    post(repo, 'existing_gap_issue.md', r.gaps)

    # Same content, renamed ids and reformatted quote, plus one new gap.
    second_run = [
        finding('a-renamed', code_excerpt='2:     y = 2 * x\n'),
        gap('b-renamed', code_excerpt='    return z + 1'),
        gap('new', code_excerpt='def g(z):'),
    ]
    r = run_triage(repo, second_run)
    assert r.outputs == outputs(has_gap_items=True)
    assert '`new`' in r.gaps and '`b-renamed`' not in r.gaps and r.inc == ''
    assert '1 other finding(s) from this run were already reported' in r.gaps

    # Nothing new at all: no issue post.
    post(repo, 'existing_gap_issue.md', r.gaps)
    r = run_triage(repo, second_run)
    assert r.outputs == NOTHING
    assert r.inc == '' and r.gaps == ''


def test_finding_is_reposted_only_when_its_severity_rises(repo):
    """Minor then serious is reposted as an escalation; serious then minor is not."""
    r = run_triage(repo, [finding('first-minor')])
    assert '## Minor: needs a manual edit' in r.inc
    post(repo, 'existing_inconsistency_issue.md', r.inc)

    # Same code, now rated serious: reposted under the serious section.
    r = run_triage(repo, [finding('now-serious', severity='serious')])
    assert r.outputs == outputs(has_inconsistency_items=True)
    assert '## Serious: possible code bug' in r.inc and '`now-serious`' in r.inc
    assert '1 finding(s) above were reported before as unverified or at a lower' in r.inc
    post(repo, 'existing_inconsistency_issue.md', r.inc)

    # Serious again under a new id, or back to minor: not posted.
    for severity in ('serious', 'minor'):
        r = run_triage(repo, [finding('again', severity=severity)])
        assert r.outputs == NOTHING, severity
        assert r.inc == ''


def test_unverified_copy_does_not_hide_the_verified_finding(repo):
    """A finding first posted as unverified is reposted once a later run verifies it."""
    suspect = {'severity': 'serious', 'doc_excerpt': 'An invented doc quote.'}
    r = run_triage(repo, [finding('suspect', **suspect)])
    assert '## Unverified inconsistencies' in r.inc and 'serious-unverified -->' in r.inc
    post(repo, 'existing_inconsistency_issue.md', r.inc)

    # Same code quote, now verified: posted under the serious section.
    r = run_triage(repo, [finding('verified', severity='serious')])
    assert '## Serious: possible code bug' in r.inc and '`verified`' in r.inc
    post(repo, 'existing_inconsistency_issue.md', r.inc)

    # A suspect copy after the verified one is not posted again.
    r = run_triage(repo, [finding('suspect-again', **suspect)])
    assert r.outputs == NOTHING and r.inc == ''


def test_duplicate_within_one_run_is_posted_once(repo):
    """Two findings quoting the same code in one run produce a single issue entry."""
    r = run_triage(repo, [finding('one'), finding('two')])
    assert r.inc.count('docs-consistency-fp:') == 1
    assert 'already reported' not in r.inc


def test_empty_analysis_sets_no_outputs(repo):
    """No findings: no PR, no issue posts, and empty issue bodies."""
    r = run_triage(repo, [])
    assert r.outputs == NOTHING
    assert r.doc == DOC_TEXT
    assert r.inc == '' and r.gaps == ''
