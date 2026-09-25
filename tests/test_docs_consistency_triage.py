"""Tests for .github/scripts/docs_consistency/triage_and_apply.py.

Runs the triage step on hand-written analysis.json files in a temporary repo tree and
checks the contract the workflow relies on. See docs/How-to/test.md.
"""

import importlib
import json
import sys
from pathlib import Path

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
    """Run main() on the given findings; return (outputs, doc text, issue body, PR body)."""
    (repo / 'analysis.json').write_text(json.dumps(findings))
    (repo / 'github_output').write_text('')
    triage.main()
    lines = (repo / 'github_output').read_text().splitlines()
    outputs = dict(line.split('=', 1) for line in lines)
    return (
        outputs,
        (repo / DOC_FILE).read_text(),
        (repo / 'gap_issue_body.md').read_text(),
        (repo / 'pr_body.md').read_text(),
    )


def test_fixless_inconsistencies_open_no_pr_and_reach_the_issue(repo):
    """Inconsistencies with no applicable fix leave has_fixes false and go to the issue."""
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
    outputs, doc, issue, pr = run_triage(repo, findings)
    assert outputs == {'has_fixes': 'false', 'has_issue_items': 'true'}
    assert doc == DOC_TEXT
    assert '`no-fix`' in issue and '`fix-text-missing`' in issue
    assert '## Inconsistencies needing a manual edit' in issue
    assert 'No fixes applied' in pr


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
    outputs, doc, issue, _ = run_triage(repo, [finding('unsafe', **fix)])
    assert doc == DOC_TEXT
    assert outputs['has_fixes'] == 'false'
    assert 'Suggested fix (not applied: ' in issue and reason in issue


def test_unique_fix_edits_only_the_quoted_passage(repo):
    """A unique old_text inside the excerpt is replaced; the other alpha=1 is kept."""
    fix = {'old_text': 'moist lapse uses alpha=1', 'new_text': 'moist lapse uses alpha=5'}
    outputs, doc, issue, pr = run_triage(repo, [finding('good', suggested_fix=fix)])
    assert outputs == {'has_fixes': 'true', 'has_issue_items': 'false'}
    assert 'The moist lapse uses alpha=5 too.' in doc
    # The dry section also contains alpha=1 and must be left alone.
    assert 'The dry lapse uses alpha=1 here.' in doc
    assert '`good`' in pr and issue == ''


def test_second_fix_on_a_changed_passage_is_not_applied(repo):
    """Once one fix rewrites a passage, a second fix quoting the old text is refused."""
    first = finding(
        'first',
        severity='serious',
        suggested_fix={'old_text': 'alpha=1 too', 'new_text': 'alpha=5 too'},
    )
    second = finding(
        'second', suggested_fix={'old_text': 'moist lapse', 'new_text': 'wet lapse'}
    )
    outputs, doc, issue, _ = run_triage(repo, [second, first])
    # Serious findings are applied first, whatever their order in analysis.json.
    assert 'alpha=5 too' in doc and 'wet lapse' not in doc
    assert outputs['has_fixes'] == 'true'
    assert 'changed by an earlier fix' in issue


@pytest.mark.parametrize(
    'overrides',
    [
        {'code_excerpt': '    y = 3 * x'},
        {'doc_excerpt': 'The model uses beta=7.'},
        {'code_file': '../../etc/passwd'},
        {'doc_excerpt': ''},
    ],
    ids=['invented-code-quote', 'invented-doc-quote', 'code-file-outside-list', 'no-doc-quote'],
)
def test_unmatched_quotes_are_filed_as_unverified_and_never_applied(repo, overrides):
    """A finding whose quotes are not in the files is posted as unverified, never fixed."""
    fix = {'old_text': 'moist lapse uses alpha=1', 'new_text': 'moist lapse uses alpha=5'}
    outputs, doc, issue, _ = run_triage(repo, [finding('bad', suggested_fix=fix, **overrides)])
    assert doc == DOC_TEXT
    assert outputs == {'has_fixes': 'false', 'has_issue_items': 'true'}
    assert '## Unverified findings' in issue
    assert 'not applied: finding unverified' in issue


def test_line_number_prefixes_are_stripped_before_matching(repo):
    """A code quote copied with 'N: ' prefixes verifies, and is posted without them."""
    quoted = '2:     y = 2 * x\n3:     return y'
    _, _, issue, _ = run_triage(repo, [finding('prefixed', code_excerpt=quoted)])
    assert '## Unverified findings' not in issue
    assert '    y = 2 * x\n    return y' in issue
    assert '2:     y' not in issue


def test_gap_may_have_an_empty_doc_excerpt(repo):
    """A gap with no related doc passage is filed as a gap, not as unverified."""
    gap = finding('gap', type='gap', doc_excerpt='')
    del gap['severity']
    _, _, issue, _ = run_triage(repo, [gap])
    assert '## Documentation gaps' in issue
    assert '## Unverified findings' not in issue
    assert '_no related passage_' in issue


def test_already_reported_findings_are_not_posted_again(repo):
    """Findings already in the open issue are skipped, even under new ids."""
    first_run = [finding('a'), finding('b', type='gap', code_excerpt='    return z + 1')]
    _, _, issue, _ = run_triage(repo, first_run)
    assert issue.count('docs-consistency-fp:') == 2
    (repo / 'existing_issue.md').write_text(issue)

    # Same content, renamed ids and reformatted quote, plus one new finding.
    second_run = [
        finding('a-renamed', code_excerpt='2:     y = 2 * x\n'),
        finding('b-renamed', type='gap', code_excerpt='    return z + 1'),
        finding('new', type='gap', code_excerpt='def g(z):'),
    ]
    outputs, _, issue, _ = run_triage(repo, second_run)
    assert outputs['has_issue_items'] == 'true'
    assert '`new`' in issue
    assert '`a-renamed`' not in issue and '`b-renamed`' not in issue
    assert '2 other finding(s) from this run were already reported' in issue

    # Nothing new at all: no issue post.
    (repo / 'existing_issue.md').write_text((repo / 'existing_issue.md').read_text() + issue)
    outputs, _, issue, _ = run_triage(repo, second_run)
    assert outputs == {'has_fixes': 'false', 'has_issue_items': 'false'}
    assert issue == ''


def test_duplicate_within_one_run_is_posted_once(repo):
    """Two findings quoting the same code in one run produce a single issue entry."""
    _, _, issue, _ = run_triage(repo, [finding('one'), finding('two')])
    assert issue.count('docs-consistency-fp:') == 1
    assert 'already reported' not in issue


def test_empty_analysis_sets_no_outputs(repo):
    """No findings: no PR, no issue post, and an empty issue body."""
    outputs, doc, issue, _ = run_triage(repo, [])
    assert outputs == {'has_fixes': 'false', 'has_issue_items': 'false'}
    assert doc == DOC_TEXT
    assert issue == ''
