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
NOTHING = {
    'has_fixes': 'false',
    'has_inconsistency_items': 'false',
    'has_gap_items': 'false',
    'has_new_serious': 'false',
}


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
        applied=[f['id'] for f in json.loads((repo / 'applied_fixes.json').read_text())],
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


# A hand-wrapped note: each sentence runs over several lines.
# A hand-wrapped note: sentences run over several lines, and lines start with math, a
# capitalised name or a lowercase word.
WRAPPED_DOC = (
    '!!! note\n'
    '    Setting T_r = 0 K disables the stratosphere: since no level\n'
    '    reaches 0 K, the criterion never triggers and\n'
    '    $T_r$ is not used by\n'
    '    SOCRATES at all. The step is skipped.\n'
    '    Nothing else runs.\n'
    '\n'
    '## Head\n'
    '| a | b |\n'
)


@pytest.mark.parametrize(
    ('old', 'reason'),
    [
        ('Setting T_r = 0 K disables the stratosphere: since no level', 'stops mid-sentence'),
        ('reaches 0 K, the criterion never triggers and', 'stops mid-sentence'),
        ('$T_r$ is not used by', 'stops mid-sentence'),
        ('$T_r$ is not used by\n', 'stops mid-sentence'),
        ('reaches 0 K, the criterion never triggers and\n    $T_r$', 'starts mid-sentence'),
        ('    SOCRATES at all.', 'starts mid-sentence'),
        ('\n    SOCRATES at all.', 'starts mid-sentence'),
    ],
    ids=[
        'stops-before-lowercase-line',
        'stops-before-math-line',
        'stops-before-capitalised-line',
        'stops-on-trailing-line-break',
        'starts-after-unfinished-line',
        'starts-with-indentation',
        'starts-with-leading-line-break',
    ],
)
def test_fix_on_part_of_a_wrapped_sentence_is_not_applied(repo, old, reason):
    """Editing some physical lines of a hand-wrapped sentence leaves the rest unchecked,
    whatever the next line starts with and however the quote is padded."""
    (repo / DOC_FILE).write_text(WRAPPED_DOC)
    excerpt = WRAPPED_DOC.split('\n', 1)[1].rstrip('\n')
    fix = {'old_text': old, 'new_text': old.upper()}
    r = run_triage(repo, [finding('half', doc_excerpt=excerpt, suggested_fix=fix)])
    assert r.doc == WRAPPED_DOC
    assert r.outputs == outputs(has_inconsistency_items=True)
    assert reason in r.inc


@pytest.mark.parametrize(
    'old',
    [
        'Setting T_r = 0 K disables the stratosphere: since no level\n'
        '    reaches 0 K, the criterion never triggers and\n'
        '    $T_r$ is not used by\n'
        '    SOCRATES at all.',
        'The step is skipped.',
        'The step is skipped.\n',
        'Nothing else runs.',
        'the criterion never triggers',
        '## Head',
        '| a | b |',
    ],
    ids=[
        'whole-wrapped-sentence',
        'sentence-ending-at-a-line-break',
        'same-with-trailing-line-break',
        'sentence-ending-before-blank-line',
        'clause-inside-one-line',
        'heading-before-table',
        'table-row-after-heading',
    ],
)
def test_fix_that_cuts_no_sentence_is_applied(repo, old):
    """Whole sentences, edits inside one line and whole block lines are still applied."""
    (repo / DOC_FILE).write_text(WRAPPED_DOC)
    excerpt = WRAPPED_DOC.split('\n', 1)[1].rstrip('\n')
    fix = {'old_text': old, 'new_text': old.upper()}
    r = run_triage(repo, [finding('whole', doc_excerpt=excerpt, suggested_fix=fix)])
    assert r.outputs == outputs(has_fixes=True)
    assert r.doc == WRAPPED_DOC.replace(old, old.upper())


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
    assert r.outputs == outputs(has_inconsistency_items=True, has_new_serious=True)
    assert 'not a minor finding' in r.inc and 'No fixes applied' in r.pr
    # The serious section comes before the minor one, whatever the input order.
    serious_at = r.inc.index('## Serious: possible code bug')
    assert serious_at < r.inc.index('`code-bug`') < r.inc.index('## Minor: needs a manual edit')
    assert r.inc.index('## Minor: needs a manual edit') < r.inc.index('`small`')


def test_serious_fix_is_listed_in_the_pr_but_not_applied(repo):
    """The PR lists serious findings' doc fixes unapplied; fixless ones stay in the issue."""
    serious_fix = {'old_text': 'gamma=3', 'new_text': 'gamma=4'}
    r = run_triage(
        repo,
        [
            finding('small', suggested_fix=SAFE_FIX),
            finding(
                'code-bug',
                severity='serious',
                doc_excerpt='Surface gamma=3 applies.',
                code_excerpt='    return y',
                suggested_fix=serious_fix,
            ),
            finding(
                'no-fix-bug',
                severity='serious',
                doc_excerpt='Surface gamma=3 applies.',
                code_excerpt='    return z + 1',
            ),
        ],
    )
    assert r.outputs == outputs(
        has_fixes=True, has_inconsistency_items=True, has_new_serious=True
    )
    # Only the minor fix is applied.
    assert 'alpha=5 too' in r.doc and 'Surface gamma=3 applies.' in r.doc
    section = r.pr.index('## Serious findings: suggested doc-side fixes, not applied')
    assert r.pr.index('`small`') < section < r.pr.index('`code-bug`')
    assert 'Suggested fix (not applied: not a minor finding' in r.pr
    assert '`no-fix-bug`' not in r.pr and '`no-fix-bug`' in r.inc
    # Only the applied fix is handed on to a later update of an open PR.
    assert r.applied == ['small']


def test_serious_fix_alone_writes_no_pr_section(repo):
    """With nothing applied no PR is opened, so the serious section is not written."""
    r = run_triage(repo, [finding('code-bug', severity='serious', suggested_fix=SAFE_FIX)])
    assert r.outputs == outputs(has_inconsistency_items=True, has_new_serious=True)
    assert 'No fixes applied' in r.pr and '## Serious findings' not in r.pr
    assert r.applied == []


def test_quote_cut_from_an_equation_is_shown_as_the_whole_equation(repo):
    """A quote partly inside $...$ or $$...$$ is widened so it renders, with the exact
    quote kept underneath; a quote outside math is shown unchanged."""
    (repo / DOC_FILE).write_text('Intro $x$ text.\n$$a + b \\tag{5}$$\nPlain words here.\n')
    r = run_triage(
        repo,
        [
            finding('display', doc_excerpt='a + b'),
            finding('inline', doc_excerpt='Intro $x', code_excerpt='    return y'),
            finding('plain', doc_excerpt='Plain words', code_excerpt='    return z + 1'),
        ],
    )
    assert '> $$a + b \\tag{5}$$' in r.inc and '_Quoted part:_ ` a + b `' in r.inc
    assert '> Intro $x$\n' in r.inc and '_Quoted part:_ ` Intro $x `' in r.inc
    assert '> Plain words\n' in r.inc and r.inc.count('_Quoted part:_') == 2


def test_dashes_between_words_are_replaced_but_doc_dashes_are_kept(repo):
    """Dashes the model writes between words are replaced in the description and the fix.
    A standalone dash such as a table's empty-cell marker is kept, and so is a doc quote."""
    em, en = '\u2014', '\u2013'
    row = f'| Overlap | {em} | {en} | Random |'
    (repo / DOC_FILE).write_text(DOC_TEXT + f'Note {em} kept.\n{row}\n')
    dashed = finding(
        'dashed',
        description=f'Doc and code disagree {em} badly {en} over 10{en}20 K.',
        suggested_fix={
            'old_text': 'moist lapse',
            'new_text': f'moist {em} wet, 10{en}20 K, lapse',
        },
    )
    table = finding(
        'table',
        doc_excerpt=row,
        code_excerpt='    return z + 1',
        suggested_fix={'old_text': f'{en} | Random |', 'new_text': f'{en} | Equivalent |'},
    )
    quote = finding(
        'quote', severity='serious', doc_excerpt=f'Note {em} kept.', code_excerpt='    return y'
    )
    r = run_triage(repo, [dashed, table, quote])
    assert 'moist, wet, 10-20 K, lapse uses alpha=1 too.' in r.doc
    # The table's placeholder survives the fix next to it.
    assert f'| Overlap | {em} | {en} | Equivalent |' in r.doc
    assert 'Doc and code disagree, badly, over 10-20 K.' in r.pr
    # The en-dash left in the PR is the table's own placeholder, quoted verbatim.
    assert r.pr.count(en) == r.pr.count(f'{en} |')
    # The only dash in the issue is the verbatim doc quote.
    assert f'> Note {em} kept.' in r.inc and r.inc.count(em) == 1


def onto(repo, branch_doc):
    """Simulate the open PR's checkout holding branch_doc and re-apply this run's fixes."""
    target = repo / 'pr'
    (target / DOC_FILE).parent.mkdir(parents=True)
    (target / DOC_FILE).write_text(branch_doc)
    (repo / 'github_output').write_text('')
    triage.apply_onto(target)
    return SimpleNamespace(
        doc=(target / DOC_FILE).read_text(),
        update=(repo / 'pr_update.md').read_text(),
        inc=(repo / 'inconsistency_issue_body.md').read_text(),
        output=(repo / 'github_output').read_text(),
    )


GAMMA_FIX = finding(
    'gamma',
    doc_excerpt='Surface gamma=3 applies.',
    code_excerpt='    return z + 1',
    suggested_fix={'old_text': 'gamma=3', 'new_text': 'gamma=4'},
)


def test_fixes_go_onto_an_untouched_open_pr_one_by_one(repo):
    """On a PR branch that has none of this run's fixes, each one is applied and listed."""
    run_triage(repo, [finding('alpha', suggested_fix=SAFE_FIX), GAMMA_FIX])
    r = onto(repo, DOC_TEXT)
    assert 'moist lapse uses alpha=5' in r.doc and 'Surface gamma=4 applies.' in r.doc
    assert '`alpha`' in r.update and '`gamma`' in r.update
    assert 'conflict' not in r.inc and r.output == ''


def test_fix_the_pr_has_is_skipped_and_a_corrected_one_is_filed_not_dropped(repo):
    """A reviewer corrected one fix on the PR and the other is already there: the first
    is filed in the issue, the second is not repeated, the PR text is untouched."""
    run_triage(repo, [finding('alpha', suggested_fix=SAFE_FIX), GAMMA_FIX])
    branch = DOC_TEXT.replace('alpha=1 too', 'alpha=6 too').replace('gamma=3', 'gamma=4')
    r = onto(repo, branch)
    assert r.doc == branch and r.update == ''
    assert '## Fixes that conflict with the open PR' in r.inc and '`alpha`' in r.inc
    assert '`gamma`' not in r.inc and r.output == 'has_inconsistency_items=true\n'


def test_fix_next_to_a_reviewer_edit_counts_as_a_conflict(repo):
    """old_text is still on the PR but the quoted sentence around it was edited: the fix
    is filed, not applied, as on main where the whole quote must still match."""
    run_triage(repo, [GAMMA_FIX])
    branch = DOC_TEXT.replace('applies', 'holds')
    r = onto(repo, branch)
    assert r.doc == branch
    assert '`gamma`' in r.inc and r.update == ''


def test_fix_that_extends_its_passage_is_not_applied_twice(repo):
    """new_text containing old_text is recognised as already in the PR."""
    grow = {'old_text': 'gamma=3 applies.', 'new_text': 'gamma=3 applies. Always.'}
    run_triage(
        repo, [finding('grow', doc_excerpt='Surface gamma=3 applies.', suggested_fix=grow)]
    )
    r = onto(repo, DOC_TEXT.replace('applies.', 'applies. Always.'))
    assert r.doc.count('Always.') == 1
    assert r.update == '' and r.output == ''


def test_fix_the_pr_offered_before_is_not_re_added_after_a_revert(repo):
    """A reviewer reverted a fix on the PR: its fingerprint is in the PR description, so
    the next run neither re-adds it nor files it, and the revert stands."""
    r = run_triage(repo, [GAMMA_FIX])
    assert '<!-- docs-consistency-fp: ' in r.pr
    # Control: without the PR description the reverted passage would get the fix again.
    assert 'Surface gamma=4 applies.' in onto(repo, DOC_TEXT).doc
    post(repo, 'existing_pr_body.md', r.pr)
    (repo / 'pr' / DOC_FILE).unlink()
    (repo / 'pr' / DOC_FILE).parent.rmdir()
    reverted = onto(repo, DOC_TEXT)
    assert reverted.doc == DOC_TEXT and reverted.update == ''
    assert 'gamma' not in reverted.inc and reverted.output == ''


def test_conflict_already_in_the_issue_is_not_filed_again(repo):
    """A conflicting fix reported in an earlier run is not posted a second time."""
    run_triage(repo, [GAMMA_FIX])
    branch = DOC_TEXT.replace('applies', 'holds')
    first = onto(repo, branch)
    assert '`gamma`' in first.inc
    post(repo, 'existing_inconsistency_issue.md', first.inc)
    (repo / 'pr' / DOC_FILE).unlink()
    (repo / 'pr' / DOC_FILE).parent.rmdir()
    # The next run starts from main again, where the fix is still unapplied.
    (repo / DOC_FILE).write_text(DOC_TEXT)
    run_triage(repo, [GAMMA_FIX])
    r = onto(repo, branch)
    assert 'Fixes that conflict' not in r.inc
    assert r.output == ''


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


@pytest.mark.parametrize(
    ('excerpt', 'verified'),
    [
        ('    y = 2 * x\n...\n    return z + 1', True),
        ('    return z + 1\n...\n    y = 2 * x', False),
        ('    y = 2 * x\n    return z + 1', False),
        ('...', False),
        ('...\n    y = 2 * x\n...', True),
        ('    y = 2 * x\n...\n...\n    return z + 1', True),
        ('    y = 2 * x\n...\n    y = 3 * x', False),
    ],
    ids=[
        'blocks-in-order',
        'blocks-out-of-order',
        'joined-without-marker',
        'marker-only',
        'markers-at-both-ends',
        'doubled-marker',
        'second-block-invented',
    ],
)
def test_code_blocks_joined_by_an_ellipsis_line_are_checked_one_by_one(repo, excerpt, verified):
    """Separate code passages joined by a '...' line verify only if each one occurs
    verbatim, in order; lines joined without the marker do not."""
    r = run_triage(repo, [finding('joined', code_excerpt=excerpt)])
    assert ('## Unverified inconsistencies' not in r.inc) is verified
    assert ('## Minor: needs a manual edit' in r.inc) is verified


def test_line_number_prefixes_are_stripped_before_matching(repo):
    """A code quote copied with 'N: ' prefixes verifies, and is posted without them."""
    quoted = '2:     y = 2 * x\n3:     return y'
    r = run_triage(repo, [finding('prefixed', code_excerpt=quoted)])
    assert '## Unverified' not in r.inc
    assert '    y = 2 * x\n    return y' in r.inc
    assert '2:     y' not in r.inc


def test_multi_line_fix_is_shown_as_a_diff_block(repo):
    """A fix spanning lines is rendered as a fenced diff, so a blank line or a line that
    starts a list cannot break out of it; a one-line fix keeps its inline form."""
    doc = DOC_TEXT + '!!! note\n    First line here\n    ends the sentence.\n'
    (repo / DOC_FILE).write_text(doc)
    old = 'First line here\n    ends the sentence.'
    new = 'First line here\n\n    - ends the `cp` sentence.'
    wrapped = finding(
        'wrapped',
        doc_excerpt=old,
        code_excerpt='    return y',
        suggested_fix={'old_text': old, 'new_text': new},
    )
    r = run_triage(repo, [wrapped, finding('short', suggested_fix=SAFE_FIX)])
    block = '```diff\n- First line here\n-     ends the sentence.\n+ First line here\n+ \n'
    assert block + '+     - ends the `cp` sentence.\n```' in r.pr
    assert '- old: ` moist lapse uses alpha=1 `' in r.pr


def test_gap_may_have_an_empty_doc_excerpt(repo):
    """A gap with no related doc passage is filed as a gap, not as unverified."""
    r = run_triage(repo, [gap('gap')])
    assert r.outputs == outputs(has_gap_items=True)
    assert '## Documentation gaps' in r.gaps and '## Unverified' not in r.gaps
    assert '_no related passage_' in r.gaps


def test_gap_without_a_doc_file_is_not_filed_as_unverified(repo):
    """A gap with no related passage may leave doc_file empty; a quote still needs one."""
    r = run_triage(
        repo,
        [
            gap('no-file', doc_file=''),
            finding('bad-file', doc_file='', code_excerpt='    return y'),
        ],
    )
    assert '`no-file`' in r.gaps and '## Unverified' not in r.gaps
    assert '**Doc**: _no related passage_' in r.gaps
    assert '## Unverified inconsistencies' in r.inc and '`bad-file`' in r.inc


def test_backticks_in_a_fix_do_not_break_its_code_span(repo):
    """Fix text containing backticks is wrapped in a longer backtick delimiter."""
    fix = {'old_text': 'moist lapse', 'new_text': 'lapse from `cp`'}
    r = run_triage(repo, [finding('ticks', suggested_fix=fix)])
    assert 'The lapse from `cp` uses alpha=1 too.' in r.doc
    assert '- new: `` lapse from `cp` ``' in r.pr and '- old: ` moist lapse `' in r.pr


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
    assert r.outputs == outputs(has_inconsistency_items=True, has_new_serious=True)
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
