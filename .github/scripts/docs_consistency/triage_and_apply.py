#!/usr/bin/env python3
"""Deterministic triage of analysis.json: verify quoted excerpts, apply literal doc fixes,
and write pr_body.md, gap_issue_body.md and serious_issue_body.md for the workflow to use.
Serious inconsistencies get their own issue, since they may point to a bug in the code.
No LLM calls: everything is a plain string operation.
"""

import hashlib
import json
import os
import re
import sys
from pathlib import Path

from analyze import DOC_FILES, SOURCE_FILES

REPO_ROOT = Path(__file__).resolve().parents[3]

SEVERITY_ORDER = {'serious': 0, 'minor': 1}
ALLOWED_DOC_FILES = set(DOC_FILES)
ALLOWED_SOURCE_FILES = set(SOURCE_FILES)

# The source listing in the prompt is prefixed with 'N: ' line numbers, which the model
# may copy into code_excerpt. Strip.
LINE_PREFIX = re.compile(r'^\d+: ', flags=re.MULTILINE)

# Hidden marker embedded in every finding posted to the issue, so later runs can tell which
# findings were already reported.
FINGERPRINT_MARKER = '<!-- docs-consistency-fp: {} -->'
FINGERPRINT_RE = re.compile(r'<!-- docs-consistency-fp: ([0-9a-f]{16}) -->')


def load_findings():
    return json.loads((REPO_ROOT / 'analysis.json').read_text())


def fingerprint(finding):
    """Hash of code_file plus the normalised code_excerpt, stable across runs.

    Whitespace within each line is collapsed and blank lines dropped, so re-indenting the
    quote or adding a trailing newline does not change the hash. Call after
    verify_excerpts, which strips line-number prefixes from code_excerpt.
    """
    lines = (' '.join(line.split()) for line in finding['code_excerpt'].splitlines())
    normalised = '\n'.join(line for line in lines if line)
    key = f'{finding.get("code_file")}\n{normalised}'
    return hashlib.sha256(key.encode()).hexdigest()[:16]


def load_reported_fingerprints(existing_file):
    """Fingerprints in an open issue's body and comments, as fetched by the workflow."""
    path = REPO_ROOT / existing_file
    return set(FINGERPRINT_RE.findall(path.read_text())) if path.exists() else set()


def drop_reported(sections, existing_file):
    """Filter each section's (finding, detail) pairs down to those not yet in the issue.

    Sections are processed in order and share one set, so a finding repeated within this
    run is kept only in the first section it appears in. Returns the filtered sections,
    the number already in the issue, and the number of within-run duplicates.
    """
    previously = load_reported_fingerprints(existing_file)
    reported = set(previously)
    filtered, already, duplicates = [], 0, 0
    for items in sections:
        kept = []
        for f, detail in items:
            fp = f['fingerprint']
            if fp in previously:
                already += 1
            elif fp in reported:
                duplicates += 1
            else:
                reported.add(fp)
                kept.append((f, detail))
        filtered.append(kept)
    return filtered, already, duplicates


def issue_body(sections, already_reported):
    """Join (heading, intro, rendered findings) sections, skipping empty ones."""
    lines = []
    for heading, intro, rendered in sections:
        if rendered:
            lines += [heading, intro, *rendered]
    if lines and already_reported:
        lines.append(
            f'_{already_reported} other finding(s) from this run were already reported in '
            'this issue and are not repeated._\n'
        )
    return '\n'.join(lines)


def verify_excerpts(finding):
    """Check that the quoted doc and code passages exist verbatim in the named files.

    Normalises code_excerpt in place (line-number prefixes stripped) so the posted quote
    is the verified one. Returns a list of problems; empty means verified.
    """
    problems = []
    doc_file = finding.get('doc_file')
    doc_excerpt = finding.get('doc_excerpt') or ''
    if doc_file not in ALLOWED_DOC_FILES:
        problems.append(f'doc_file {doc_file!r} is not one of the reviewed doc files')
    elif doc_excerpt:
        if doc_excerpt not in (REPO_ROOT / doc_file).read_text():
            problems.append(f'doc_excerpt not found verbatim in {doc_file}')
    elif finding.get('type') != 'gap':
        problems.append('doc_excerpt is empty')

    code_file = finding.get('code_file')
    code_excerpt = LINE_PREFIX.sub('', finding.get('code_excerpt') or '')
    finding['code_excerpt'] = code_excerpt
    if code_file not in ALLOWED_SOURCE_FILES:
        problems.append(f'code_file {code_file!r} is not one of the reviewed source files')
    elif not code_excerpt.strip():
        problems.append('code_excerpt is empty')
    elif code_excerpt not in (REPO_ROOT / code_file).read_text():
        problems.append(f'code_excerpt not found verbatim in {code_file}')

    for p in problems:
        print(f"WARNING: finding '{finding.get('id')}': {p}", file=sys.stderr)
    return problems


def apply_fix(finding):
    """Try to apply a literal doc fix in place. Returns None if applied, else the reason."""
    fix = finding.get('suggested_fix')
    if not fix:
        return 'no safe literal fix suggested'
    doc_path = REPO_ROOT / finding['doc_file']
    # Re-read each time: an earlier fix in this run may have changed the file.
    text = doc_path.read_text()
    old, new = fix['old_text'], fix['new_text']
    excerpt = finding.get('doc_excerpt') or ''
    count = text.count(old) if old else 0
    if finding.get('severity') != 'minor':
        # A serious disagreement may be a bug in the code, not the doc
        reason = 'not a minor finding: a human must decide whether the doc or the code is wrong'
    elif not old:
        reason = 'old_text is empty'
    elif old == new:
        reason = 'old_text and new_text are identical'
    elif old not in excerpt:
        reason = 'old_text does not lie inside doc_excerpt'
    elif excerpt not in text:
        reason = 'doc_excerpt no longer matches the doc (changed by an earlier fix)'
    elif count != 1:
        reason = f'old_text occurs {count} times in the doc, not exactly once'
    else:
        # old_text occurs once and lies inside doc_excerpt, which is in the doc, so the
        # single occurrence is necessarily the one inside the quoted passage.
        doc_path.write_text(text.replace(old, new))
        return None
    print(
        f"WARNING: finding '{finding['id']}': {reason}; leaving unfixed for manual review",
        file=sys.stderr,
    )
    return reason


def render_finding(finding, fix_status=None, problems=None, with_fingerprint=False):
    """Render one finding as markdown.

    fix_status is None when the fix was applied, else the reason it was not.
    problems lists the excerpt checks the finding failed, if any.
    with_fingerprint embeds the hidden dedupe marker (for issue posts).
    """
    heading = f'### `{finding["id"]}`'
    if finding.get('severity'):
        heading += f' — {finding["severity"]}'
    lines = [heading]
    if with_fingerprint:
        lines.append(FINGERPRINT_MARKER.format(finding['fingerprint']))
    if problems:
        lines += [
            '',
            '> [!WARNING]',
            '> The quotes in this finding could not be matched to the files, so it was not '
            'auto-fixed and its claims are unverified:',
        ]
        lines += [f'> - {p}' for p in problems]
        lines.append('')
    lines += [finding['description'], '']
    doc_excerpt = finding.get('doc_excerpt') or ''
    if doc_excerpt:
        lines.append(f'**Doc** (`{finding["doc_file"]}`):')
        lines += [f'> {line}' for line in doc_excerpt.splitlines()]
    else:
        lines.append(f'**Doc** (`{finding["doc_file"]}`): _no related passage_')
    lines += [
        '',
        f'**Code** (`{finding["code_file"]}`, {finding["code_lines"]}):',
        '```',
        finding['code_excerpt'],
        '```',
    ]
    fix = finding.get('suggested_fix')
    if fix:
        if fix_status is None:
            label = 'Applied fix'
        else:
            label = f'Suggested fix (not applied: {fix_status})'
        lines += [
            '',
            f'**{label}**',
            f'- old: `{fix["old_text"]}`',
            f'- new: `{fix["new_text"]}`',
        ]
    elif finding['type'] == 'inconsistency':
        lines += ['', '_No safe literal fix suggested — needs a manual edit._']
    lines.append('')
    return '\n'.join(lines)


def main():
    findings = load_findings()

    # Verify every excerpt against the unmodified files before any fix is applied.
    verified, unverified = [], []
    for f in findings:
        problems = verify_excerpts(f)
        f['fingerprint'] = fingerprint(f)
        if problems:
            unverified.append((f, problems))
        else:
            verified.append(f)

    inconsistencies = [f for f in verified if f['type'] == 'inconsistency']
    gaps = [f for f in verified if f['type'] == 'gap']
    inconsistencies.sort(key=lambda f: SEVERITY_ORDER.get(f.get('severity'), 99))
    unverified.sort(key=lambda fp: SEVERITY_ORDER.get(fp[0].get('severity'), 99))

    applied, serious, unapplied = [], [], []
    for f in inconsistencies:
        reason = apply_fix(f)
        if reason is None:
            applied.append(f)
        elif f.get('severity') == 'serious':
            serious.append((f, reason))
        else:
            unapplied.append((f, reason))

    pr_body_path = REPO_ROOT / 'pr_body.md'
    if applied:
        plural = 'y' if len(applied) == 1 else 'ies'
        pr_body_path.write_text(
            '## Documentation consistency check\n\n'
            f'Automated weekly check applied {len(applied)} fix(es) for minor '
            f'inconsistenc{plural} between `docs/Explanations/model.md` and the source code. '
            'Every claim below quotes the exact doc passage and code so it can be checked '
            'directly. Please verify each one before merging. Fixes are literal text '
            'replacements. Serious inconsistencies are never auto-fixed, since the code may '
            'be the side that is wrong. They are filed in the `docs-code-conflict` issue, '
            'and other inconsistencies without an applicable fix in the `docs-gaps` issue.\n\n'
            '## Checklist\n\n'
            '- [ ] I have verified each finding below against the actual doc and code\n'
            '- [ ] I have reviewed (and corrected if needed) every applied fix\n'
            '- [ ] I have marked this PR ready for review, which starts CI\n\n'
            '## Findings\n\n' + '\n'.join(render_finding(f) for f in applied)
        )
    else:
        pr_body_path.write_text(
            '## Documentation consistency check\n\nNo fixes applied this run.\n'
        )

    # Each issue dedupes against its own posts, so a finding first filed as minor is still
    # escalated if a later run rates it serious.
    (serious,), serious_already, serious_dups = drop_reported(
        [serious], 'existing_serious_issue.md'
    )
    # Verified findings come first, so they claim a fingerprint before an unverified copy.
    (unapplied, gaps, unverified), gap_already, gap_dups = drop_reported(
        [unapplied, [(f, None) for f in gaps], unverified], 'existing_issue.md'
    )
    gaps = [f for f, _ in gaps]

    serious_body = issue_body(
        [
            (
                '## Serious inconsistencies\n',
                'Automated weekly check found places where the doc and the code disagree in a '
                'way that would change model output or physical interpretation. For each, '
                'decide which side is correct, then fix the code or the doc.\n',
                [
                    render_finding(f, fix_status=reason, with_fingerprint=True)
                    for f, reason in serious
                ],
            ),
        ],
        serious_already,
    )
    (REPO_ROOT / 'serious_issue_body.md').write_text(serious_body)

    gap_body = issue_body(
        [
            (
                '## Inconsistencies needing a manual edit\n',
                'Automated weekly check found doc passages that contradict the source code '
                'but could not be fixed automatically. Each finding states why its fix was '
                'not applied.\n',
                [
                    render_finding(f, fix_status=reason, with_fingerprint=True)
                    for f, reason in unapplied
                ],
            ),
            (
                '## Documentation gaps\n',
                'Automated weekly check found source-code behaviour with no corresponding '
                'documentation. These are not auto-fixed — writing new model-description '
                'prose needs a human who can vouch for the physics.\n',
                [render_finding(f, fix_status='gap', with_fingerprint=True) for f in gaps],
            ),
            (
                '## Unverified findings\n',
                'These findings quote doc or code text that does not occur verbatim in the '
                'named file. These findings are suspect: check manually. None of their fixes '
                'were applied.\n',
                [
                    render_finding(
                        f,
                        fix_status='finding unverified',
                        problems=problems,
                        with_fingerprint=True,
                    )
                    for f, problems in unverified
                ],
            ),
        ],
        gap_already,
    )
    (REPO_ROOT / 'gap_issue_body.md').write_text(gap_body)

    summary = (
        f'Docs consistency triage: {len(applied)} fix(es) applied, '
        f'{len(serious)} serious inconsistenc(ies), '
        f'{len(unapplied)} inconsistenc(ies) need a manual edit, {len(gaps)} gap(s), '
        f'{len(unverified)} unverified finding(s)'
    )
    if unverified:
        summary += ' (' + ', '.join(f'`{f.get("id")}`' for f, _ in unverified) + ')'
    summary += (
        f' newly filed; {serious_already + gap_already} already reported, '
        f'{serious_dups + gap_dups} duplicate(s) within this run'
    )
    print(summary)
    step_summary = os.environ.get('GITHUB_STEP_SUMMARY')
    if step_summary:
        with open(step_summary, 'a') as f:
            f.write(summary + '\n')

    results = {
        'has_fixes': 'true' if applied else 'false',
        'has_serious_items': 'true' if serious_body else 'false',
        'has_issue_items': 'true' if gap_body else 'false',
    }
    github_output = os.environ.get('GITHUB_OUTPUT')
    if github_output:
        with open(github_output, 'a') as f:
            for key, value in results.items():
                f.write(f'{key}={value}\n')
    else:
        for key, value in results.items():
            print(f'{key}={value}')


if __name__ == '__main__':
    main()
