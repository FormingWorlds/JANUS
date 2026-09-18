#!/usr/bin/env python3
"""Deterministic triage of analysis.json: apply literal doc fixes, and write
pr_body.md / gap_issue_body.md for the workflow to use. No LLM calls: everything is a plain string operation.
"""
import json
import os
import sys
from pathlib import Path

from analyze import DOC_FILES

REPO_ROOT = Path(__file__).resolve().parents[3]

SEVERITY_ORDER = {"serious": 0, "minor": 1}
ALLOWED_DOC_FILES = set(DOC_FILES)


def load_findings():
    return json.loads((REPO_ROOT / "analysis.json").read_text())


def apply_fix(finding):
    """Try to apply a literal doc fix in place. Returns True if applied."""
    fix = finding.get("suggested_fix")
    if not fix:
        return False
    doc_file = finding.get("doc_file")
    if doc_file not in ALLOWED_DOC_FILES:
        print(
            f"WARNING: finding '{finding['id']}' targets doc_file={doc_file!r}, which is not "
            f"one of the reviewed doc files {sorted(ALLOWED_DOC_FILES)}; leaving unfixed for manual review",
            file=sys.stderr,
        )
        return False
    doc_path = REPO_ROOT / doc_file
    text = doc_path.read_text()
    old, new = fix["old_text"], fix["new_text"]
    if old not in text:
        print(
            f"WARNING: old_text for finding '{finding['id']}' not found verbatim "
            f"in {finding['doc_file']}; leaving unfixed for manual review",
            file=sys.stderr,
        )
        return False
    doc_path.write_text(text.replace(old, new, 1))
    return True


def render_finding(finding, applied):
    heading = f"### `{finding['id']}`"
    if finding.get("severity"):
        heading += f" — {finding['severity']}"
    lines = [
        heading,
        finding["description"],
        "",
        f"**Doc** (`{finding['doc_file']}`):",
        f"> {finding['doc_excerpt']}",
        "",
        f"**Code** (`{finding['code_file']}`, {finding['code_lines']}):",
        "```",
        finding["code_excerpt"],
        "```",
    ]
    fix = finding.get("suggested_fix")
    if fix:
        label = "Applied fix" if applied else "Suggested fix (verbatim match not found — needs manual edit)"
        lines += [
            "",
            f"**{label}**",
            f"- old: `{fix['old_text']}`",
            f"- new: `{fix['new_text']}`",
        ]
    elif finding["type"] == "inconsistency":
        lines += ["", "_No safe literal fix suggested — needs a manual edit._"]
    lines.append("")
    return "\n".join(lines)


def main():
    findings = load_findings()
    inconsistencies = [f for f in findings if f["type"] == "inconsistency"]
    gaps = [f for f in findings if f["type"] == "gap"]
    inconsistencies.sort(key=lambda f: SEVERITY_ORDER.get(f.get("severity"), 99))

    applied_count = 0
    rendered = []
    for f in inconsistencies:
        applied = apply_fix(f)
        applied_count += int(applied)
        rendered.append(render_finding(f, applied))

    pr_body_path = REPO_ROOT / "pr_body.md"
    if inconsistencies:
        plural = "y" if len(inconsistencies) == 1 else "ies"
        pr_body_path.write_text(
            "## Documentation consistency check\n\n"
            f"Automated weekly check found {len(inconsistencies)} inconsistenc{plural} "
            "between `docs/Explanations/model.md` and the source code, ordered serious → minor. "
            "Every claim below quotes the exact doc passage and code so it can be checked directly — "
            f"please verify each one before merging. {applied_count} fix(es) were applied automatically "
            "as literal text replacements; the rest need a manual edit.\n\n"
            "## Checklist\n\n"
            "- [ ] I have verified each finding below against the actual doc and code\n"
            "- [ ] I have reviewed (and corrected if needed) every applied fix\n"
            "- [ ] I have addressed or dismissed every fix-less finding\n\n"
            "## Findings\n\n" + "\n".join(rendered)
        )
    else:
        pr_body_path.write_text("## Documentation consistency check\n\nNo inconsistencies found this run.\n")

    gap_body_path = REPO_ROOT / "gap_issue_body.md"
    if gaps:
        gap_lines = [
            "Automated weekly check found source-code behaviour with no corresponding documentation. "
            "These are not auto-fixed — writing new model-description prose needs a human who can vouch "
            "for the physics.\n",
        ]
        gap_lines += [render_finding(f, applied=False) for f in gaps]
        gap_body_path.write_text("\n".join(gap_lines))
    else:
        gap_body_path.write_text("")

    results = {
        "has_fixes": "true" if inconsistencies else "false",
        "has_gaps": "true" if gaps else "false",
    }
    github_output = os.environ.get("GITHUB_OUTPUT")
    if github_output:
        with open(github_output, "a") as f:
            for key, value in results.items():
                f.write(f"{key}={value}\n")
    else:
        for key, value in results.items():
            print(f"{key}={value}")


if __name__ == "__main__":
    main()
