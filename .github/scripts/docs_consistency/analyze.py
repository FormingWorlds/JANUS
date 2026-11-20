#!/usr/bin/env python3
"""Compare JANUS's model-description docs against the source code via a headless
Claude Code CLI call authenticated with CLAUDE_CODE_OAUTH_TOKEN, then let a second,
smaller call reword the suggested fixes so they read clearly, and write the findings
to analysis.json at the repo root.
"""

import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent

DOC_FILES = [
    'docs/Explanations/model.md',
]

SOURCE_FILES = [
    'src/janus/modules/compute_moist_adiabat.py',
    'src/janus/modules/dry_adiabat_setup.py',
    'src/janus/modules/dry_adiabat_timestep.py',
    'src/janus/modules/find_tropopause.py',
    'src/janus/modules/relative_humidity.py',
    'src/janus/modules/set_stratosphere.py',
    'src/janus/modules/solve_pt.py',
    'src/janus/modules/spectral_planck_surface.py',
    'src/janus/utils/atmosphere_column.py',
    'src/janus/utils/cp_funcs.py',
    'src/janus/utils/GeneralAdiabat.py',
    'src/janus/utils/height.py',
    'src/janus/utils/phys.py',
    'src/janus/utils/RayleighSpectrum.py',
    'src/janus/utils/water_tables.py',
]

MODEL = os.environ.get('CLAUDE_MODEL', 'claude-opus-5-5')
EFFORT = os.environ.get('CLAUDE_EFFORT', 'medium')

# The doc and source text are inlined into the prompt. The call gets no tools,
# no settings files and no MCP servers, and runs from an empty directory so no
# CLAUDE.md is seen.
ISOLATION_FLAGS = ['--tools', '', '--setting-sources', '', '--strict-mcp-config']

MAX_TURNS = 3
TIMEOUT_S = int(os.environ.get('CLAUDE_TIMEOUT_S', '1200'))
# The readability call sees only the fixes; both calls together stay below the job's
# 30 min limit.
READABILITY_TIMEOUT_S = int(os.environ.get('CLAUDE_READABILITY_TIMEOUT_S', '300'))

FINDINGS_SCHEMA = {
    'type': 'object',
    'properties': {
        'findings': {
            'type': 'array',
            'items': {
                'type': 'object',
                'properties': {
                    'id': {
                        'type': 'string',
                        'description': "short unique slug, e.g. 'lapse-rate-sign'",
                    },
                    'type': {'type': 'string', 'enum': ['inconsistency', 'gap']},
                    'severity': {
                        'type': 'string',
                        'enum': ['serious', 'minor'],
                        'description': 'required for type=inconsistency, omit for type=gap',
                    },
                    'doc_file': {'type': 'string'},
                    'doc_excerpt': {
                        'type': 'string',
                        'description': 'verbatim quote from doc_file',
                    },
                    'code_file': {'type': 'string'},
                    'code_lines': {'type': 'string', 'description': "e.g. 'L42-L58'"},
                    'code_excerpt': {
                        'type': 'string',
                        'description': 'verbatim quote from code_file',
                    },
                    'description': {
                        'type': 'string',
                        'description': "why these disagree, or what's missing",
                    },
                    'suggested_fix': {
                        'type': 'object',
                        'description': 'omit entirely if no safe literal fix can be given',
                        'properties': {
                            'old_text': {
                                'type': 'string',
                                'description': (
                                    'verbatim substring of doc_excerpt that occurs exactly '
                                    'once in doc_file'
                                ),
                            },
                            'new_text': {'type': 'string', 'description': 'replacement text'},
                        },
                        'required': ['old_text', 'new_text'],
                    },
                },
                'required': [
                    'id',
                    'type',
                    'doc_file',
                    'doc_excerpt',
                    'code_file',
                    'code_lines',
                    'code_excerpt',
                    'description',
                ],
            },
        },
    },
    'required': ['findings'],
}

READABILITY_SCHEMA = {
    'type': 'object',
    'properties': {
        'fixes': {
            'type': 'array',
            'items': {
                'type': 'object',
                'properties': {
                    'id': {'type': 'string', 'description': 'id of the finding'},
                    'new_text': {
                        'type': 'string',
                        'description': 'the reworded new_text, or the given one unchanged',
                    },
                },
                'required': ['id', 'new_text'],
            },
        },
    },
    'required': ['fixes'],
}


# What a reworded fix must keep from the original: its numbers, inline and display
# math, code spans and footnote references, in any order.
KEPT_PATTERNS = [
    re.compile(r'\d+(?:\.\d+)?'),
    re.compile(r'\$\$.+?\$\$|\$.+?\$', flags=re.DOTALL),
    re.compile(r'(`+).+?\1', flags=re.DOTALL),
    re.compile(r'\[\^[^\]]+\]'),
]


def read_numbered(rel_path):
    text = (REPO_ROOT / rel_path).read_text()
    numbered = '\n'.join(f'{i + 1}: {line}' for i, line in enumerate(text.splitlines()))
    return f'--- {rel_path} ---\n{numbered}'


def read_plain(rel_path):
    return f'--- {rel_path} ---\n{(REPO_ROOT / rel_path).read_text()}'


def build_prompt():
    template = (SCRIPT_DIR / 'prompt_template.md').read_text()
    # Docs go in unnumbered: any prefix the model copied from a numbered listing
    # would make an true quote fail to match.
    docs_blob = '\n\n'.join(read_plain(p) for p in DOC_FILES)
    source_blob = '\n\n'.join(read_numbered(p) for p in SOURCE_FILES)
    return template.replace('{{DOCS}}', docs_blob).replace('{{SOURCE}}', source_blob)


def run_claude(prompt, schema, timeout_s):
    """Run one isolated CLI call and return its structured output.

    Raises RuntimeError with the reason when the call fails or returns no output
    matching the schema.
    """
    # The prompt (docs + 15 source files inlined) can be well over argv's
    # OS-level size limit, so it goes in over stdin rather than as an argument.
    try:
        with tempfile.TemporaryDirectory() as empty_dir:
            proc = subprocess.run(
                [
                    'claude',
                    '--print',
                    '--output-format',
                    'json',
                    '--json-schema',
                    json.dumps(schema),
                    '--model',
                    MODEL,
                    '--effort',
                    EFFORT,
                    '--max-turns',
                    str(MAX_TURNS),
                    '--permission-mode',
                    'dontAsk',
                    *ISOLATION_FLAGS,
                ],
                input=prompt,
                capture_output=True,
                text=True,
                cwd=empty_dir,
                timeout=timeout_s,
            )
    except subprocess.TimeoutExpired:
        raise RuntimeError(f'claude CLI did not finish within {timeout_s} s') from None

    if proc.returncode != 0:
        raise RuntimeError(
            f'claude CLI exited {proc.returncode}\n'
            f'stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}'
        )

    try:
        envelope = json.loads(proc.stdout)
    except json.JSONDecodeError:
        raise RuntimeError(
            f'claude CLI did not print valid JSON on stdout:\n{proc.stdout}'
        ) from None

    if envelope.get('is_error') or envelope.get('subtype') != 'success':
        raise RuntimeError(f'claude CLI run did not succeed: {envelope}')

    structured = envelope.get('structured_output')
    key = next(iter(schema['required']))
    if not isinstance(structured, dict) or key not in structured:
        raise RuntimeError(
            f'claude CLI response had no structured_output.{key} matching the schema; '
            f'raw result text was: {envelope.get("result")!r}'
        )
    return structured


def build_readability_prompt(fixed):
    template = (SCRIPT_DIR / 'readability_template.md').read_text()
    items = [
        {
            'id': f['id'],
            'doc_excerpt': f['doc_excerpt'],
            'code_excerpt': f['code_excerpt'],
            'description': f['description'],
            'old_text': f['suggested_fix']['old_text'],
            'new_text': f['suggested_fix']['new_text'],
        }
        for f in fixed
    ]
    return template.replace('{{FIXES}}', json.dumps(items, indent=2))


# Text structure a reworded fix must keep: blank lines between paragraphs, and the
# markers that start a list item, heading, table row, admonition or quote.
BLANK_LINE_RE = re.compile(r'\n[ \t]*\n')
LINE_MARKER_RE = re.compile(r'^[ \t]*([-*+]|\d+\.|#{1,6}|\||!!!|>)[ \t]', flags=re.MULTILINE)


def continuation_indents(text):
    """The indentation of every non-blank line after the first."""
    return {
        line[: len(line) - len(line.lstrip())] for line in text.split('\n')[1:] if line.strip()
    }


def keeps_content(original, rewrite, indent=''):
    """Whether rewrite keeps the numbers, math, code spans, footnote references and
    structure of original.

    Rewrapping may add or remove lines, but every new line must use an indentation the
    original uses (indent, that of the doc line the fix starts on, when the original is
    one line), so text inside a note stays inside it. Paragraph breaks and line-start
    markers must stay the same.
    """
    if not continuation_indents(rewrite) <= (continuation_indents(original) or {indent}):
        return False
    if len(BLANK_LINE_RE.findall(original)) != len(BLANK_LINE_RE.findall(rewrite)):
        return False
    if sorted(LINE_MARKER_RE.findall(original)) != sorted(LINE_MARKER_RE.findall(rewrite)):
        return False
    return all(
        sorted(m[0] for m in pattern.finditer(original))
        == sorted(m[0] for m in pattern.finditer(rewrite))
        for pattern in KEPT_PATTERNS
    )


def apply_rewrites(findings, rewrites):
    """Put the reworded new_text back on the findings it belongs to.

    A rewrite is used only when its id names exactly one finding with a fix, its text
    changes more than white space, does not undo the fix, and keeps the content of the
    original new_text (see keeps_content). Returns the number of fixes that changed.
    """
    by_id = {}
    for f in findings:
        if f.get('suggested_fix'):
            by_id.setdefault(f['id'], []).append(f)
    changed = 0
    for r in rewrites:
        matches = by_id.get(r['id'], [])
        new = r['new_text']
        if len(matches) != 1:
            continue
        fix = matches[0]['suggested_fix']
        excerpt = matches[0].get('doc_excerpt', '')
        start = max(excerpt.find(fix['old_text']), 0)
        line = excerpt[excerpt.rfind('\n', 0, start) + 1 : start]
        indent = line[: len(line) - len(line.lstrip())]
        words = new.split()
        if not words or words in (fix['old_text'].split(), fix['new_text'].split()):
            continue
        if not keeps_content(fix['new_text'], new, indent):
            print(f'Kept the original wording of {r["id"]}: the rewording changed its content')
            continue
        fix['new_text'] = new
        changed += 1
    return changed


def main():
    if not os.environ.get('CLAUDE_CODE_OAUTH_TOKEN'):
        print('CLAUDE_CODE_OAUTH_TOKEN is not set', file=sys.stderr)
        sys.exit(1)

    try:
        findings = run_claude(build_prompt(), FINDINGS_SCHEMA, TIMEOUT_S)['findings']
    except RuntimeError as err:
        print(err, file=sys.stderr)
        sys.exit(1)

    # The readability check only rewords; if it fails, the fixes keep their wording.
    fixed = [f for f in findings if f.get('suggested_fix')]
    if fixed:
        try:
            rewrites = run_claude(
                build_readability_prompt(fixed), READABILITY_SCHEMA, READABILITY_TIMEOUT_S
            )['fixes']
            changed = apply_rewrites(findings, rewrites)
            print(f'Readability check reworded {changed} of {len(fixed)} fixes')
        except Exception as err:  # any failure here keeps the findings as they are
            print(f'Readability check skipped: {err!r}', file=sys.stderr)

    out_path = REPO_ROOT / 'analysis.json'
    out_path.write_text(json.dumps(findings, indent=2))
    print(f'Wrote {len(findings)} findings to {out_path}')


if __name__ == '__main__':
    main()
