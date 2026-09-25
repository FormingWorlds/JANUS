#!/usr/bin/env python3
"""Compare JANUS's model-description docs against the source code via one
headless Claude Code CLI call authenticated with CLAUDE_CODE_OAUTH_TOKEN,
and write the resulting findings to analysis.json at the repo root.
"""

import json
import os
import subprocess
import sys
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

# Tools the CLI could otherwise use; the doc and source text are already
# inlined into the prompt below, so none of these are needed and disallowing
# them keeps this call to a single read-only, side-effect-free turn.
DISALLOWED_TOOLS = 'Bash,Edit,Write,NotebookEdit,WebFetch,WebSearch,Read,Glob,Grep'

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


def main():
    if not os.environ.get('CLAUDE_CODE_OAUTH_TOKEN'):
        print('CLAUDE_CODE_OAUTH_TOKEN is not set', file=sys.stderr)
        sys.exit(1)

    prompt = build_prompt()

    # The prompt (docs + 15 source files inlined) can be well over argv's
    # OS-level size limit, so it goes in over stdin rather than as an argument.
    proc = subprocess.run(
        [
            'claude',
            '--print',
            '--output-format',
            'json',
            '--json-schema',
            json.dumps(FINDINGS_SCHEMA),
            '--model',
            MODEL,
            '--effort',
            EFFORT,
            '--permission-mode',
            'dontAsk',
            '--disallowedTools',
            DISALLOWED_TOOLS,
        ],
        input=prompt,
        capture_output=True,
        text=True,
    )

    if proc.returncode != 0:
        print(
            f'claude CLI exited {proc.returncode}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}',
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        envelope = json.loads(proc.stdout)
    except json.JSONDecodeError:
        print(f'claude CLI did not print valid JSON on stdout:\n{proc.stdout}', file=sys.stderr)
        sys.exit(1)

    if envelope.get('is_error') or envelope.get('subtype') != 'success':
        print(f'claude CLI run did not succeed: {envelope}', file=sys.stderr)
        sys.exit(1)

    structured = envelope.get('structured_output')
    if not isinstance(structured, dict) or 'findings' not in structured:
        print(
            'claude CLI response had no structured_output.findings matching the schema; '
            f'raw result text was: {envelope.get("result")!r}',
            file=sys.stderr,
        )
        sys.exit(1)

    findings = structured['findings']

    out_path = REPO_ROOT / 'analysis.json'
    out_path.write_text(json.dumps(findings, indent=2))
    print(f'Wrote {len(findings)} findings to {out_path}')


if __name__ == '__main__':
    main()
