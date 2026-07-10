#!/usr/bin/env bash
# Test-structure gate for the JANUS test suite (nested "Shape A" layout:
# sources under src/janus/, tests under tests/<subdir>/test_<file>.py).
#
# Two checks:
#
#   (a) Module-level marker (HARD FAIL). Every tests/**/test_*.py file must
#       declare a module-level `pytestmark` that contains exactly one tier
#       marker: unit / smoke / integration / slow. A file with no tier
#       marker is invisible to CI's marker-filtered runs.
#
#   (b) Source mirror (WARN ONLY). A test in a subdirectory,
#       tests/<subdir>/test_<stem>.py, is expected to mirror the source
#       src/janus/<subdir>/<stem>.py. JANUS has many sources and only a
#       subset are tested, so a missing mirror is a warning, not a failure.
#       Root-level cross-cutting harness tests (integration / slow tests that
#       exercise several modules at once) do not mirror a single source and
#       are tolerated.
#
# Run from repository root:
#
#     bash tools/validate_test_structure.sh
#
# Exit 1 only when a module-level tier marker is missing.

set -e

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

if [ ! -d tests ]; then
    echo "ERROR: tests/ not found in $REPO_ROOT" >&2
    exit 2
fi

python3 - <<'PY'
import ast
import pathlib
import sys

REPO = pathlib.Path('.').resolve()
TESTS = REPO / 'tests'
SRC = REPO / 'src' / 'janus'
TIER = {'unit', 'smoke', 'integration', 'slow'}

# Root-level cross-cutting harness tests: integration / slow tests that
# exercise multi-module behaviour and intentionally do not mirror a single
# source file. The mirror check tolerates these without warning.
KNOWN_HARNESS = {
    'test_instellation.py',
    'test_runaway_greenhouse.py',
    'test_code.py',
    'test_constants.py',
}


def mark_names_from_nodes(nodes):
    out = []
    for d in nodes:
        node = d.func if isinstance(d, ast.Call) else d
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Attribute):
            if (
                isinstance(node.value.value, ast.Name)
                and node.value.value.id == 'pytest'
                and node.value.attr == 'mark'
            ):
                out.append(node.attr)
    return out


def module_tier_markers(tree):
    tiers = []
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id == 'pytestmark':
                val = node.value
                items = val.elts if isinstance(val, (ast.List, ast.Tuple)) else [val]
                for name in mark_names_from_nodes(items):
                    if name in TIER:
                        tiers.append(name)
    return tiers


hard_failures = []
warnings = []
total = 0

for path in sorted(TESTS.rglob('test_*.py')):
    if any(part.startswith('.') or part == '__pycache__' for part in path.parts):
        continue
    total += 1
    rel = path.relative_to(REPO)
    try:
        tree = ast.parse(path.read_text(encoding='utf-8'))
    except SyntaxError as exc:
        hard_failures.append(f'{rel}:{exc.lineno}: SyntaxError {exc.msg}')
        continue

    # (a) Module-level tier marker (hard fail).
    tiers = module_tier_markers(tree)
    unique_tiers = sorted(set(tiers))
    if not unique_tiers:
        hard_failures.append(
            f'{rel}: no module-level pytestmark tier marker '
            '(need one of unit / smoke / integration / slow)'
        )
    elif len(unique_tiers) > 1:
        hard_failures.append(
            f'{rel}: module-level pytestmark declares multiple tier markers '
            f'{unique_tiers}; exactly one is required.'
        )

    # (b) Source mirror (warn only).
    parts = path.relative_to(TESTS).parts
    if len(parts) >= 2:
        subdir = pathlib.Path(*parts[:-1])
        stem = path.name[len('test_'):]
        source = SRC / subdir / stem
        if not source.exists():
            warnings.append(
                f'{rel}: no mirrored source src/janus/{subdir}/{stem}'
            )
    elif path.name not in KNOWN_HARNESS:
        warnings.append(
            f'{rel}: root-level test with no mirrored source '
            '(cross-cutting harness test?)'
        )

if warnings:
    print(f'[!] Mirror warnings ({len(warnings)}):')
    for w in warnings:
        print(f'    {w}')
    print()

if hard_failures:
    print(f'Marker validation FAILED on {len(hard_failures)} file(s):')
    for f in hard_failures:
        print(f'  {f}')
    sys.exit(1)

print(f'Marker validation OK: {total} test file(s), all carry a module-level tier marker.')
PY
