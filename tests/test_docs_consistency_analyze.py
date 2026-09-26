"""Tests for .github/scripts/docs_consistency/analyze.py.

Mocks the Claude CLI subprocess and checks the contract an unattended weekly run relies
on. See docs/How-to/test.md.
"""

import importlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

SCRIPT_DIR = Path(__file__).resolve().parents[1] / '.github' / 'scripts' / 'docs_consistency'
sys.path.insert(0, str(SCRIPT_DIR))
analyze = importlib.import_module('analyze')

FINDINGS = [{'id': 'x', 'type': 'gap'}]
GOOD = {'subtype': 'success', 'is_error': False, 'structured_output': {'findings': FINDINGS}}


@pytest.fixture
def run_cli(tmp_path, monkeypatch, capsys):
    """Run analyze.main() against a fake CLI result; return (exit code, stderr, call)."""
    # REPO_ROOT is a module-level constant computed at import, so patch it directly.
    monkeypatch.setattr(analyze, 'REPO_ROOT', tmp_path)
    monkeypatch.setattr(analyze, 'build_prompt', lambda: 'the prompt')
    monkeypatch.setenv('CLAUDE_CODE_OAUTH_TOKEN', 'token')

    def run(result):
        call = {}

        def fake_run(cmd, **kwargs):
            call.update(cmd=cmd, cwd_listing=os.listdir(kwargs['cwd']), **kwargs)
            if isinstance(result, Exception):
                raise result
            return result

        monkeypatch.setattr(analyze.subprocess, 'run', fake_run)
        try:
            analyze.main()
            code = 0
        except SystemExit as e:
            code = e.code
        return code, capsys.readouterr().err, call

    return run


def completed(stdout, returncode=0, stderr=''):
    """A finished CLI process with the given output."""
    return subprocess.CompletedProcess(['claude'], returncode, stdout, stderr)


@pytest.mark.parametrize(
    ('result', 'message'),
    [
        (subprocess.TimeoutExpired('claude', 1200), 'did not finish within 1200 s'),
        (completed('', returncode=1, stderr='Invalid API key'), 'claude CLI exited 1'),
        (completed('not json at all'), 'did not print valid JSON'),
        (completed(json.dumps({**GOOD, 'is_error': True})), 'run did not succeed'),
        (completed(json.dumps({**GOOD, 'subtype': 'error_max_turns'})), 'run did not succeed'),
        (completed(json.dumps({'subtype': 'success', 'result': 'hi'})), 'no structured_output'),
        (
            completed(json.dumps({**GOOD, 'structured_output': {'items': []}})),
            'no structured_output',
        ),
    ],
    ids=[
        'timeout',
        'non-zero-exit',
        'invalid-json',
        'error-envelope',
        'max-turns-reached',
        'structured-output-missing',
        'findings-key-missing',
    ],
)
def test_every_failure_exits_nonzero_and_writes_nothing(run_cli, tmp_path, result, message):
    """A failed or unusable CLI run fails the job loudly and leaves no analysis.json."""
    code, err, _ = run_cli(result)
    assert code == 1
    assert message in err
    assert not (tmp_path / 'analysis.json').exists()


def test_missing_token_fails_before_calling_the_cli(run_cli, monkeypatch, tmp_path):
    """Without CLAUDE_CODE_OAUTH_TOKEN the run stops before starting the CLI."""
    monkeypatch.delenv('CLAUDE_CODE_OAUTH_TOKEN')
    code, err, call = run_cli(completed(json.dumps(GOOD)))
    assert code == 1 and 'CLAUDE_CODE_OAUTH_TOKEN is not set' in err
    assert call == {} and not (tmp_path / 'analysis.json').exists()


def test_good_response_writes_findings_from_an_isolated_cli_call(run_cli, tmp_path):
    """A successful run writes the findings, and the CLI was started isolated with a timeout."""
    code, _, call = run_cli(completed(json.dumps(GOOD)))
    assert code == 0
    assert json.loads((tmp_path / 'analysis.json').read_text()) == FINDINGS
    cmd = call['cmd']
    assert cmd[cmd.index('--tools') + 1] == '' and '--strict-mcp-config' in cmd
    assert cmd[cmd.index('--setting-sources') + 1] == ''
    assert call['cwd_listing'] == [] and call['cwd'] != str(tmp_path)
    assert call['timeout'] == analyze.TIMEOUT_S and call['input'] == 'the prompt'
