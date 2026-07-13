# JANUS AI Agent Guidelines

**Trust these instructions.** Only search if information is incomplete or found to be in error.

**Identity & Mission**: You are an expert Scientific Software Engineer working on the JANUS module of the PROTEUS ecosystem.

## High-Level Instructions

> ### Rule files you MUST read on every session
>
> JANUS keeps its Claude-Code rule files under `.github/.claude/rules/` (NOT the conventional repo-root `.claude/`, which is gitignored and so cannot be shared with collaborators). Claude Code does NOT auto-discover the rules at this unusual path. Read them explicitly at the start of every session and any time you open a related file:
>
> - [`.github/.claude/rules/janus-tests.md`](.claude/rules/janus-tests.md) -- test quality deep-dive: anti-happy-path patterns, discriminating-value guards, physics-invariant tiering, validation certification markers, adversarial-review trigger, `cp_mode` default flip, `satvps` argument-swap, hypothesis seed stability, pipeline solver intermediate-state assertions. **Required reading before editing any file under `tests/**` or `src/janus/**`.**
> - [`.github/.claude/rules/janus-code-review.md`](.claude/rules/janus-code-review.md) -- review-pass gate, domain-aware physics review (spectral unit boundaries, per-band normalization, dry-adiabat exponent safety, PROTEUS-coupling contract). **Required reading before any code review pass.**
>
> These two files plus this one are the canonical sources of truth for testing rigor and review criteria. Together they enforce JANUS's extreme-rigor stance on physics validity, anti-happy-path testing, and validation certification.

1. **Always** read the two rule files above plus the Testing Standards section below before any code change.
2. **Always** inform the user that you are reading in this file by printing a message at the start of your response: "(Read in copilot-instructions.md...)"
3. When creating a PR, **always** follow the PR template (`.github/pull_request_template.md`) and ensure all sections are filled out with relevant information.
4. **Claude-specific**: `CLAUDE.md` is a symlink to this file. Session learnings, plans, and memories live in `~/.claude/projects/<repo>/memory/`; they do NOT live in this repository.

## Ecosystem Context

JANUS is the 1D convective / radiative atmosphere module of the PROTEUS ecosystem. It generates the temperature-pressure structure of a planetary atmosphere and its outgoing radiation, and is called by the main [PROTEUS](https://github.com/FormingWorlds/PROTEUS) coupled atmosphere-interior framework through the `atmos_clim` step. JANUS is also usable standalone for atmosphere-structure studies.

PROTEUS is a coupled atmosphere-interior framework with a modular architecture:

- **[PROTEUS](https://github.com/FormingWorlds/PROTEUS)** (main repository): Core coupling framework and orchestration
- **[AGNI](https://github.com/nichollsh/AGNI)**: Radiative-convective atmospheric energy module (Julia)
- **[SOCRATES](https://github.com/FormingWorlds/SOCRATES)**: Spectral radiative transfer code (Fortran)
- **[CALLIOPE](https://github.com/FormingWorlds/CALLIOPE)**: Volatile in-/outgassing and thermodynamics module (Python)
- **[JANUS](https://github.com/FormingWorlds/JANUS)**: 1D convective atmosphere module (Python)
- **[MORS](https://github.com/FormingWorlds/MORS)**: Stellar evolution module (Python)
- **[ARAGOG](https://github.com/FormingWorlds/aragog)**: Interior thermal evolution module based on T-P formalism (Python)
- **[SPIDER](https://github.com/FormingWorlds/SPIDER)**: Interior thermal evolution module based on T-S formalism (C)
- **[VULCAN](https://github.com/FormingWorlds/VULCAN)**: Atmospheric chemistry module (Python)
- **[ZEPHYRUS](https://github.com/FormingWorlds/ZEPHYRUS)**: Atmospheric escape module (Python)
- **[BOREAS](https://github.com/FormingWorlds/BOREAS)**: Hydrodynamic atmospheric escape module (Python)
- **[Obliqua](https://github.com/FormingWorlds/Obliqua)**: Tidal evolution module (Julia)

**Project Type**: Scientific simulation module (Python).

**Languages**: Python, with a compiled Fortran radiative-transfer dependency (SOCRATES).

**Distribution**: PyPI package `fwl-janus`; import package `janus`.

**Target Runtime**: Python 3.11 to 3.13 on Linux / macOS (3.12 primary; CI runs the matrix 3.11, 3.12, 3.13).

## Build & Validation

### Environment Setup

**Prerequisites**:

1. Python 3.11 to 3.13 (via conda / miniforge or system).
2. A Fortran / C build toolchain (`gfortran`, `gcc`, `make`) and NetCDF-Fortran development headers, for SOCRATES.
3. Git.

**Developer Install**:

```bash
# 1. Install and build SOCRATES (Fortran radiative transfer), then point RAD_DIR at it
git clone git@github.com:FormingWorlds/JANUS.git
cd JANUS
bash tools/get_socrates.sh /path/to/socrates
export RAD_DIR=/path/to/socrates

# 2. Install JANUS as editable with the develop extras (pulls in fwl-mors, pytest, coverage)
pip install -e ".[develop]"
pre-commit install -f

# 3. Point FWL_DATA at a data directory and download spectral + stellar data
export FWL_DATA=/path/to/fwl_data
janus download spectral
janus download stellar
```

JANUS has one compiled dependency, SOCRATES, built via `tools/get_socrates.sh`. Both `RAD_DIR` (compiled SOCRATES tree) and `FWL_DATA` (spectral / stellar data) must be set before running or testing. When working alongside PROTEUS, the recommended pattern is editable-install via `pip install -e JANUS/.` from the PROTEUS tree so changes propagate without re-install.

### Test Commands

**Run all tests**:

```bash
pytest
```

**Run by category** (the marker tiers this framework uses):

```bash
pytest -m "(unit or smoke) and not skip"     # What PR checks run
pytest -m unit                                # Fast unit tests (<100ms each)
pytest -m smoke                               # Real SOCRATES binary, low resolution, 1 step
pytest -m integration                         # Full pipeline coupling (nightly)
pytest -m slow                                # Full physics validation (nightly)
```

**With coverage** (matches CI):

```bash
coverage run -m pytest
coverage report
coverage html
```

**Coverage thresholds** (in `pyproject.toml`; auto-ratcheting, never manually decreased, capped at 90 by `tools/update_coverage_threshold.py`):

- Fast gate (`[tool.janus.coverage_fast]`, unit + smoke, every PR): ratcheting toward **90%** (the PROTEUS-ecosystem ceiling).
- Full gate (`[tool.coverage.report]`, unit + smoke + integration + slow, nightly): **90%**.

**Validate test structure**:

```bash
bash tools/validate_test_structure.sh
```

**Test quality lint** (blocking on PRs):

```bash
python tools/check_test_quality.py --check
```

### Lint Commands

**Always run before committing**:

```bash
ruff check src/ tests/        # Check for issues
ruff check --fix src/ tests/  # Auto-fix issues
ruff format src/ tests/       # Format code
```

**Pre-commit hook** (runs automatically on commit):

```bash
pre-commit install -f
```

### Validation Pipeline

**CI runs on PRs** (`.github/workflows/tests.yaml`):

1. **Build SOCRATES**: checkout and compile the Fortran radiative-transfer dependency (cached across runs).
2. **Build JANUS**: `pip install -e .[develop]` across the Python 3.11 / 3.12 / 3.13 matrix.
3. **Tests with coverage**: `coverage run -m pytest` with `RAD_DIR` and `FWL_DATA` exported.
4. **Coverage report**: totals surfaced in the job summary.

As the test framework is adopted, the PR checks add: the fast coverage gate (`[tool.janus.coverage_fast].fail_under` against unit + smoke coverage), the structure check (`bash tools/validate_test_structure.sh`), the blocking test-quality lint (`python tools/check_test_quality.py --check`), and the coverage ratchet guard that rejects any PR lowering `[tool.coverage.report].fail_under` below `min(base_ref, 90.0)`.

**All must pass** before merge. Coverage thresholds auto-ratchet upward (never decrease).

## Project Layout

### Key Directories

- `src/janus/` - Main Python source code (nested layout)
  - `__init__.py`, `cli.py`, `set_socrates_env.py`, `socrates.py` - Top-level package files (CLI entry, environment resolution, SOCRATES bridge)
  - `modules/` - Atmosphere-physics steps: `solve_pt.py` (the `RadConvEqm` and `MCPA_CBL` pipeline entries), `dry_adiabat_setup.py`, `moist_adjustment_H2O.py`, `spectral_planck_surface.py`, `compute_moist_adiabat.py`, `find_tropopause.py`, `set_stratosphere.py`, and the plotting-only `plot_emission_spectrum.py` / `plot_flux_balance.py`
  - `utils/` - Support code: `phys.py` (thermodynamic constants, Clausius-Clapeyron saturation vapour pressure, Planck function), `height.py` (hydrostatic height integration, inverse-square gravity), `cp_funcs.py` (NIST Shomate heat capacity), `atmosphere_column.py`, `data.py`, `logs.py`, and plotting-only `ClimateGraphicsMPL.py`

- `tests/` - Test suite. New per-source tests mirror source at `tests/<subdir>/test_<file>.py` (e.g. `tests/utils/test_phys.py`, `tests/modules/test_spectral_planck_surface.py`). The current cross-cutting files (`test_constants.py`, `test_code.py`) and pipeline-level files (`test_runaway_greenhouse.py`, `test_instellation.py`) are the documented exceptions. Shared helpers live in `tests/helpers/`.

- `tools/` - Build / utility scripts
  - `get_socrates.sh` - Clone and build SOCRATES
  - `check_test_quality.py` - AST linter (blocking on PRs)
  - `update_coverage_threshold.py` - One-way coverage ratchet (capped at 90)
  - `validate_test_structure.sh` - Module-level marker validator
  - `generate_test_badges.py` - Physics-invariant / reference-pinned coverage badges

- `docs/` - Documentation (Zensical; Diátaxis structure)
  - `Explanations/model.md`, `How-to/installation.md`, `How-to/test.md`, `Tutorials/first_run.md`, `Reference/publications.md`
  - `Validation/<file>.md` - Per-source-file inventory of `@pytest.mark.reference_pinned` tests (created when the first such test lands)

### Configuration Files

- `pyproject.toml` - Package metadata, pytest config, coverage thresholds (fast + full gates), ruff rules.
- `mkdocs.yml` - Documentation configuration (used by Zensical).
- `.github/workflows/` - CI / CD pipelines
  - `tests.yaml` - PR validation (SOCRATES build + JANUS build + pytest + coverage)
  - `docs.yaml` - Documentation build
  - `publish.yaml` - PyPI release on tag

### Entry Points

- **CLI** (`janus = "janus.cli:cli"`): `janus download spectral`, `janus download stellar`, `janus download socrates`, `janus env`.
- **Python API**: `from janus.modules.solve_pt import RadConvEqm, MCPA_CBL` for the atmosphere-structure pipeline; PROTEUS calls JANUS through its `atmos_clim` wrapper.

## Testing Standards

JANUS is scientific simulation code, so the test suite is held to physics-grade rigor. The rules below are the contract; the deep-dive (anti-happy-path patterns, discriminating-value guards, certification markers, adversarial-review trigger, `cp_mode` default flip, `satvps` argument-swap, hypothesis seed stability, pipeline solver intermediate-state assertions) lives in [`.github/.claude/rules/janus-tests.md`](.claude/rules/janus-tests.md). Read that file before editing any test file or any source file under `src/janus/**`. The two files must be kept in sync; if you change one, mirror the change in the other.

### Structure

- New per-source tests mirror source: `src/janus/<subdir>/<file>.py` -> `tests/<subdir>/test_<file>.py`. Cross-cutting tests (`test_constants.py`, `test_code.py`) and pipeline-level tests (`test_runaway_greenhouse.py`, `test_instellation.py`) are the exception, not the rule.
- Framework: `pytest` exclusively in the `tests/` directory.

### Markers and the module-level marker rule

Tier markers, with their CI surface and per-test wall-time budgets:

| Marker | What it tests | Speed budget | When CI runs it |
|---|---|---|---|
| `@pytest.mark.unit` | Python logic, heavy physics mocked | < 100 ms per test | Every PR (`unit and not skip`) |
| `@pytest.mark.smoke` | Real SOCRATES binary, low resolution, 1 step | < 30 s per test | Every PR (`smoke and not skip`) |
| `@pytest.mark.integration` | Full pipeline coupling | Minutes per test | Nightly only |
| `@pytest.mark.slow` | Full physics validation | Up to hours per test | Nightly only |
| `@pytest.mark.skip` | Placeholder, deliberately disabled | n/a | Never |

**Mandatory module-level marker** (no exceptions): every test file begins with

```python
pytestmark = [pytest.mark.<tier>, pytest.mark.timeout(<budget>)]
```

with timeouts: 30 s for unit, 60 s for smoke, 300 s for integration, 3600 s for slow. Per-function markers are additive but do not replace the module-level marker. CI runs `pytest -m "(unit or smoke) and not skip"`; tests without a tier marker are invisible to CI. The `pytest-timeout` ceiling is a defensive net against future regressions that introduce a hang.

### Physics validity

Every unit test on a **physics source** must assert at least one of the invariant families below. The six physics sources (path under `src/janus/`) are:

```
utils/phys.py                        Clausius-Clapeyron saturation vapour pressure; Planck function
utils/height.py                      hydrostatic height integration; inverse-square gravity
utils/cp_funcs.py                    NIST Shomate heat capacity cpv(vol, tmp)
modules/spectral_planck_surface.py   surface Planck emission -> Stefan-Boltzmann limit
modules/moist_adjustment_H2O.py      moist convective adjustment
modules/dry_adiabat_setup.py         Poisson dry adiabat / potential temperature
```

- **Conservation / balance**: hydrostatic balance along the column; band-integrated surface flux equal to the grey Stefan-Boltzmann limit `(1 - albedo) sigma ts**4`; column array-length / ordering integrity across an adjustment.
- **Positivity / boundedness**: T > 0, P > 0, saturation vapour pressures non-negative and finite, Planck emission non-negative, `cpv` strictly positive, mass / mole fractions in [0, 1].
- **Monotonicity or symmetry**: saturation vapour pressure increasing with T; gravity decreasing with radius; dry-adiabat T increasing with pressure; integrated height monotone increasing upward; Planck increasing with T at fixed frequency.
- **Pinned numeric value with a discrimination guard**: a closed-form value or published table entry pinned via `pytest.approx`, accompanied by explicit assertions that wrong-formula / wrong-mode / wrong-argument results would differ from the correct one by more than the tolerance.

Utility sources (`__init__.py`, `_version.py`, `set_socrates_env.py`, and plotting-only modules under `utils/` and `modules/`) are **exempt** from the physics-invariant requirement but still subject to the anti-happy-path rules.

Tag every test that asserts a physical invariant with `@pytest.mark.physics_invariant`. Per-source-file granularity: each of the six physics files needs at least one such test in its mirrored `tests/<subdir>/test_<file>.py`.

### Reference-pinned validation

Tag tests that pin against a published benchmark, an analytical limit, or a cross-implementation cross-check with `@pytest.mark.reference_pinned`. Each of the six physics files must have at least one such test. Anchors include the Smithsonian saturation-pressure tables (`satvpw`), the NIST Shomate coefficients (`cpv`), the `satvps(T0) = e0` and Poisson `T = ts` at `p = ps` identities, the grey Stefan-Boltzmann surface limit, and the pinned runaway-greenhouse / instellation OLR references. The specific anchor is chosen by the test author and recorded in `docs/Validation/<file>.md` (created when the first reference_pinned test for that source lands). The `--reference-pinned-status` mode of the linter reports the punch list of physics sources still missing a reference_pinned test.

### Anti-happy-path rules (every new test)

Every new test function MUST include:

1. **At least one edge case** (boundary value, empty input, extreme physical parameter).
2. **At least one path that exercises the error contract** (documented exception, guard return, graceful clamp). If the function under test has no validation, exercise the limit-input behavior (a subsaturated column, `p = ps` on the Poisson adiabat, `T = T0` at the Clausius-Clapeyron reference) and assert the mathematical invariant.
3. **Assertion values that are NOT trivially derivable from the implementation**: discriminating numeric pins or property-based assertions (monotonicity, conservation) preferred over point checks.

**Forbidden patterns** (flagged by `tools/check_test_quality.py`):

- Single-assert test functions.
- Standalone weak assertions (`assert result is not None`, `assert result > 0`, `assert len(result) > 0`, `assert isinstance(result, dict)`) as the only meaningful check.
- Tests with no function-level docstring.
- Tests using `==` adjacent to float literals.
- Tests asserting on a fixture's implicit default.

### Float and numerical comparison

NEVER use `==` for floats. Use `pytest.approx(val, rel=1e-5)` or `np.testing.assert_allclose(actual, expected, rtol=..., atol=...)`. For pinned numeric values, include a **discrimination guard**: a follow-up `assert` showing the wrong-formula / wrong-mode / wrong-argument value would differ from the correct one by more than the tolerance. See `janus-tests.md` Section 2 for the canonical pattern.

### Mocking discipline

- Default to `unittest.mock` for ALL external calls in unit tests: the SOCRATES binary, netCDF file I/O, spectral-file reads, `mors` calls, network.
- Mock at the narrowest scope: a specific function, not a whole module.
- A mocked physics function must return **physically plausible** values; a mock that returns `0.0` or `1.0` for everything can mask real bugs.
- NEVER mock the function under test.
- Smoke / integration / slow tiers use the real SOCRATES binary.

### Optional-dependency imports

Any test that imports an optional dependency (`hypothesis`, `mors`) MUST call `pytest.importorskip('<dep>')` at module top. `mors` (the `fwl-mors` stellar-evolution package) is pulled in transitively by the helpers under `tests/helpers/`; the `pip install --no-deps` CI image lacks it, so any test file importing it will otherwise fail to collect.

### Module-level constants and `monkeypatch`

When the source under test reads an env var into a module-level constant at import time (e.g. `FWL_DATA_DIR` in `janus.utils.data`, `RAD_DIR` in `janus.set_socrates_env`), `monkeypatch.setenv` alone is not sufficient: the constant is frozen at import. Patch the constant directly with `monkeypatch.setattr('janus.utils.data.FWL_DATA_DIR', tmp_path, raising=False)`, and patch both the env var and the constant when downstream code re-reads the env var.

### Voice rule for test artifacts

The repo-wide voice rule (zero AI-process disclosure in any public artifact) applies to test code with the same strictness as to source. Scope: test-skip reasons, test-file / function docstrings, test-function / class names, parametrize ids, log-capture assertions, **commit messages on test-touching commits, pull-request titles and bodies on test-touching PRs**, GitHub Actions job / step names, inline `src/janus/**` comments, and shipped log strings. Out of scope: the rule documents themselves (this file, `janus-tests.md`, `janus-code-review.md`, `docs/How-to/test.md`) may legitimately name the procedures they define, though the public docs page applies the rule silently and does not enumerate the banned phrases.

Banned phrases inside in-scope artifacts: "audit", "review pass", "adversarial review", AI-roadmap labels (`Phase X`, `Stage X.Y`, `Iteration N`, `T1.x`, `Group A/B/C/D` when AI-organized), `claude-config/...` paths, "Generated with Claude", AI-tool names, em-dashes, en-dashes (except bibliographic page ranges).

Write the OUTCOME, never the PROCESS.

### Speed and determinism

- Unit tests: < 100 ms wall-time each.
- Aggressively mock the SOCRATES binary and file I/O in unit tests.
- Set seeds for any randomness: `np.random.seed(42)`, `random.seed(42)`. Hypothesis tests use `@settings(derandomize=True)` or an explicit `--hypothesis-seed` (see `janus-tests.md` Section 16).
- Use `tmp_path` (pytest fixture) for temporary files.

### Documentation per test

- File-level docstring: name the source under test, list the invariants and contract clauses the file exercises, link to `docs/How-to/test.md`.
- Function-level docstring: state the physical scenario or contract clause being verified. Required (lint-enforced).
- Inline comments: explain **why** a specific input range was chosen.

### Independent review trigger

A pull request that adds or substantially modifies > 50 lines of test code across all its commits triggers an independent review pass before merge. The denominator is PR-level (`git diff origin/main...HEAD -- 'tests/**'`); splitting into many sub-50-line commits does not dodge the trigger. The reviewer cites the anti-happy-path rule, the discrimination-guard requirement, and the physics-invariant tier.

### Tooling

- Validate test structure: `bash tools/validate_test_structure.sh`
- Test-quality lint: `python tools/check_test_quality.py --check`
- Baseline regeneration (after a deliberate sweep): `python tools/check_test_quality.py --baseline`
- Reference-pinned status: `python tools/check_test_quality.py --reference-pinned-status`
- Coverage ratchet (one-way, capped at 90): `python tools/update_coverage_threshold.py`
- Test badges: `python tools/generate_test_badges.py`
- Format: `ruff format src/ tests/`
- Lint: `ruff check src/ tests/`

### Coverage architecture

JANUS uses two gates with explicit sub-targets:

| Gate | Tests included | Target | Enforced |
|---|---|---|---|
| Fast gate (`tool.janus.coverage_fast.fail_under`) | unit + smoke | Ratcheting toward **90%** | Every PR |
| Full gate (`tool.coverage.report.fail_under`) | unit + smoke + integration + slow | **90%** | Nightly |

Both gates ratchet toward 90, capped at 90 (`tools/update_coverage_threshold.py` enforces the 90.0 ecosystem ceiling); neither may be manually decreased. The CI guard in `tests.yaml` rejects any PR that lowers `[tool.coverage.report].fail_under` below `min(base_ref, 90.0)`.

## Safety & Determinism

- **Randomness**: explicitly set seeds in tests.
- **Files**: do not generate tests that produce large netCDF output files; use `tempfile` or mocks.

## Code Quality

**Style** (enforced by ruff):

- Line length < 96 chars.
- Variables / functions: `snake_case`.
- Constants: `UPPER_CASE`.
- Type hints: standard Python.
- Docstrings: brief descriptions of physical scenarios.

**Pre-commit**: runs `ruff check --fix` automatically. Fix issues before committing.

## Common Workflows

### Making a Code Change

1. **Create branch**: `git checkout -b <initials>/<short-description>`.
2. **Make changes** in `src/janus/`.
3. **Write / update tests** in `tests/<subdir>/test_<file>.py` (mirror structure).
4. **Run tests locally**: `pytest -m "(unit or smoke) and not skip"` (needs `RAD_DIR` and `FWL_DATA` set).
5. **Check coverage**: `coverage run -m pytest && coverage report`.
6. **Lint**: `ruff check --fix src/ tests/ && ruff format src/ tests/`.
7. **Validate structure**: `bash tools/validate_test_structure.sh`.
8. **Test quality**: `python tools/check_test_quality.py --check`.
9. **Commit**: plain-language subject, first-person voice, no AI-process disclosure.
10. **Push**: CI runs automatically on PR.

### Adding a New Physics Source

1. Create `src/janus/<subdir>/<file>.py`.
2. Create `tests/<subdir>/test_<file>.py` with module-level `pytestmark`.
3. Add at least one `@pytest.mark.physics_invariant` test asserting one of the four invariant families.
4. Plan a `@pytest.mark.reference_pinned` test (anchor: paper, analytical limit, or cross-check); create `docs/Validation/<file>.md` when it lands.
5. Run the full PR checks locally.

### Debugging Test Failures

```bash
pytest -v --showlocals                                # Verbose with local variables
pytest -x                                             # Stop at first failure
pytest tests/<subdir>/test_<file>.py::test_function   # Run specific test
pytest --pdb                                          # Drop into debugger on failure
```

## Documentation References

- **Testing rules**: `.github/.claude/rules/janus-tests.md`, `.github/.claude/rules/janus-code-review.md`
- **Test how-to**: `docs/How-to/test.md`
- **Installation**: `docs/How-to/installation.md`
- **Model concepts**: `docs/Explanations/model.md`
- **First run**: `docs/Tutorials/first_run.md`

## Project memory and session learnings

Session-specific knowledge (debugging logs, design rationale, sprint focus, ADR drafts) lives outside this repository, in the Claude memory tree under `~/.claude/projects/<project>/memory/`. The memory tree is per-user, sync-ready across machines, and not exposed in public commit history.

What still lives in this repository:

- Architectural decisions that affect every contributor: this file (`.github/copilot-instructions.md`).
- Test and review rules: `.github/.claude/rules/janus-tests.md` and `.github/.claude/rules/janus-code-review.md`.
- Per-PR rationale: PR descriptions.
- Per-commit rationale: commit messages.
- Module-level scientific validation: `docs/Validation/<file>.md` (created when the first `@pytest.mark.reference_pinned` test for that source lands).

Do not introduce a new in-repo "memory" or "decisions log" file. The four channels above are the contract.

---

## Quick Reference

```bash
# Setup
bash tools/get_socrates.sh /path/to/socrates
export RAD_DIR=/path/to/socrates
export FWL_DATA=/path/to/fwl_data
pip install -e ".[develop]"
pre-commit install -f

# Test
pytest -m "(unit or smoke) and not skip"
coverage run -m pytest && coverage report

# Lint
ruff check --fix src/ tests/
ruff format src/ tests/

# Validate
bash tools/validate_test_structure.sh
python tools/check_test_quality.py --check

# Serve docs locally
pip install -e '.[docs]'
zensical serve
```

**Remember**: Trust these instructions. Only search if information is incomplete or found to be in error.

---

> **⚠️ FILE SIZE LIMIT: This file must stay below 750 lines.** Enforced by pre-commit hook (`tools/check_file_sizes.sh`). File located at `.github/copilot-instructions.md`.
