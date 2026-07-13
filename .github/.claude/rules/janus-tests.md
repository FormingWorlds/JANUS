---
description: JANUS test quality deep-dive. Anti-happy-path patterns, discriminating-value guards, physics-invariant tiering, validation certification markers, adversarial-review trigger. Extends the Testing Standards section in `.github/copilot-instructions.md`.
---

# JANUS Test Quality Rules

This file is the canonical deep-dive on test quality. The high-level summary lives in [`.github/copilot-instructions.md`](../../copilot-instructions.md) under "Testing Standards". The two files MUST stay in sync. If you change one, mirror the change in the other.

> **Discovery note.** JANUS keeps its Claude-Code rule files under `.github/.claude/rules/` (not the conventional repo-root `.claude/`) so they can be tracked in git and shared across collaborators. Claude does NOT auto-discover them at this path; the repo-root `CLAUDE.md` (symlinked to `.github/copilot-instructions.md`) names this file and `janus-code-review.md` explicitly so AI tooling and human readers know to load them. **When opening or editing any file under `tests/**` or `src/janus/**`, read this file first.**

Sister rule files:

- [`.github/copilot-instructions.md`](../../copilot-instructions.md): high-level rules, applied repo-wide.
- [`.github/.claude/rules/janus-code-review.md`](janus-code-review.md): review-pass gate, domain-aware code review (spectral unit boundaries, per-band normalization, dry-adiabat exponent safety, PROTEUS-coupling contract). Test-marker discipline lives there too.

JANUS is a 1D convective / radiative atmosphere module and the test suite is held to physics-grade rigor. Tests exist to catch real bugs. A test that asserts the wrong thing, or that passes for the wrong reason, is worse than no test because it generates false confidence. The rules below codify what "real test" means here.

---

## 1. Anti-happy-path rules (every new test)

Every new test function MUST include:

1. **At least one edge case**: a boundary value (`p = ps` on the dry adiabat, `T = T0` at the Clausius-Clapeyron reference point, `albedo_s = 0` or `1` for the surface Planck flux, a subsaturated column for the moist adjustment, `nu -> 0` or very large `T` for the Planck function), an empty input, or an extreme physical parameter.
2. **At least one path that exercises the error contract**:
   - If the function under test has documented validation (raises on non-positive temperature, refuses an unknown `cp_mode`), test that the error fires AND that no side effect ran.
   - If the function has no validation (closed-form thermodynamics: Clausius-Clapeyron, Poisson adiabat, inverse-square gravity), exercise the **limit-input behavior** (a subsaturated column is a degenerate fixed point of the moist adjustment; `p = ps` is the identity point of the Poisson relation) and assert the corresponding mathematical invariant.
   - "No validation in source therefore no error test" is not an exemption; the limit-input substitute is.
3. **Assertion values NOT trivially derivable from the implementation**: discriminating numeric pins (see Section 2 below) or property-based assertions (monotonicity, conservation, symmetry, boundedness).

### Forbidden patterns

These are flagged by `tools/check_test_quality.py` and rejected at PR time.

- **Single-assert test functions**. Two or more assertions per test; the second usually pins the invariant the first hand-waves over. Exception: a single assertion of a hard-fail invariant (hydrostatic closure within `1e-9`) is acceptable if the test is the only test of that invariant in the file.
- **Weak assertions when they stand alone as the sole meaningful check in the test.** The shapes are:
  - `assert result is not None`
  - `assert result > 0`
  - `assert len(result) > 0`
  - `assert isinstance(result, dict)`
  - `assert result is None` where the function returns `None` implicitly

  Required carve-out: the three-class discrimination guard (Section 2) uses `assert val > 0` as the sign-error guard and `assert lo < val < hi` as the scale-error guard alongside a primary `pytest.approx(...)` pin. Those secondary lines look like weak assertions in isolation; they are NOT flagged when paired with a stronger primary assertion in the same test. The linter applies the carve-out automatically: weak shapes are flagged only when the test has exactly one `assert` statement (`len(asserts) == 1`) and that assertion is itself the weak shape.
- **Tests with no function-level docstring**. The docstring states which physical scenario or contract clause is being verified.
- **`==` adjacent to a float literal**. Use `pytest.approx(val, rel=...)` or `np.testing.assert_allclose(actual, expected, rtol=..., atol=...)`. Comparing two floats with `==` is a known flake source even for "exact" identities like 0.0 (-0.0 vs +0.0, NaN propagation).
- **Tests asserting on a fixture's implicit default**: e.g. `assert fixture_returning_none() is None`. This is trivially true. Delete the test; do not strengthen it by adding more `is None` assertions.

---

## 2. Discriminating test values

The test contract is: a regression that introduces a plausible bug must fail the test. "Plausible bug" means off-by-one exponent, wrong sign, swapped factor of 2, missing factor of pi, dimensionally-wrong unit, **wrong-mode / wrong-argument selection** (`cp_mode = "constant"` vs `"T-dependent"`, a `satvps` molecular-weight / latent-heat argument swap). Pick input values where the wrong-formula result is far from the correct one.

### Bad / good examples

| Pattern | Bad (any-exponent-passes) | Good (discriminates) |
|---|---|---|
| Clausius-Clapeyron `satvps(T, T0, e0, M, L)` | Test at `T = T0` only (the exponential collapses to `e0`; any `L`, any `M` passes) | Test at `T` well below `T0` (e.g. `T = T0 - 30`) so the exponential factor is far from 1; a factor slip in `L` or `M` (a doubled or halved `L M / Rstar` coefficient) moves the result outside tolerance. A clean `M` / `L` swap is symmetric in the coefficient and cannot be caught by a value pin (see Section 16); it is a call-site review item. |
| Planck function `B(nu, T)` | Test at one `(nu, T)` pair | Test at `T = 300 K` AND `T = 1500 K` so the `T**3` vs `T**4` slope is resolved; a dropped factor changes the ratio |
| Poisson dry adiabat `T = ts (p / ps)**Rcp` | Test at `p = ps` (identity point, `T = ts` for any `Rcp`) | Test at `p = ps / 10` where `(p/ps)**Rcp` discriminates the exponent sign and magnitude |
| Heat capacity `cpv(vol, tmp, cp_mode=...)` | Test one `cp_mode` at one `T` | Test `cp_mode = "constant"` vs `"T-dependent"` at high `T` where the Shomate curvature separates the two paths beyond tolerance |

### Discrimination guard (REQUIRED for pinned-value tests)

When a test pins a numeric value, include explicit assertions that the wrong-formula result would differ from the correct one for **each plausible bug class**. "Each plausible bug class" means at minimum:

1. **Exponent or factor error** (off-by-one exponent, missing factor of 2 / pi). `abs(val - wrong_value)` discriminates; a monotonicity check across two inputs also catches it.
2. **Sign error** (`-x` vs `+x`). `abs()` hides this; assert the sign explicitly with `val > 0` or `val < 0`.
3. **Unit-conversion error** (Pa vs bar, dyn/cm^2 vs Pa, K vs C, cm^-1 vs m^-1). Pin the absolute scale with the unit named in the comment.
4. **Wrong-mode / wrong-argument selection** (`cp_mode = "constant"` vs `"T-dependent"`; a `satvps` argument swap; a per-band vs per-wavenumber normalization slip). When the function dispatches by a mode string or takes several same-typed float arguments, the discrimination guard MUST include a value that distinguishes the chosen path from a sibling path.

**Carve-out for conservation-style invariants.** When the primary assertion IS a conservation closure (hydrostatic balance, energy balance, band-integrated flux equal to the Stefan-Boltzmann limit), the equality form `lhs == pytest.approx(rhs)` already discriminates exponent / factor errors by construction. The exponent guard is satisfied by the conservation equality itself; sign and scale guards remain mandatory.

Canonical pattern:

```python
def test_satvpw_at_steam_point_matches_one_atmosphere():
    """Pin satvpw at the 373.16 K steam point against the 1 atm reference."""
    val = phys.satvpw(373.16)
    expected = 101514.5  # Pa; Smithsonian liquid-water curve, 0.19% above 1 atm
    assert val == pytest.approx(expected, rel=1e-4)
    # Unit / scale guard: the Smithsonian formula computes in dyn/cm^2 and the
    # source multiplies by 0.1 to reach Pa. A forgotten conversion leaves the
    # value near 1.015e6 (ten times too large); pin the magnitude.
    assert 9.0e4 < val < 1.2e5
    # Sign / positivity guard: a saturation vapour pressure is strictly positive.
    assert val > 0
    # Exponent guard: satvpw is monotone increasing in T, so a value 20 K cooler
    # must be strictly smaller. A flipped temperature ratio would invert this.
    assert phys.satvpw(353.16) < val
```

The guard lines are mandatory whenever the test's primary assertion is a `pytest.approx` against a hand-calculated or published value. For functions that dispatch by mode (`cpv` with `cp_mode`), the guard MUST additionally pin a value that distinguishes the selected mode from its sibling (Section 2 rule 4). Property-based assertions (monotonicity, conservation, symmetry) do not need a separate guard because they are already discriminating across the input space.

---

## 3. Physics-invariant assertions (tiered)

### When required

Every unit test on a **physics source** must assert at least one of the four invariants below. The physics sources, by their path under `src/janus/`, are:

```
utils/phys.py                        thermodynamic constants; Clausius-Clapeyron
                                     saturation vapour pressure (satvps, satvpw,
                                     satvpi); Planck function B(nu, T)
utils/height.py                      hydrostatic height integration; inverse-square
                                     gravity g = G m / r**2
utils/cp_funcs.py                    NIST Shomate heat capacity cpv(vol, tmp)
modules/spectral_planck_surface.py   surface Planck emission, band-integrated to the
                                     Stefan-Boltzmann limit (1 - albedo) sigma ts**4
modules/moist_adjustment_H2O.py      moist convective adjustment; subsaturated column
                                     returns zero temperature tendency
modules/dry_adiabat_setup.py         Poisson dry adiabat / potential temperature
                                     T = ts (p / ps)**(R/cp)
```

Per-source-file granularity: each of the six physics files needs at least one `@pytest.mark.physics_invariant` test and at least one `@pytest.mark.reference_pinned` test in its mirrored `tests/<subdir>/test_<file>.py`. Granularity is per source file, not per directory: `utils/phys.py` and `utils/height.py` each need their own pinned test even though both live under `utils/`.

Utility sources are exempt from the physics-invariant requirement but still subject to all anti-happy-path rules:

```
src/janus/__init__.py           (re-exports)
src/janus/_version.py           (version string)
src/janus/set_socrates_env.py   (environment resolution, no physical quantity)
plotting-only modules under utils/ and modules/ that solely render figures
   (e.g. utils/ClimateGraphicsMPL.py, modules/plot_emission_spectrum.py,
    modules/plot_flux_balance.py)
```

### The four invariant families

1. **Conservation / balance**
   - Hydrostatic balance: the integrated height step is consistent with `dP = -rho g dz` along the column.
   - Energy balance: the band-integrated surface flux equals the grey Stefan-Boltzmann limit `(1 - albedo) sigma ts**4` as the spectral grid widens.
   - Column integrity: pressure, temperature, and height arrays stay the same length and stay ordered after an adjustment step.
2. **Positivity / boundedness**
   - `T > 0` Kelvin everywhere, `P > 0` Pa everywhere.
   - Saturation vapour pressures (`satvps`, `satvpw`, `satvpi`) non-negative and finite (no `nan` / `inf`).
   - Planck emission `B(nu, T)` non-negative; band-integrated surface flux non-negative.
   - Heat capacity `cpv` strictly positive; mass / mole fractions in `[0, 1]`.
3. **Monotonicity or symmetry**
   - Saturation vapour pressure increasing with temperature (Clausius-Clapeyron).
   - Gravity decreasing with radius (inverse-square law).
   - Dry-adiabat temperature increasing with pressure at fixed `ts`, `ps`.
   - Integrated height monotone increasing upward through the column.
   - Planck `B(nu, T)` increasing with `T` at fixed `nu`.
4. **Pinned numeric value with a discrimination guard**: see Section 2. Acceptable as the sole invariant when a closed-form result or published table value is the contract (`satvps(T0) = e0`; `satvpw(373.16) ~= 101514.5 Pa`; the runaway-greenhouse OLR references).

Property-based assertions (monotonicity, conservation, symmetry, boundedness) are preferred over point-value pins when both are possible. They hold for any valid input and so catch bugs across the entire input space.

### Validation certification markers

Two markers track validation quality independently of line coverage:

- **`@pytest.mark.physics_invariant`** -- this test asserts at least one of the four invariants. Tag every qualifying test in a physics-source test file. `tools/check_test_quality.py` warns when a physics-source test asserts no invariant and is not tagged.
- **`@pytest.mark.reference_pinned`** -- this test pins behavior against a **published benchmark** (paper, figure, table; cite explicitly in the test docstring, e.g. the Smithsonian Meteorological Tables for `satvpw`, the NIST Shomate coefficients for `cpv`), an **analytical limit** (`satvps(T0) = e0` identity, the Poisson identity `T = ts` at `p = ps`, the grey Stefan-Boltzmann surface limit), or a **cross-implementation cross-check** (the JANUS OLR against the pinned runaway-greenhouse and instellation reference values in the existing suite).
  - **Per-source-file**: each of the six physics source files must have at least one `reference_pinned` test in its mirrored `tests/<subdir>/test_<file>.py`. Anchor type is one of {published benchmark, analytical limit, cross-implementation cross-check}; the specific paper or limit is chosen by the test author and recorded in `docs/Validation/<file>.md`.
  - **Tracking**: each physics source gets a page at `docs/Validation/<file>.md`, created when the first reference_pinned test for that source lands. The page records: the source under test, the reference cited, the test ids carrying the marker, and the date of last comparison against the source.
  - **Status report**: `python tools/check_test_quality.py --reference-pinned-status` reports the physics source files missing a `reference_pinned` test. This is the punch list for follow-up validation work.

Both markers are registered in `pyproject.toml` under `[tool.pytest.ini_options] markers`. They do not gate CI on their own; their coverage is a separate KPI surfaced in the PR summary comment.

---

## 4. Mocking discipline

- Default to `unittest.mock` for ALL external calls in unit tests: the SOCRATES radiative-transfer binary, netCDF file I/O, spectral-file reads, `mors` stellar-spectrum calls, HTTP / OSF downloads, subprocess.
- Mock at the **narrowest scope**: patch the specific function (`unittest.mock.patch('janus.modules.solve_pt.some_helper')`), not the whole module.
- A mocked physics function MUST return **physically plausible** values. A mock that returns `0.0` or `1.0` for everything will mask sign / clamp / fallback bugs.
- NEVER mock the function under test. If you're tempted to, the test is asking the wrong question.
- Smoke tests run the real SOCRATES binary on a low-resolution atmosphere over a single step; integration and slow tests run the full pipeline (`RadConvEqm`, `MCPA_CBL`). The rules in this file still apply to those tiers, but the mocking discipline is relaxed because the real call is the contract.

---

## 5. Optional-dependency imports

Any test that imports an optional dependency MUST call `pytest.importorskip` at module top so `pip install --no-deps` CI runs do not fail collection:

```python
import pytest

hypothesis = pytest.importorskip('hypothesis')
# ... or for a helper that pulls in the stellar-evolution dependency:
pytest.importorskip('mors')
```

Optional deps recognized by the linter (`OPTIONAL_DEPS` constant in `tools/check_test_quality.py`):

- `hypothesis` (used in property-based / fuzz tests; lives in `[develop]` extras).
- `mors` (the `fwl-mors` stellar-evolution package, pulled in transitively by the test helpers under `tests/helpers/`; absent from the `pip install --no-deps` PR image, so it must always be guarded if imported into a test file).

The lint script enforces this. Rule key `missing_importorskip`: any module-top `import <optional_dep>` or `from <optional_dep> import ...` that is not preceded by a module-scope `pytest.importorskip('<optional_dep>')` is flagged.

---

## 6. Module-level constants and `monkeypatch`

When the source under test reads an env var or a class-level default into a **module-level constant** at import time, `monkeypatch.setenv` alone is not sufficient. The constant is frozen at the original import.

Pattern:

```python
# Source: src/janus/set_socrates_env.py
RAD_DIR = os.environ.get('RAD_DIR', ...)  # frozen at import
```

```python
# Test (wrong):
monkeypatch.setenv('RAD_DIR', str(tmp_path))          # too late; constant already bound

# Test (right):
monkeypatch.setattr('janus.set_socrates_env.RAD_DIR', tmp_path, raising=False)
```

The same trap applies to `FWL_DATA`, which `src/janus/utils/data.py` resolves into the module constant `FWL_DATA_DIR` at import; patch `janus.utils.data.FWL_DATA_DIR` directly, not just `setenv('FWL_DATA', ...)`. A related pattern is a function-default flip such as `cpv(vol, tmp, cp_mode='constant')`: mutating `cpv.__defaults__` to force a mode is **single-threaded only**; never apply it inside library code (it has process-wide visibility and races under pytest-xdist).

When in doubt, do both the env-var monkeypatch and the constant monkeypatch. The lint script does NOT currently flag this pattern (it would require source-side analysis to know which constants are env-derived); this is a discipline rule enforced via the >50 LOC review trigger and the recurring-trap table in Section 16.

---

## 7. Marker discipline and timeouts

### Module-level marker is mandatory

Every test file MUST begin with:

```python
import pytest

pytestmark = [pytest.mark.<tier>, pytest.mark.timeout(<budget>)]
```

with budgets:

- `unit` -> `timeout(30)` (target wall-time per test is `< 100 ms`; the 30 s cap is a defensive net).
- `smoke` -> `timeout(60)` (target `< 30 s`).
- `integration` -> `timeout(300)`.
- `slow` -> `timeout(3600)`.

PR CI runs `pytest -m "(unit or smoke) and not skip"`. Tests without the tier marker are invisible to CI and shipped untested. The lint script blocks any file missing the module-level `pytestmark`.

### Per-function markers

Per-function `@pytest.mark.<tier>` markers are **additive**, not a replacement for the module-level marker. They are useful when a file's tests span multiple tiers (rare; prefer separate files).

### Timeout is a safety net, not a target

The `timeout` ceiling exists so a future regression that introduces a hang (real SOCRATES call, infinite convective loop, network retry) surfaces as a specific-test failure rather than a generic job timeout. Current unit test wall times are 100x below the ceiling; if you find yourself needing the full 30 s for a unit test, something has gone wrong and you should reduce scope or move the test to a slower tier.

---

## 8. Float and numerical comparison

- NEVER use `==` for floats. Use `pytest.approx(val, rel=1e-5)` or `np.testing.assert_allclose(actual, expected, rtol=..., atol=...)`.
- State the tolerance rationale in a comment when the choice is non-obvious. E.g. "`rtol=1e-5` because the runaway-greenhouse OLR references are pinned to five significant figures in the existing suite".
- For pinned numeric values, include a **discrimination guard** (Section 2).
- For property-based assertions, use `pytest.approx` against the exact symbolic relation, with the tightest tolerance the implementation can hit (typically `rel=1e-12` for closed-form algebra such as the Poisson identity; looser for pipeline outputs that route through SOCRATES).

---

## 9. Voice rule for test artifacts

The repo-wide voice rule (zero AI-process disclosure in any public artifact) applies to test code with the same strictness as to source. The voice rule is **scoped** to public artifacts other contributors and external readers see; it does NOT apply to the rule documents themselves (this file, `janus-code-review.md`, `copilot-instructions.md`), which legitimately name the procedures they define.

In scope (the voice rule is BANNED here):

- Test-skip reasons (`@pytest.mark.skip(reason='...')`).
- Test-file docstrings.
- Test-function and test-class names.
- Test-function docstrings.
- Parametrize ids (`@pytest.mark.parametrize('name', [...], ids=[...])`).
- Log-capture assertions (regex against `caplog.records`).
- Commit messages on test-touching commits (subject AND body).
- **Pull-request titles and bodies on test-touching PRs**.
- GitHub Actions job names and step names that ship to the PR Checks tab.
- Inline source comments and docstrings on `src/janus/**`.
- Log strings that ship with the repo.
- **All public-facing documentation** (anything under `docs/`, the repo README, CONTRIBUTING.md, tutorials, wiki pages). Public docs apply the rule silently; they do NOT enumerate the banned phrases, name the voice rule, advertise the existence of `.github/.claude/` rule infrastructure, or cross-reference `.github/.claude/rules/*.md` files. A user docs page that describes the testing infrastructure must do so without naming the AI-process rules that produced it.

Out of scope (these may NAME the procedures they define):

- This file (`janus-tests.md`).
- `janus-code-review.md`.
- `copilot-instructions.md`.

Banned phrases inside the in-scope artifacts: "audit", "review pass", "adversarial review", "Phase X" (when "X" is an AI-organized roadmap label, not a real project phase), "T1.x", "Group A/B/C/D" (when AI-organized work groups), `claude-config/...` paths, "Generated with Claude", AI-tool names, em-dashes, en-dashes (except in bibliographic page ranges within citations), process meta-commentary ("after careful analysis").

Write the OUTCOME (what the test verifies; what the PR achieves) never the PROCESS (how the rule was derived; which review caught what). First-person Tim voice. Going-forward only, no history rewrite.

---

## 10. Fixture and parameter conventions

- Use SI units in test parameters unless the function under test explicitly expects other units (bar, wavenumber in cm^-1). Name the unit in the parameter comment when it is not SI.
- Use `@pytest.mark.parametrize` when the same logic spans multiple physical regimes (sub-runaway, Simpson-Nakajima plateau, post-runaway; Earth-like vs hot super-Earth; near-side vs far-side instellation). Each parametrize id must read like a physical scenario, not a tuple of numbers.
- Set seeds for any randomness:
  ```python
  np.random.seed(42)
  random.seed(42)
  ```
  Hypothesis tests use `@settings(derandomize=True)` or an explicit `--hypothesis-seed` to keep replays stable across versions (see Section 16 trap).
- Use `tmp_path` (pytest fixture) for temporary files. Do not produce large netCDF outputs in the test path.

---

## 11. Documentation per test

- **File-level docstring**: name the source file under test (`Tests for src/janus/<subdir>/<file>.py`), list the invariants and contract clauses the file exercises, link to `docs/How-to/test.md`. Required.
- **Function-level docstring**: state the physical scenario or contract clause in plain language. Required (lint-enforced).
- **Inline comments**: explain **why** a specific input range was chosen ("`T = 300 K` and `T = 1500 K` so the Planck `T**3` vs `T**4` slope difference is resolved well above tolerance").

---

## 12. Naming

- Test names describe behavior, not the called function: `test_satvpw_monotonic_with_temperature`, NOT `test_satvpw`.
- Test names use snake_case and read as full sentences.
- Group related tests in classes (`class TestSatVP:`) when they share setup; use the class to thread a single fixture through several scenarios.
- Test file names mirror source: `src/janus/<subdir>/<file>.py` -> `tests/<subdir>/test_<file>.py`. Documented exceptions to the strict mirror:
  - **Cross-cutting utility / constant tests** (`tests/test_constants.py`, `tests/test_code.py`): tests that span multiple source files (physical constants, planet database, logger, CLI) or test package-level concerns.
  - **Pipeline-level tests** (`tests/test_runaway_greenhouse.py`, `tests/test_instellation.py`): end-to-end runs of `RadConvEqm` and `MCPA_CBL` that pin OLR / flux references against the whole coupled column rather than a single source file. These are the smoke / integration anchors and their reference values feed the `reference_pinned` inventory.
  - **Topical sub-files of a large physics source**: when a physics source exceeds ~500 LOC and its tests split into independent topics that would not benefit from consolidation, topical sub-files are acceptable alongside the primary `tests/<subdir>/test_<file>.py`. The primary file must still exist and carry at least one `reference_pinned` and one `physics_invariant` test; the topical sub-files cover the remaining surface.

---

## 13. Adversarial review trigger

A pull request that adds or substantially modifies **> 50 lines of test code across all its commits** triggers an independent review pass before merge. This is a discipline rule, not CI-automated: the author runs the review pass via a `code-reviewer` agent before pushing the final test-touching commit. The denominator is PR-level, not per-commit: `git diff origin/main...HEAD -- 'tests/**'` is the source of truth. Splitting one large change into 49 + 49 + 49 line commits does NOT dodge the trigger.

The reviewer's mandate:

- Cite the anti-happy-path rule (Section 1) and the discrimination-guard requirement (Section 2).
- Flag single-assert tests, weak `is not None` patterns, missing module-level marker, missing `physics_invariant` tag on a physics-source test, missing `reference_pinned` tag on a per-source benchmark test, dead tests (passes for the wrong reason), tests that mock the function under test.
- Verify discriminating values: re-derive the expected value from a plausible wrong formula (a swapped `satvps` argument, a dropped `/1000.0` band-width factor, a flipped Poisson exponent) and assert the test fails with that wrong formula.
- Verify physics-source coverage of the four invariants: which of the four does this test exercise? If none, why is the test in `tests/<subdir>/test_<physics_file>.py`?

The reviewer is a separate process from the test author. For Claude-Code workflow this means spawning a `proteus-review` skill or a `code-reviewer` agent with the test files in scope; the review must complete and surface findings before the test commit is pushed.

The reviewer's findings are addressed in a follow-up commit (not amended into the test commit). The follow-up subject line is in plain language describing the OUTCOME ("sharpen satvp assertions to pin the Clausius-Clapeyron exponent sign", NOT "address review findings").

---

## 14. Tooling

The repo provides:

- `bash tools/validate_test_structure.sh` -- structural check (marker presence, file naming, source-to-test mirror).
- `python tools/check_test_quality.py --check` -- CI mode: AST scan for the forbidden patterns in Section 1 and the marker requirement in Section 7. Fails the PR if violations exceed the baseline.
- `python tools/check_test_quality.py --baseline` -- after a deliberate sweep, regenerates `tools/test_quality_baseline.json`. Only run when you have intentionally reduced violations. Setting `JANUS_TEST_QUALITY_ALLOW_REGRESS=1` lets a single run pass above the baseline for an intentional, separately-tracked change; do not set it to silence an accidental regression.
- `python tools/check_test_quality.py --reference-pinned-status` -- prints physics source files missing a `reference_pinned` test.
- `python tools/update_coverage_threshold.py` -- ratchet the fast PR gate upward when measured coverage exceeds the current `fail_under`. Capped at the 90% ecosystem ceiling.
- `python tools/generate_test_badges.py` -- regenerate the physics-invariant and reference-pinned coverage badges surfaced in the PR summary.
- `ruff check src/ tests/` and `ruff format src/ tests/` -- run before commit.

The lint script is wired into PR CI (`tests.yaml`). The step runs in **blocking** mode: any regression above the baseline fails the PR.

---

## 15. Coverage strategy (operator's view)

JANUS uses two coverage gates with explicit sub-targets. The fast gate is for PR cycle time; the full nightly gate is the long-running KPI.

| Gate | Tests | Target | When |
|---|---|---|---|
| Fast gate (`tool.janus.coverage_fast`) | unit + smoke | ratcheting toward **90%** (PROTEUS-ecosystem ceiling) | Every PR |
| Full gate (`tool.coverage.report`) | unit + smoke + integration + slow | **90%** | Nightly |

The ratchet is one-way (`tools/update_coverage_threshold.py`), capped at 90%. Never manually decrease the threshold. The CI guard in `.github/workflows/tests.yaml` rejects any PR that lowers `[tool.coverage.report].fail_under` below `min(base_ref, 90.0)`.

What this means for adding tests:

- A new closed-form helper in a utility module: a unit test is sufficient.
- A new function in a physics source: a unit test (counts toward both gates), plus a `physics_invariant` tag if it qualifies. If the function feeds a published benchmark, a `reference_pinned` test goes with it (counts toward the per-source-file inventory in `docs/Validation/<file>.md`).
- A new end-to-end column run: a smoke-tier test that calls the real SOCRATES binary on a low-resolution atmosphere and pins against a reference OLR / flux.

---

## 16. Failure modes to recognize on review

These are real patterns that have shipped in the past. The lint script catches some of them mechanically; reviewers catch the rest.

| Pattern | Example | Why it slipped | Fix |
|---|---|---|---|
| **`cp_mode` default flip propagation** | A test pins `cpv(vol, T)` at high `T` against the default `cp_mode = "constant"`. The default flips to `"T-dependent"` (Shomate); the constant-mode value no longer matches but the test still passes against the stale reference. | The reference value is tied to the mode choice but the test does not name the mode in its docstring or assert against a mode-discriminating value. | Discrimination guard (Section 2 rule 4) must include the alternative-mode value so a regression that silently dispatches to the wrong mode fails the test. Cite the mode AND the Shomate source in the docstring. |
| **`satvps` argument swap** | `satvps(T, T0, e0, MolecularWeight, LatentHeat)` takes several same-typed floats; a caller swaps `MolecularWeight` and `LatentHeat`. | The exponent coefficient is `LatentHeat * MolecularWeight / Rstar`, symmetric under the swap, so the swapped call returns an identical value at every `T`. No value pin can catch a clean swap; it is invisible to the test suite and must be caught at the call site. | Catch the swap during code review by confirming argument order at each call site. A value pin at `T != T0` still earns its place, but for a different bug: it pins the `L M / Rstar` coefficient, so it catches an individual factor slip in `L` or `M` (a doubled or halved coefficient) that the `T = T0` identity cannot. |
| **Hypothesis seed and version stability** | `@given(...)` test passes on hypothesis 6.0 with the default seed strategy; on hypothesis 6.100 the strategy produces a different sequence and the test surfaces a previously-hidden flake or stops covering the previously-hidden bug class. | Hypothesis seed semantics changed between minor releases; the test author relied on implicit determinism. | Add `pytest.importorskip('hypothesis')` at module top. Use `@settings(derandomize=True)` or pass `--hypothesis-seed=<fixed>` in CI. Document the chosen seed in the test docstring. |
| **Pipeline solver intermediate-state types** | The moist adjustment or `solve_pt` loop produces a `nan` or negative intermediate temperature; the silently-clamped final OLR looks plausible but is wrong. The unit test only checks the final OLR. | The output check is too late; the bug lives in an intermediate convective step. | Tests of the adjustment / solver path must assert intermediate state is finite and physical at each step: `assert np.all(np.isfinite(profile))` and `assert np.all(profile > 0)`. The source-side defense (clamp against `nan` / negative T) belongs in the module, but tests verify the defense actually fires. |
| **Silent skip in helper** | `def _band_for(idx): ...; if actual is None: continue` masks a broken band lookup | Helper hides a real failure as a no-op | Hard assertion: `assert actual is not None, ...` |
| **Log-line-only assertion** | Test captures a log line and asserts on its text; a regression that changes the code path but keeps the log still passes | Logs are not the contract | Capture the call kwarg (the flux value, the temperature) and assert on the value passed in |
| **Module-level constant patched only via env var** | `monkeypatch.setenv('RAD_DIR', ...)` on a source that read it at import time | Constants are frozen at import; setenv is too late | `monkeypatch.setattr('mod.CONST', ...)` in addition to setenv |
| **Optional dep imported unconditionally** | `import mors` at module top | `pip install --no-deps` build skips the optional install | `pytest.importorskip('mors')` at module top |
| **Per-band vs per-wavenumber normalization slip** | A test of the surface Planck flux omits the `/1000.0` band-width normalization and still passes because it only checks the sign | The band-width factor is a scale error the test does not probe | Pin the band-integrated flux against `(1 - albedo) sigma ts**4` in the wide-grid limit; the scale guard catches a dropped or doubled normalization |
| **Stale marker after refactor** | File moved into `tests/utils/` without re-applying the module-level `pytestmark` | CI marker filter passed because of per-function markers; coverage tier became invisible | Restore module-level `pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]` |
| **Trivially-true on implicit None** | `def fixture(): pass`; `def test_x(fixture): assert fixture is None` | Fixture returned None implicitly; test passes for the wrong reason | Delete the test |

When you spot a new variant of these, add it here.

---

## 17. Sister rules (cross-link)

- `.github/copilot-instructions.md` "Testing Standards" -- the high-level summary readers without `tests/**` context see first.
- `.github/.claude/rules/janus-code-review.md` "Test marker discipline" -- the review-pass gate that backs up the rules in this file. Also contains domain-aware physics checks (spectral unit boundaries, per-band normalization, dry-adiabat exponent safety, PROTEUS-coupling contract) that apply when reviewing the **source** code that tests cover.

Any change to the rule set: update both files in the same commit and call out the cross-reference in the commit body.
