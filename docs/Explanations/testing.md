# Testing suite

[![codecov](https://img.shields.io/codecov/c/github/FormingWorlds/JANUS?label=coverage&logo=codecov&color=brightgreen)](https://app.codecov.io/gh/FormingWorlds/JANUS)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/JANUS/tests.yaml?branch=main&label=Unit%20Tests&color=brightgreen)](https://github.com/FormingWorlds/JANUS/actions/workflows/tests.yaml)
[![Integration Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/JANUS/nightly.yml?branch=main&label=Integration%20Tests&color=brightgreen)](https://github.com/FormingWorlds/JANUS/actions/workflows/nightly.yml)
[![tests](https://img.shields.io/endpoint?url=https://proteus-framework.org/JANUS/badges/tests-total.json)](https://proteus-framework.org/testing)

Tests verify that the code does what was written; physical correctness is judged by data, not by tests.
The suite catches regressions in the saturation-vapour-pressure relations, the Planck emission, the hydrostatic height integration, the heat-capacity fits, the moist and dry adiabats, and the coupled radiative-convective solve, but it cannot certify that those formulae match nature.
That judgement belongs to the validation runs and the published comparisons against measured atmospheric structure.

This page describes the suite as a whole.
Contributors writing or modifying tests should read it together with the [Run and build tests](../How-to/test.md) how-to.

## Test quality contract

Five layers enforce test rigor across the suite:

1. A four-marker tier scheme (`unit`, `smoke`, `integration`, `slow`) selects what runs in the PR gate versus the nightly.
2. Two validation markers (`physics_invariant`, `reference_pinned`) tag tests that carry physical meaning beyond pure code coverage.
3. A per-source mirroring rule pairs every physics source with a same-named test file under `tests/<subdir>/`.
4. An AST linter (`tools/check_test_quality.py`) rejects weak-test patterns on every PR.
5. A coverage ratchet capped at 90 % keeps the gate moving upward over time.

Layers 1, 4, and 5 are blocking on PRs.
Layers 2 and 3 are advisory: the linter reports gaps but does not fail the build.

## The four-marker tier scheme

Every test in the suite carries exactly one tier marker, applied at module level (`pytestmark = [pytest.mark.X, pytest.mark.timeout(N)]`).

| Marker | What it tests | Per-test budget | CI surface |
|---|---|---|---|
| `unit` | Python logic, individual helpers, saturation and heat-capacity fits, surface-flux and adiabat formulae. The SOCRATES binary and file I/O are mocked. | < 100 ms | PR + nightly |
| `smoke` | The real SOCRATES binary on a minimal configuration (low resolution, one step). | < 30 s | PR + nightly |
| `integration` | The full `RadConvEqm` / `MCPA_CBL` pipeline coupling. | minutes | Nightly only |
| `slow` | Long parameter sweeps and full physics validation. | up to an hour | Nightly only |
| `skip` | Placeholder, deliberately disabled. | n/a | Never |

Live counts per tier are shown by the `tests` badge at the top of this page (the total) and by `pytest -m <tier> --collect-only -q` locally.

Tests without a tier marker are invisible to CI.
The PR gate runs `pytest -m "(unit or smoke) and not skip"`; the nightly runs everything in `(unit or smoke or integration or slow) and not skip`.

## Module-level marker and timeout

Every test file declares its tier and its wall-time ceiling at module top:

```python
pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]
```

| Tier | Timeout |
|---|---|
| `unit` | 30 s |
| `smoke` | 60 s |
| `integration` | 300 s |
| `slow` | 3600 s |

The timeout is a defensive ceiling, not a target.
A unit test that takes 25 s of wall time has either picked the wrong tier or has a leak somewhere; the ceiling catches future regressions that introduce a hang.
Per-function markers are additive and do not replace the module-level marker.

## Physics-invariant tiering

A unit test on any physics source must assert at least one of four invariant families:

- **Conservation or balance**: hydrostatic balance along the column; band-integrated surface flux equal to the grey Stefan-Boltzmann limit `(1 - albedo) sigma ts**4`; column array-length and ordering integrity across an adjustment.
- **Positivity or boundedness**: `T > 0`, `P > 0`, saturation vapour pressures non-negative and finite, Planck emission non-negative, `cpv` strictly positive, mole fractions in `[0, 1]`.
- **Monotonicity or symmetry**: saturation vapour pressure increasing with `T`; gravity decreasing with radius; dry-adiabat temperature increasing with pressure; integrated height increasing upward; Planck emission increasing with `T` at fixed frequency.
- **Pinned numeric value with a discrimination guard**: a closed-form value or published table entry pinned via `pytest.approx`, plus an explicit follow-up assertion showing that the most plausible wrong formula would differ from the correct one by more than the tolerance.

Tests that meet one or more of these are tagged `@pytest.mark.physics_invariant`.
The marker is per-function, not module-level: structural tests in the same file (an ordering or pass-through-assignment check) should not carry it.

The six physics sources are:

| Source | Physics |
|---|---|
| `utils/phys.py` | Clausius-Clapeyron saturation vapour pressure; Planck function |
| `utils/height.py` | hydrostatic height integration; inverse-square gravity |
| `utils/cp_funcs.py` | NIST Shomate heat capacity `cpv(vol, tmp)` |
| `modules/spectral_planck_surface.py` | surface Planck emission and its Stefan-Boltzmann limit |
| `modules/moist_adjustment_H2O.py` | moist convective adjustment |
| `modules/dry_adiabat_setup.py` | Poisson dry adiabat / potential temperature |

Utility sources (`__init__.py`, `set_socrates_env.py`, plotting-only modules under `utils/` and `modules/`) are exempt from the invariant requirement but remain subject to the anti-happy-path rules below.

## Reference-pinned validation

Tests that pin behaviour against an external anchor are tagged `@pytest.mark.reference_pinned`.
The anchor is one of:

- a **published benchmark** (cite paper, table, or figure),
- an **analytical limit** (for example, the Stefan-Boltzmann black-body limit or the Poisson potential-temperature identity),
- a **cross-implementation check** (a second, independent code path evaluating the same quantity).

Each of the six physics sources carries at least one reference-pinned test, recorded on the matching `docs/Validation/<file>.md` page.
The current anchors:

| Source | Anchor | Test |
|---|---|---|
| [`utils/phys.py`](../Validation/phys.md) | `satvps(T0) = e0` Clausius-Clapeyron identity; Smithsonian `satvpw` steam point at 1 atm; CODATA fundamental-constant self-consistency; molar masses cross-checked against the gas objects | `tests/utils/test_phys.py::test_satvps_returns_reference_pressure_at_reference_temperature`, `::test_satvpw_matches_one_atm_at_steam_point`, `::test_planck_radiance_positive_and_increases_with_temperature`, `::test_fundamental_constants_are_codata_self_consistent`, `::test_molar_masses_agree_between_lookup_table_and_gas_objects` |
| [`utils/height.py`](../Validation/height.md) | Earth inverse-square gravity; isothermal scale-height integral | `tests/utils/test_height.py::test_surface_gravity_matches_earth_inverse_square`, `::test_integrate_heights_recovers_isothermal_scale_height` |
| [`utils/cp_funcs.py`](../Validation/cp_funcs.md) | NIST-JANAF Shomate water-vapour heat capacity at 500 K | `tests/utils/test_cp_funcs.py::test_water_heat_capacity_matches_nist_shomate_at_500K` |
| [`modules/spectral_planck_surface.py`](../Validation/spectral_planck_surface.md) | grey Stefan-Boltzmann surface limit | `tests/modules/test_spectral_planck_surface.py::test_band_integrated_surface_flux_matches_stefan_boltzmann` |
| [`modules/moist_adjustment_H2O.py`](../Validation/moist_adjustment_H2O.md) | water triple point (273.15 K at 611.657 Pa) and the constant-latent-heat Clausius-Clapeyron inversion | `tests/modules/test_moist_adjustment_H2O.py::test_dew_point_target_matches_clausius_clapeyron_reference` |
| [`modules/dry_adiabat_setup.py`](../Validation/dry_adiabat_setup.md) | Poisson potential-temperature identity | `tests/modules/test_dry_adiabat_setup.py::test_dry_adiabat_conserves_potential_temperature` |
| [`modules/solve_pt.py`](../Validation/solve_pt.md) (pipeline) | pinned runaway-greenhouse OLR and instellation flux references | `tests/test_runaway_greenhouse.py::test_runaway_greenhouse`, `tests/test_instellation.py::test_instellation` |

The marker is not the same thing as physical correctness: a reference-pinned test certifies that *this implementation* reproduces *that anchor*; it does not certify that the anchor is the right physics for every atmospheric regime.

## Per-source-to-test mirroring

New physics tests mirror source at `tests/<subdir>/test_<file>.py`:

| Source | Test |
|---|---|
| `src/janus/utils/phys.py` | `tests/utils/test_phys.py` |
| `src/janus/utils/height.py` | `tests/utils/test_height.py` |
| `src/janus/utils/cp_funcs.py` | `tests/utils/test_cp_funcs.py` |
| `src/janus/modules/spectral_planck_surface.py` | `tests/modules/test_spectral_planck_surface.py` |
| `src/janus/modules/moist_adjustment_H2O.py` | `tests/modules/test_moist_adjustment_H2O.py` |
| `src/janus/modules/dry_adiabat_setup.py` | `tests/modules/test_dry_adiabat_setup.py` |

Cross-cutting and pipeline tests are the documented exception, not the rule:

- `tests/test_constants.py`, `tests/test_code.py`: physical constants and utility / CLI code that span many sources.
- `tests/test_runaway_greenhouse.py`, `tests/test_instellation.py`: the full `RadConvEqm` / `MCPA_CBL` pipeline, which composes the adiabat, adjustment, and radiative steps together and runs the real SOCRATES binary (integration tier).

When a new physics source is added, its per-source test file is created at the same time; the matching `docs/Validation/<file>.md` page is added when the first reference-pinned test for that source lands.

## AST test-quality linter

`tools/check_test_quality.py` walks the test files as an AST and enforces the anti-happy-path rules:

| Rule | What it flags |
|---|---|
| `missing_module_pytestmark` | Test file with no module-level `pytestmark` (a tier marker is required). |
| `missing_docstring` | Test function with no docstring. |
| `single_assert` | Function with exactly one assertion (a single assert rarely discriminates the correct formula from plausible wrong ones). |
| `no_assertions` | Function with zero assertions (only valid for a test that exercises an exception path with `pytest.raises`). |
| `weak_assert_*` | Standalone `assert result is not None`, `assert result > 0`, `assert len(result) > 0`, etc., as the sole meaningful check. A weak assertion *alongside* a strong primary one (the sign guard in a discrimination pattern) is not flagged. |
| `float_eq_literal` | `==` adjacent to a numeric literal in a test body (use `pytest.approx`). |
| `missing_importorskip` | An optional dependency (`hypothesis`, `mors`) imported at module top without a preceding `pytest.importorskip('<name>')`. The PR image installs with `--no-deps`; without `importorskip`, collection fails. |

The linter runs in two modes:

- `python tools/check_test_quality.py --baseline` walks the suite and writes the per-rule violation counts to `tools/test_quality_baseline.json`.
  This is the floor.
  Regenerate the baseline only after a deliberate sweep that reduced violations.
- `python tools/check_test_quality.py --check` (CI mode) walks the suite, compares the current counts to the baseline, and exits non-zero if any rule's count exceeds the baseline.
  The PR workflow runs this and blocks on regression.

An advisory mode reports gaps without failing CI:

- `python tools/check_test_quality.py --reference-pinned-status`: lists physics sources whose matching test file has no `@pytest.mark.reference_pinned` test.

## Local commands

```console
pytest -m unit                              # fast unit tests
pytest -m smoke                             # minimal-config SOCRATES tests
pytest -m integration                       # full pipeline coupling
pytest -m slow                              # sweeps and full validation
pytest -m "(unit or smoke) and not skip"    # PR-gate selection
pytest -m "not skip"                        # everything that should ever run
```

Coverage:

```console
pytest --cov=janus --cov-report=term -m "not skip"
pytest --cov=janus --cov-report=html -m "not skip"   # htmlcov/
```

Lint and structure:

```console
bash tools/validate_test_structure.sh         # module-level marker validator
python tools/check_test_quality.py --check     # AST linter against baseline
python tools/check_test_quality.py --reference-pinned-status
```

## Public-facing badges versus internal taxonomy

Public-facing badges (README, project website) collapse `smoke + integration + slow` into a single `Integration Tests` category, because a four-way taxonomy is confusing to non-developer readers.
The four-marker internal scheme remains for CI granularity: the PR gate runs `(unit or smoke)`, the nightly runs everything, and the test-count badge fetches the JSON files written into the documentation site during the docs deploy.

## Badge system

The documentation deploy (`.github/workflows/docs.yaml`) regenerates three JSON files, `tests-{total,unit,integration}.json`, from `pytest --collect-only` and writes them into the published site under `badges/`.
Shields.io fetches them live from the site via the endpoint URL embedded in the test-count badge.
The counts refresh on every documentation build, so they track the suite without running it.

## Coverage gates

Two gates are declared in `pyproject.toml`:

| Gate | Tests included | Threshold key | Where it runs |
|---|---|---|---|
| Fast | `unit + smoke` | `[tool.janus.coverage_fast].fail_under` | PR `Run unit + smoke tests` step |
| Full | `unit + smoke + integration + slow` | `[tool.coverage.report].fail_under` | Nightly `coverage report` |

Both gates ratchet toward 90 %, capped at 90 % (`tools/update_coverage_threshold.py` enforces `ECOSYSTEM_CEILING = 90.0`); neither may be manually decreased.
The PR gate has a pre-flight step that fetches the base branch's `pyproject.toml` and rejects any PR that drops `[tool.coverage.report].fail_under` below `min(base, 90.0)`.
A one-time ratchet down to the 90 % ceiling is allowed; any drop below 90 % is blocked.

## PR validation pipeline

`.github/workflows/tests.yaml` runs on every PR over `ubuntu-latest` and Python 3.12.
The step sequence:

1. **Build SOCRATES**: checkout and compile the Fortran radiative-transfer dependency (cached across runs).
2. **Build JANUS**: `pip install -e .[develop]`.
3. **Validate test markers** (`bash tools/validate_test_structure.sh`): rejects any test file without a module-level `pytestmark`.
4. **Run test-quality lint** (`python tools/check_test_quality.py --check`): blocking; rejects regression against `tools/test_quality_baseline.json`.
5. **Pre-flight fail_under ratchet check**: rejects any PR that drops `[tool.coverage.report].fail_under` below `min(base, 90.0)`.
6. **Run unit + smoke tests**: `pytest -m "(unit or smoke) and not skip" --cov=janus --cov-fail-under=${FAST_FAIL_UNDER}`, where `FAST_FAIL_UNDER` is read live from `[tool.janus.coverage_fast].fail_under`.

Nightly (`.github/workflows/nightly.yml`) runs the full suite, uploads coverage to Codecov, and enforces the full-gate `fail_under`.

## Canonical specification

The repository-wide rules that every PROTEUS-ecosystem submodule follows are at [proteus-framework.org/PROTEUS/Explanations/ecosystem_testing_standard/](https://proteus-framework.org/PROTEUS/Explanations/ecosystem_testing_standard/).
