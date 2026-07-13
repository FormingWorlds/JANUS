# Run and build tests

JANUS uses [pytest](https://docs.pytest.org/) for its test suite. Tests are
sorted into four tiers by cost, verify physical invariants and pinned reference
values, and mirror the source tree. This page covers running the suite and
adding a new test. For the full contract (tiers, markers, coverage gates, the
linter), see the [Testing suite](../Explanations/testing.md) overview.

## Prerequisites

Install the develop extras and make sure `RAD_DIR` and `FWL_DATA` are set:

```bash
pip install -e ".[develop]"
echo $RAD_DIR    # should point to your compiled SOCRATES tree
echo $FWL_DATA   # should point to your data directory
```

Even the mocked unit tests import `janus`, which resolves the SOCRATES
environment at import time, so `RAD_DIR` must point to a built SOCRATES tree for
collection to succeed. See the [installation guide](installation.md) for the
SOCRATES build.

## Running the tests

From the root of the JANUS repository:

```bash
pytest -m "(unit or smoke) and not skip"    # what the PR checks run
pytest -m unit                              # fast unit tests only
pytest -m smoke                             # real SOCRATES, low resolution, one step
pytest -m integration                       # full pipeline coupling (nightly)
pytest -m slow                              # full physics validation (nightly)
pytest -m "not skip"                        # everything that should ever run
```

Run a single file or a single test:

```bash
pytest tests/utils/test_phys.py
pytest tests/modules/test_dry_adiabat_setup.py::test_dry_adiabat_conserves_potential_temperature
```

With coverage:

```bash
pytest --cov=janus --cov-report=term -m "not skip"
pytest --cov=janus --cov-report=html -m "not skip"   # writes htmlcov/
```

## How the tests are organised

Tiers are selected by a module-level marker at the top of every test file:

```python
import pytest

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]
```

| Tier | Runs | Timeout |
|---|---|---|
| `unit` | Python logic, SOCRATES and file I/O mocked | 30 s |
| `smoke` | real SOCRATES binary, one step, low resolution | 60 s |
| `integration` | full `RadConvEqm` / `MCPA_CBL` pipeline | 300 s |
| `slow` | long sweeps and full validation | 3600 s |

New per-source tests mirror the source path: `src/janus/<subdir>/<file>.py` maps
to `tests/<subdir>/test_<file>.py`. The cross-cutting files
(`tests/test_constants.py`, `tests/test_code.py`) and the pipeline files
(`tests/test_runaway_greenhouse.py`, `tests/test_instellation.py`) are the
documented exceptions.

## Adding a test

1. Create `tests/<subdir>/test_<file>.py` and give it the module-level
   `pytestmark` for its tier.
2. Write a file-level docstring naming the source under test, and a function-level
   docstring for every test naming the physical scenario it checks.
3. On a physics source, assert at least one physical invariant (conservation or
   balance, positivity or boundedness, monotonicity or symmetry, or a pinned
   numeric value with a discrimination guard) and tag the test
   `@pytest.mark.physics_invariant`.
4. Give each physics source at least one `@pytest.mark.reference_pinned` test that
   pins against a published benchmark, an analytical limit, or an independent code
   path, and record the anchor on a `docs/Validation/<file>.md` page. The existing
   pages under [Validation anchors](../Validation/phys.md) are the template.
5. Never compare floats with `==`; use `pytest.approx(value, rel=...)` or
   `np.testing.assert_allclose`. On a pinned value, add a follow-up assertion
   showing the most plausible wrong formula would differ by more than the
   tolerance.
6. If the file imports an optional dependency (`hypothesis`, `mors`), call
   `pytest.importorskip('<name>')` at the top of the module before importing it.

## Local checks before opening a PR

```bash
ruff check --fix src/ tests/ && ruff format src/ tests/
bash tools/validate_test_structure.sh          # module-level marker validator
python tools/check_test_quality.py --check      # anti-happy-path linter (blocking)
python tools/check_test_quality.py --reference-pinned-status
```

The same steps run in `.github/workflows/tests.yaml` on every pull request.
